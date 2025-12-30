
%==========================================================================
% run_Multicell_Test_v4.m
%   Cenário 9 Células - Foco eMBB + D2D URLLC Underlay
%
%   Modelo de banda:
%   - reuse_factor = 3
%   - nRBs = 273 * reuse_factor  => 3 conjuntos de 273 PRBs
%   - Cada MBS enxerga 273 PRBs (~100 MHz) do seu conjunto de reuso.
%==========================================================================

clear; close all; clc;
%addpath(genpath(pwd)); % Garante acesso às subpastas

%% 0. Parâmetros globais de Monte Carlo
nRuns = 1;                    % <-- número de execuções (experimentos)
showPerRunConvergence = false; % se true, plota convergência de cada execução

%% 1. Parâmetros de Simulação (Cenário Real - Reuso 3)
fprintf('Definindo parâmetros da simulação (Cenário 9 Células - Foco eMBB)...\n');

% --- Topologia e Usuários ---
params.nCUE       = 90;   % Total de CUEs
params.nD2D       = 45;   % Total de Pares D2D (Underlay)
params.num_cells  = 9;    % Grid 3x3 (FU-2024)

params.reuse_factor = 3;                         % Reuse-3
params.nRBs         = 273 * params.reuse_factor; % 3 conjuntos de 273 PRBs
% => Cada MBS enxerga exatamente 273 PRBs (~100 MHz)

params.eNBradius  = 500; 
params.maxD2DDist = 50;  

% --- Potência e Interferência ---
params.txPowerCUE = 21;   % 23 dBm
params.txPowerD2D = 12;   % 0-20 dBm
params.gammaCUE   = 10; 
params.gammaD2D   = 10; 
params.InterfMargin_dB       = 3.0; 
params.Max_D2D_per_PRB       = 2;    % Underlay: Até 2 D2Ds por PRB ocupado
params.Interf_CCI_media_dBm  = -110.0; 

% --- PERFIL 1: CUEs eMBB (Prioridade Máxima) ---
params.eMBB_ratio         = 0.5;     % 50% dos CUEs são eMBB
params.eMBB_Req_Rate      = 50e6;    % 50 Mbps
params.eMBB_Req_BER       = 1e-4;    
params.eMBB_Req_Lat       = 50e-3;   % 50 ms

% Tráfego eMBB (Ajustado para Estabilidade)
params.eMBB_Packet_size   = 1000;    
params.eMBB_lambda        = 40000;   % 40 Mbps de carga (para 50 Mbps de serviço)
params.eMBB_Prop_latency  = 1e-6;    
params.eMBB_Latency_SB    = 30e-3;   

% --- PERFIL 2: D2Ds e CUEs Comuns (Best Effort / Básico) ---
params.Common_Req_Rate    = 0.5e6;   % 500 kbps (Voz/Dados básicos)
params.Common_Req_BER     = 1e-3;    
params.Common_Req_Lat     = 100e-3;  % 100 ms
params.Common_lambda      = 400;     % 400 pacotes/s (~400 kbps) - Suportável por 1 PRB

% --- PERFIL 3: D2D URLLC (Agentes Críticos) ---
params.URLLC_Req_Rate     = 1e6;     % 1 Mbps
params.URLLC_Req_BER      = 1e-6;    % 10^-6
params.URLLC_Req_Lat      = 20e-3;   % 20 ms
params.URLLC_Packet_size  = 256*8;   % 256 bytes -> 2048 bits
params.URLLC_lambda       = 250;     % 250 pkts/s -> ~0.5 Mbps

% --- PARÂMETROS DAS MHAs ---
params.maxIter    = 100;  % Iterações do Loop Externo
params.colonySize = 20;  
params.limit      = 5;   
params.NLsize     = 5;  

params.colonySize_PC = 5;   
params.maxIter_PC    = 2;   % GWO Rápido

%% 2. Parâmetros Físicos
phys = struct();
phys.freq           = 3.5e9;   
phys.noiseDensity   = -174;  
phys.plCell.a       = 128.1; phys.plCell.b = 37.6; 
phys.plD2D.a        = 148;   phys.plD2D.b = 40;   
phys.shadowStdCell  = 12;     % Reduzido levemente para estabilidade
phys.shadowStdD2D   = 10;    
phys.gainTx_eNB     = 14;    
phys.gainRx_UE      = 0;     
phys.bwRB           = 360e3; 
params.physParams   = phys;

% Lista de MHAs avaliadas
algorithms = {'ABC', 'OABC', 'WOA'};
nAlgos     = numel(algorithms);

% Pré-alocação dos vetores de desempenho por MHA
convFit_all = struct();
convObj_all = struct();
bestOF_all  = struct();
thput_all   = struct();

for a = 1:nAlgos
    name = algorithms{a};
    convFit_all.(name) = zeros(nRuns, params.maxIter);
    convObj_all.(name) = zeros(nRuns, params.maxIter);
    bestOF_all.(name)  = zeros(nRuns, 1);
    thput_all.(name)   = zeros(nRuns, 1);
end

tic;

%% 3. Loop de Execuções (Monte Carlo)
for runIdx = 1:nRuns
    fprintf('\n============================================================\n');
    fprintf('                    EXECUÇÃO %d / %d\n', runIdx, nRuns);
    fprintf('============================================================\n');

    % ETAPA 0: Setup
    problem = setup_Multicell_Scenario(params);

    % ETAPA 1: Associação (UA)
    problem = perform_User_Association(problem, params);

    % ETAPA 2: Filtro de Candidatos (CPS)
    V_viabilidade = filter_candidates_cps(problem, params);

    % Diagnóstico da Etapa 2
    aprovados_cue = sum(any(V_viabilidade(1:params.nCUE, :), 2));
    aprovados_d2d = sum(any(V_viabilidade(params.nCUE+1:end, :), 2));
    fprintf('Etapa 2: %d/%d CUEs e %d/%d D2Ds aprovados.\n', ...
        aprovados_cue, params.nCUE, aprovados_d2d, params.nD2D);

    %% 5. EXECUÇÃO DA FASE 3 (JOINT RA+PC)
    fprintf('\n--- 5. EXECUÇÃO DA FASE 3 (JOINT RA+PC) ---\n');
    fprintf('Algoritmos: ABC, OABC e WOA (Comparativo)\n');

    for a = 1:nAlgos
        name = algorithms{a};
        fprintf('\n--- Execução %d - Algoritmo %s ---\n', runIdx, name);
        rng('shuffle');

        switch name
            case 'ABC'
                [bestAssign, bestOF, convFit, convObj] = ...
                    RRA_ABC_Multicell(problem, V_viabilidade, params);
            case 'OABC'
                [bestAssign, bestOF, convFit, convObj] = ...
                    RRA_OABC_Multicell(problem, V_viabilidade, params);
            case 'WOA'
                [bestAssign, bestOF, convFit, convObj] = ...
                    RRA_WOA_Multicell(problem, V_viabilidade, params);
        end

        % Garante que vetores sejam linha
        convFit = convFit(:).';
        convObj = convObj(:).';

        % Guarda convergência (truncando se necessário)
        Tfit = min(params.maxIter, numel(convFit));
        Tobj = min(params.maxIter, numel(convObj));
        convFit_all.(name)(runIdx, 1:Tfit) = convFit(1:Tfit);
        convObj_all.(name)(runIdx, 1:Tobj) = convObj(1:Tobj);

        % Métricas finais (throughput físico total)
        X_best = bestAssign;
        P_best = optimal_power_control(X_best, problem, params);
        [~, R_vec, ~, ~] = calculate_qos_metrics(X_best, P_best, problem, params);
        real_throughput_mbps = sum(R_vec) / 1e6;

        bestOF_all.(name)(runIdx) = bestOF;
        thput_all.(name)(runIdx)  = real_throughput_mbps;

        % Guarda resultados da última execução para relatório/gráficos
        if runIdx == nRuns
            results_last.(name).bestAssign = bestAssign;
            results_last.(name).bestOF     = bestOF;
            results_last.(name).convFit    = convFit;
            results_last.(name).convObj    = convObj;
        end

        % Plots individuais de convergência (opcional)
        if showPerRunConvergence
            figure(100 + a); hold on;
            plot(convObj/1e6, 'LineWidth', 1);
            xlabel('Iteração'); ylabel('Throughput (Mbps)');
            title(sprintf('Convergência %s - execução %d', name, runIdx));
            grid on;
        end
    end % loop algoritmos

    % Guarda o cenário da última execução para plots de topologia/viabilidade
    if runIdx == nRuns
        problem_last       = problem;
        V_viabilidade_last = V_viabilidade;
    end

end % loop runs

%% 4. Estatísticas globais (média e desvio padrão por MHA)
fprintf('\n=====================================================================================================\n');
fprintf('                 ESTATÍSTICAS MONTE CARLO (nRuns = %d)\n', nRuns);
fprintf('=====================================================================================================\n');

std_convFit = struct(); % desvio padrão por iteração (fitness)
std_convObj = struct(); % desvio padrão por iteração (objetivo/throughput)

for a = 1:nAlgos
    name = algorithms{a};

    mean_OF = mean(bestOF_all.(name));
    std_OF  = std(bestOF_all.(name));

    mean_TP = mean(thput_all.(name));
    std_TP  = std(thput_all.(name));

    fprintf('\n>> %s:\n', name);
    fprintf('   Fitness final  : média = %10.4f | desvio padrão = %10.4f\n', mean_OF, std_OF);
    fprintf('   Throughput [Mb]: média = %10.4f | desvio padrão = %10.4f\n', mean_TP, std_TP);

    % STD das curvas de convergência (por iteração) – ficam salvos para uso futuro
    std_convFit.(name) = std(convFit_all.(name), 0, 1);
    std_convObj.(name) = std(convObj_all.(name), 0, 1);
end

%% 5. RELATÓRIO DETALHADO DE QoS (apenas se nRuns == 1)
if nRuns == 1 
    fprintf('\n=====================================================================================================\n');
    fprintf('                          RELATÓRIO FINAL DE ALOCAÇÃO E QoS (execução única)\n');
    fprintf('=====================================================================================================\n');

    for a = 1:nAlgos
        name = algorithms{a};
        res  = results_last.(name);

        X_best = res.bestAssign;
        P_best = optimal_power_control(X_best, problem_last, params);
        [~, R_vec, BER_vec, L_vec] = calculate_qos_metrics(X_best, P_best, problem_last, params);

        % --- CÁLCULO DO THROUGHPUT REAL (FÍSICO) ---
        real_throughput_mbps = sum(R_vec) / 1e6;

        fprintf('\n>> ALGORITMO: %s \n', name);
        fprintf('   Fitness (Otimizador):     %10.2f (Inclui penalidades matemáticas)\n', res.bestOF);
        fprintf('   Throughput Real (Físico): %10.2f Mbps (Soma das taxas reais)\n', real_throughput_mbps);

        % --- [1] CUEs eMBB (Prioridade) ---
        idx_eMBB = find(problem_last.is_eMBB_CUE);
        count_ok = 0;

        fprintf('   [1] CUEs eMBB (Requisito: > %.0f Mbps, < %.0f ms)\n', ...
            params.eMBB_Req_Rate/1e6, params.eMBB_Req_Lat*1000);
        fprintf('   | ID  | Taxa (Mbps) |    BER     | Lat (ms) | Pwr (dBm) | PRBs |   STATUS    |\n');
        fprintf('   |-----|-------------|------------|----------|-----------|------|-------------|\n');

        for k = idx_eMBB.'
            val_R = R_vec(k)/1e6; 
            val_L = L_vec(k)*1000; 
            val_P = 10*log10(P_best(k))+30;
            nPRB  = sum(X_best(1,:) == k);

            ok = (R_vec(k) >= params.eMBB_Req_Rate) && ...
                 (L_vec(k) <= params.eMBB_Req_Lat)   && ...
                 (BER_vec(k) <= params.eMBB_Req_BER);

            if ok
                status='ATENDIDO'; 
                count_ok = count_ok+1; 
            else
                status='FALHOU'; 
                if R_vec(k) < params.eMBB_Req_Rate, status=[status '(Tx)']; end
                if L_vec(k) > params.eMBB_Req_Lat,  status=[status '(Lt)']; end
                if BER_vec(k) > params.eMBB_Req_BER, status=[status '(Be)']; end
            end

            fprintf('   | %3d |   %7.2f   |  %1.2e  |  %5.2f   |   %5.2f   |  %3d | %-11s |\n', ...
                k, val_R, BER_vec(k), val_L, val_P, nPRB, status);
        end
        fprintf('   => Sucesso eMBB: %.1f%% (%d/%d)\n', ...
            (count_ok/length(idx_eMBB))*100, count_ok, length(idx_eMBB));

        % --- [2] CUEs COMUNS ---
        idx_Common = find(~problem_last.is_eMBB_CUE & (1:params.nCUE).');
        valid_common = [];
        for k = idx_Common.'
            if problem_last.assoc_CUE(k) > 0 && any(V_viabilidade_last(k,:))
                valid_common = [valid_common; k]; 
            end
        end

        count_active = 0;
        if ~isempty(valid_common)
            fprintf('\n   [2] CUEs Comuns Admitidos (Melhor Esforço)\n');
            fprintf('   | ID  | Taxa (Mbps) |    BER     | Lat (ms) | Pwr (dBm) | PRBs |   STATUS    |\n');
            fprintf('   |-----|-------------|------------|----------|-----------|------|-------------|\n');
            for k = valid_common.'
                val_R = R_vec(k)/1e6; 
                val_L = L_vec(k)*1000; 
                val_P = 10*log10(P_best(k))+30;
                nPRB  = sum(X_best(1,:) == k);

                if nPRB > 0
                    status = 'ATIVO'; 
                    count_active = count_active + 1;
                else
                    status = 'SEM PRB'; 
                end

                fprintf('   | %3d |   %7.2f   |  %1.2e  |  %5.2f   |   %5.2f   |  %3d | %-11s |\n', ...
                    k, val_R, BER_vec(k), val_L, val_P, nPRB, status);
            end
            fprintf('   => CUEs Comuns Ativos: %d/%d admitidos\n', count_active, length(valid_common));
        end

        % --- [3] D2D PAIRS (URLLC / Underlay) ---
        idx_D2D = (params.nCUE + 1 : params.nCUE + params.nD2D).';
        valid_d2d = [];
        for k = idx_D2D.'
            d_idx = k - params.nCUE;
            if problem_last.assoc_D2D(d_idx) > 0 && any(V_viabilidade_last(k,:))
                valid_d2d = [valid_d2d; k]; 
            end
        end

        count_d2d_active = 0;
        if ~isempty(valid_d2d)
            fprintf('\n   [3] Pares D2D Admitidos (Underlay)\n');
            fprintf('   | ID  | Taxa (Mbps) |    BER     | Lat (ms) | Pwr (dBm) | PRBs (Reuso) |   STATUS    |\n');
            fprintf('   |-----|-------------|------------|----------|-----------|--------------|-------------|\n');
            for k = valid_d2d.'
                val_R = R_vec(k)/1e6; 
                val_P = 10*log10(P_best(k))+30;
                val_L = L_vec(k)*1000; 
                nPRB  = sum(sum(X_best(2:3,:) == k));

                if nPRB > 0
                    status = 'REUSANDO';
                    count_d2d_active = count_d2d_active + 1;
                else
                    status = 'SEM PRB'; 
                end

                fprintf('   | %3d |   %7.2f   |  %1.2e  |  %5.2f   |   %5.2f   |      %3d     | %-11s |\n', ...
                    k, val_R, BER_vec(k), val_L, val_P, nPRB, status);
            end
            fprintf('   => D2Ds Ativos: %d/%d admitidos\n', count_d2d_active, length(valid_d2d));
        end
        fprintf('-----------------------------------------------------------------------------------------------------\n');
    end
end

%% 6. VISUALIZAÇÃO FINAL (Gráficos)
fprintf('\n--- Gerando Gráficos... ---\n');

% FIG 1: Topologia 9 Células (Grid Hexagonal + Usuários) da última execução
figure(1); clf; hold on;
title('Fig 1: Topologia 9 Células (Honeycomb) e Usuários (última execução)');
xlabel('X (m)'); ylabel('Y (m)'); axis equal; grid on;

% 1. Plota Células
R = problem_last.eNBradius;
for b = 1:problem_last.nCells
    plot_hexagon(problem_last.pos_MBSs(b,1), problem_last.pos_MBSs(b,2), R);
end

% 2. Plota MBSs
h_mbs = plot(problem_last.pos_MBSs(:,1), problem_last.pos_MBSs(:,2), 'sk', ...
    'MarkerFaceColor','k', 'MarkerSize',8, 'DisplayName','MBS');

% 3. Plota Usuários
idx_e = find(problem_last.is_eMBB_CUE); 
idx_c = find(~problem_last.is_eMBB_CUE);

if ~isempty(idx_c)
    h_cue_c = plot(problem_last.pos_CUEs(idx_c,1), problem_last.pos_CUEs(idx_c,2), ...
        'ob', 'MarkerSize',5, 'DisplayName','CUE Comum');
else
    h_cue_c = [];
end

if ~isempty(idx_e)
    h_cue_e = plot(problem_last.pos_CUEs(idx_e,1), problem_last.pos_CUEs(idx_e,2), ...
        'p', 'MarkerSize',11, 'MarkerFaceColor','c', 'MarkerEdgeColor','b', ...
        'DisplayName','CUE eMBB (Comandante)');
else
    h_cue_e = [];
end

% 4. Plota D2Ds Tx/Rx
h_d2d_tx = plot(problem_last.pos_D2DTx(:,1), problem_last.pos_D2DTx(:,2), ...
    '^r', 'MarkerSize',6, 'MarkerFaceColor','r', 'DisplayName','D2D Tx (Agente)');
h_d2d_rx = plot(problem_last.pos_D2DRx(:,1), problem_last.pos_D2DRx(:,2), ...
    'or', 'MarkerSize',6, 'DisplayName','D2D Rx');

legend([h_mbs, h_cue_c, h_cue_e, h_d2d_tx, h_d2d_rx], 'Location','bestoutside');
hold off;

% FIG 3: Matriz de Viabilidade da última execução
figure(3); clf;
imagesc(V_viabilidade_last);
colormap(flipud(gray)); 
title(sprintf('Fig 3: Matriz de Viabilidade (Fase 2) - Aprovados: %d/%d (última execução)', ...
    sum(any(V_viabilidade_last, 2)), params.nCUE+params.nD2D));
xlabel('PRB Index'); ylabel('User Index (CUEs... D2Ds...)');

% FIG 4: Convergência MÉDIA (Função Utilidade / Fitness)
figure(4); clf; hold on;
for a = 1:nAlgos
    name = algorithms{a};
    mean_convFit = mean(convFit_all.(name), 1);
    plot(mean_convFit, 'LineWidth',1.5, 'DisplayName', name);
end
xlabel('Iteração'); ylabel('Fitness (Com Penalidades)');
title(sprintf('Fig 4: Convergência MÉDIA do Fitness (%d execuções)', nRuns));
legend('show', 'Location', 'best'); grid on;

% FIG 5: Convergência MÉDIA (Função Objetivo / Throughput Real)
figure(5); clf; hold on;
for a = 1:nAlgos
    name = algorithms{a};
    mean_convObj = mean(convObj_all.(name), 1)/1e6;
    plot(mean_convObj, 'LineWidth',1.5, 'DisplayName', name);
end
xlabel('Iteração'); ylabel('Throughput Total MÉDIO (Mbps)');
title(sprintf('Fig 5: Função Objetivo - Throughput REAL MÉDIO do Sistema (%d execuções)', nRuns));
legend('show', 'Location', 'best'); grid on;

fprintf('Visualização Concluída.\n');
toc;
