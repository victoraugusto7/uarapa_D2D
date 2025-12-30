%==========================================================================
% run_Multicell_Test.m (v3)
%==========================================================================
clear; close all; clc;
tic;

%% 1. Parâmetros de Simulação (CENÁRIO REAL - DISSERTAÇÃO)
fprintf('Definindo parâmetros da simulação (Cenário Real - Dissertação)...\n');

%% 1. Parâmetros de Simulação (CENÁRIO HÍBRIDO: eMBB + URLLC)
fprintf('Definindo parâmetros da simulação (eMBB + URLLC)...\n');

% --- Parâmetros do Cenário ---
params.nCUE       = 34;   % Total de CUEs
params.nD2D       = 9;   % Total de Pares D2D
params.nRBs       = 273*7; % 1911 PRBs (Reuso 3)
params.reuse_factor = 3;

params.eNBradius  = 500; 
params.maxD2DDist = 50;  

% --- Parâmetros de Potência e QoS ---
params.txPowerCUE = 23;   % 23 dBm
params.txPowerD2D = 20;   % 20 dBm
params.gammaCUE   = 10; 
params.gammaD2D   = 10; 

% Interferência Média (Fase 2)
params.Interf_CCI_media_dBm = -110.0; 
params.InterfMargin_dB  = 3.0; 
params.Max_D2D_per_PRB  = 2;

% === PERFIL 1: CUEs (eMBB - Comandantes) ===
params.eMBB_ratio         = 0.5;     % 50% dos CUEs são eMBB
params.eMBB_Req_Rate      = 50e6;    % 50 Mbps
params.eMBB_Req_BER       = 1e-4;    % Confiabilidade Padrão
params.eMBB_Req_Lat       = 50e-3;   % 50 ms

% Modelo de Tráfego eMBB (Vídeo)
params.eMBB_Packet_size   = 1000;    % 1000 bits
params.eMBB_lambda        = 40000;   % 40k pacotes/s (Carga MÉDIA)
params.eMBB_Prop_latency  = 1e-6;    
params.eMBB_Latency_SB    = 30e-3;   % Latência de Backhaul (Server-to-BS)

% === PERFIL 2: D2Ds (URLLC - Agentes) ===
% Nota: D2Ds são todos URLLC neste cenário
params.URLLC_Req_Rate     = 1e6;     % 1 Mbps (Baixa Taxa)
params.URLLC_Req_BER      = 1e-6;    % Alta Confiabilidade (10^-6)
params.URLLC_Req_Lat      = 20e-3;   % 20 ms (Baixa Latência)

% Modelo de Tráfego URLLC (Controle/Telemetria)
params.URLLC_Packet_size  = 1000;    % 1000 bits (Carga = 800kbps)
params.URLLC_lambda       = 800;     % 800 pacotes/s
params.URLLC_Prop_latency = 1e-6;    % Propagação curta (D2D)
params.URLLC_Proc_latency = 1e-3;    % 1 ms (Processamento nos dispositivos)
% (Não há Latency_SB para D2D direto)

% --- PARÂMETROS DAS MHAs ---
params.maxIter    = 5; 
params.colonySize = 20;  
params.limit      = 5;   
params.NLsize     = 10;  

params.colonySize_PC = 20;   
params.maxIter_PC    = 5;   
% -----------------------------------------------------------------
% -----------------------------------------------------------------
%% 2. Parâmetros Físicos
% ... (esta seção permanece exatamente a mesma) ...
phys = struct();
phys.freq           = 3.5e9;   
phys.noiseDensity   = -174;  
phys.plCell.a       = 128.1; phys.plCell.b = 37.6; 
phys.plD2D.a        = 148;   phys.plD2D.b = 40;   
phys.shadowStdCell  = 10;    
phys.shadowStdD2D   = 12;    
phys.gainTx_eNB     = 14;    
phys.gainRx_UE      = 0;     
phys.bwRB           = 360e3; 
params.physParams   = phys;

%% 3. Execução das Etapas 0, 1 e 2
% ETAPA 0: Criar o Cenário Multicelular (7 células)
problem = setup_Multicell_Scenario(params);

% ETAPA 1: Associação de Usuários (UA)
problem = perform_User_Association(problem, params);

% ETAPA 2: Seleção de Candidatos (CPS) - (CHAMADA ATUALIZADA)
V_viabilidade = filter_candidates_cps(problem, params);

toc;
% 5. EXECUÇÃO DA FASE 3 (JOINT RA+PC)
fprintf('\n--- 5. EXECUÇÃO DA FASE 3 (JOINT RA+PC) ---\n');

% 5.1 Parâmetros de PC
fprintf('Parâmetros de PC: Agentes=%d, Iterações=%d\n', params.colonySize_PC, params.maxIter_PC);

% --- Execução dos 3 Algoritmos ---
algorithms = {'ABC', 'OABC', 'GWO'};
results = struct();

for algo_name = algorithms
    name = algo_name{1};
    fprintf('\n--- 5.%s Executando %s ---\n', name(1), name);
    rng('shuffle');

    if strcmp(name, 'ABC')
        % Assumindo que RRA_ABC_Multicell.m foi adaptado
        [bestAssign, bestOF, conv] = RRA_ABC_Multicell(problem, V_viabilidade, params);
    elseif strcmp(name, 'OABC')
        % Assumindo que RRA_OABC_Multicell.m foi adaptado
        [bestAssign, bestOF, conv] = RRA_OABC_Multicell(problem, V_viabilidade, params);
    elseif strcmp(name, 'GWO')
        % Usando a shell RRA_GWO_Multicell.m
        [bestAssign, bestOF, conv] = RRA_GWO_Multicell(problem, V_viabilidade, params);
    end

    results.(name).bestAssign = bestAssign;
    results.(name).bestOF = bestOF;
    results.(name).conv = conv;
end

% --- Teste da Etapa 2 (CPS) ---
fprintf('Etapa 2 (CPS):\n');
nTotalUsers = params.nCUE + params.nD2D;
fprintf('  - Matriz de Viabilidade V(k,r) gerada: %d usuários x %d PRBs\n', ...
    size(V_viabilidade, 1), size(V_viabilidade, 2));

% Contar aprovados
aprovados_cue = sum(any(V_viabilidade(1:params.nCUE, :), 2));
aprovados_d2d = sum(any(V_viabilidade(params.nCUE+1:end, :), 2));
total_aprovados = aprovados_cue + aprovados_d2d;

% Calcular Taxas (%)
taxa_aprov_total = (total_aprovados / nTotalUsers) * 100;
taxa_aprov_cue   = (aprovados_cue / params.nCUE) * 100;
taxa_aprov_d2d   = (aprovados_d2d / params.nD2D) * 100;

fprintf('  - Total de Usuários Aprovados: %d / %d (%.1f%%)\n', ...
    total_aprovados, nTotalUsers, taxa_aprov_total);
fprintf('      -> CUEs Aprovados: %d / %d (%.1f%%)\n', aprovados_cue, params.nCUE, taxa_aprov_cue);
fprintf('      -> D2Ds Aprovados: %d / %d (%.1f%%)\n', aprovados_d2d, params.nD2D, taxa_aprov_d2d);

% --- ANÁLISE DETALHADA DE DESEMPENHO (QoS eMBB + URLLC) ---
fprintf('\n=====================================================================================================\n');
fprintf('                          RELATÓRIO DETALHADO DE QoS                                                 \n');
fprintf('=====================================================================================================\n');

for algo_name = algorithms
    name = algo_name{1};
    res = results.(name);
    
    % 1. Recuperar a Melhor Solução
    X_best = res.bestAssign;
    
    % 2. Recalcular Potências e Métricas Detalhadas
    P_best = optimal_power_control(X_best, problem, params);
    [~, R_vec, BER_vec, L_vec] = calculate_qos_metrics(X_best, P_best, problem, params);
    
    fprintf('\n>> ALGORITMO: %s (Throughput Total: %.2f Mbps)\n', name, res.bestOF/1e6);
    
    % --- RELATÓRIO 1: CUEs eMBB (Comandantes) ---
    idx_eMBB = find(problem.is_eMBB_CUE);
    num_eMBB = length(idx_eMBB);
    served_count = 0;
    
    Req_Rate_e = params.eMBB_Req_Rate / 1e6;
    
    fprintf('   -------------------------------- [ CUEs eMBB - Comandantes ] ----------------------------------\n');
    fprintf('   | ID  | Taxa (Mbps) |    BER     | Lat (ms) | Pwr (dBm) | PRBs |  Req Rate |   STATUS    |\n');
    fprintf('   |-----|-------------|------------|----------|-----------|------|-----------|-------------|\n');
    
    for k = idx_eMBB'
        val_R   = R_vec(k) / 1e6;
        val_BER = BER_vec(k);
        val_Lat = L_vec(k) * 1000;
        val_P   = 10*log10(P_best(k)) + 30;
        num_PRB = sum(X_best(1,:) == k); % Conta PRBs na linha 1 (CUE)

        pass_R = (R_vec(k) >= params.eMBB_Req_Rate);
        pass_B = (val_BER <= params.eMBB_Req_BER);
        pass_L = (val_Lat <= params.eMBB_Req_Lat * 1000);
        
        if pass_R && pass_B && pass_L
            status = 'ATENDIDO';
            served_count = served_count + 1;
        else
            status = 'FALHOU';
            if ~pass_R, status = [status '(Tx)']; end
            if ~pass_B, status = [status '(BER)']; end
            if ~pass_L, status = [status '(Lat)']; end
        end
        
        fprintf('   | %3d |   %7.2f   |  %1.2e  |  %5.2f   |   %5.2f   |  %3d |   > %2.0f    | %-11s |\n', ...
            k, val_R, val_BER, val_Lat, val_P, num_PRB, Req_Rate_e, status);
    end
    fprintf('   => eMBB Atendidos: %.1f%% (%d/%d)\n', (served_count/num_eMBB)*100, served_count, num_eMBB);
    
    % --- RELATÓRIO 2: D2Ds (URLLC - Agentes) ---
    % Os D2Ds estão nos índices globais após os CUEs (nCUE + 1 até nTotal)
    idx_D2D_start = problem.nCUE + 1;
    idx_D2D_end   = problem.nCUE + problem.nD2D;
    
    Req_Rate_u = params.URLLC_Req_Rate / 1e6;
    served_d2d = 0;
    
    fprintf('\n   --------------------------------- [ D2Ds URLLC - Agentes ] ------------------------------------\n');
    fprintf('   | ID  | Taxa (Mbps) |    BER     | Lat (ms) | Pwr (dBm) | PRBs |  Req Rate |   STATUS    |\n');
    fprintf('   |-----|-------------|------------|----------|-----------|------|-----------|-------------|\n');

    for k = idx_D2D_start:idx_D2D_end
        val_R   = R_vec(k) / 1e6;
        val_BER = BER_vec(k);
        val_Lat = L_vec(k) * 1000;
        val_P   = 10*log10(P_best(k)) + 30;
        
        % Conta PRBs nas linhas 2 e 3 da Matriz X
        d2d_prbs = sum(sum(X_best(2:3, :) == k)); 

        pass_R = (R_vec(k) >= params.URLLC_Req_Rate);
        pass_B = (val_BER <= params.URLLC_Req_BER);
        pass_L = (val_Lat <= params.URLLC_Req_Lat * 1000);
        
        if pass_R && pass_B && pass_L
            status = 'ATENDIDO';
            served_d2d = served_d2d + 1;
        else
            status = 'FALHOU';
            if ~pass_R, status = [status '(Tx)']; end
            if ~pass_B, status = [status '(BER)']; end
            if ~pass_L, status = [status '(Lat)']; end
        end
        
        fprintf('   | %3d |   %7.2f   |  %1.2e  |  %5.2f   |   %5.2f   |  %3d |   > %2.0f    | %-11s |\n', ...
            k, val_R, val_BER, val_Lat, val_P, d2d_prbs, Req_Rate_u, status);
    end
    fprintf('   => D2D URLLC Atendidos: %.1f%% (%d/%d)\n', (served_d2d/problem.nD2D)*100, served_d2d, problem.nD2D);
    
    fprintf('   -----------------------------------------------------------------------------------------------------\n');
end
fprintf('\n===========================================================================================\n');
fprintf('\n===========================================================================================\n');
%% 6. VISUALIZAÇÃO FINAL DOS RESULTADOS
fprintf('\n--- Gerando Gráficos Finais... ---\n');

% =========================================================================
% FIGURA 1: CENÁRIO FÍSICO E ASSOCIAÇÃO
% =========================================================================
figure(1); clf;
hold on;
title('Fig 1: Topologia da Rede e Associação (Fase 1)');
xlabel('Posição X (m)'); ylabel('Posição Y (m)'); 
axis equal; grid on;

% 1a. Plotar Células (Hexágonos)
R_cell = problem.eNBradius;
for b = 1:problem.nCells
    plot_hexagon(problem.pos_MBSs(b,1), problem.pos_MBSs(b,2), R_cell);
end

% 1b. Plotar MBSs
h1 = plot(problem.pos_MBSs(:,1), problem.pos_MBSs(:,2), 'sk', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'DisplayName', 'MBS');

% 1c. Separar usuários para legenda correta
idx_cue_c = find(problem.is_eMBB_CUE == 0);
idx_cue_e = find(problem.is_eMBB_CUE == 1);
idx_d2d = 1:problem.nD2D; % Todos D2Ds são comuns no Cenário 1

% 1d. Plotar CUEs
h2 = plot(problem.pos_CUEs(idx_cue_c, 1), problem.pos_CUEs(idx_cue_c, 2), ...
    'ob', 'MarkerSize', 5, 'DisplayName', 'CUE Comum');
h3 = plot(problem.pos_CUEs(idx_cue_e, 1), problem.pos_CUEs(idx_cue_e, 2), ...
    'pb', 'MarkerSize', 9, 'MarkerFaceColor', 'cyan', 'LineWidth', 1.5, ...
    'DisplayName', 'CUE eMBB (Comandante)');

% 1e. Plotar D2Ds (Tx e Rx)
h4 = plot(problem.pos_D2DRx(idx_d2d, 1), problem.pos_D2DRx(idx_d2d, 2), ...
    'or', 'MarkerSize', 5, 'DisplayName', 'D2D-Rx');
h5 = plot(problem.pos_D2DTx(idx_d2d, 1), problem.pos_D2DTx(idx_d2d, 2), ...
    '^r', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', 'D2D-Tx');

legend([h1 h2 h3 h4 h5], 'Location', 'bestoutside');
hold off;

% =========================================================================
% FIGURA 2: BALANCEAMENTO DE CARGA (HISTOGRAMAS)
% =========================================================================
figure(2); clf;
sgtitle('Fig 2: Distribuição de Usuários por Célula (Fase 1: UA)');

subplot(1, 2, 1);
histogram(problem.assoc_CUE, 'BinMethod', 'integer', 'FaceColor', 'b');
title(sprintf('CUEs (Total: %d)', params.nCUE));
xlabel('MBS ID'); ylabel('Quantidade'); grid on;

subplot(1, 2, 2);
histogram(problem.assoc_D2D, 'BinMethod', 'integer', 'FaceColor', 'r');
title(sprintf('D2Ds (Total: %d)', params.nD2D));
xlabel('MBS ID'); ylabel('Quantidade'); grid on;

% =========================================================================
% FIGURA 3: MATRIZ DE VIABILIDADE (FILTRO CPS)
% =========================================================================
figure(3); clf;
imagesc(V_viabilidade);
colormap(flipud(gray)); % Branco = 0 (Reprovado), Preto = 1 (Aprovado)
title(sprintf('Fig 3: Matriz de Viabilidade (Fase 2) - Aprovados: %d/%d', ...
    sum(any(V_viabilidade, 2)), params.nCUE+params.nD2D));
xlabel('Índice do PRB (r)'); 
ylabel('Índice do Usuário (k)');
colorbar;

% =========================================================================
% FIGURA 4: CONVERGÊNCIA (FASE 3)
% =========================================================================
figure(4); clf;
hold on;
if exist('results', 'var')
    l1 = plot(results.ABC.conv / 1e6, 'b-', 'LineWidth', 1.5, 'DisplayName', 'ABC');
    l2 = plot(results.OABC.conv / 1e6, 'm-', 'LineWidth', 1.5, 'DisplayName', 'OABC');
    l3 = plot(results.GWO.conv / 1e6, 'k--', 'LineWidth', 1.5, 'DisplayName', 'GWO');
    
    xlabel('Iteração (Loop Externo)');
    ylabel('Throughput Total do Sistema [Mbps]');
    title('Fig 4: Convergência Comparativa (Joint RA+PC)');
    legend([l1 l2 l3], 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'Dados de Convergência não encontrados.', ...
        'HorizontalAlignment', 'center');
end
hold off;

fprintf('\n--- Visualização Concluída. ---\n');