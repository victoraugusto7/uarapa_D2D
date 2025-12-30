%==========================================================================
% run_Multicell_Sweep.m
%   Varre cenário multicell alterando nCUE, n_eMBB e nD2D e gera gráficos:
%   1) Throughput total x nº de CUEs
%   2) CUEs eMBB atendidos x nº de CUEs eMBB
%   3) Throughput total x nº de pares D2D
%   4) D2D URLLC atendidos x nº de pares D2D
%
%   Modelo de banda:
%   - reuse_factor = 3
%   - nRBs = 273 * reuse_factor  => 3 conjuntos de 273 PRBs
%   - Cada MBS enxerga 273 PRBs (~100 MHz) do seu conjunto de reuso.
%
%   *** Versão média/STD ***
%   - Executa nRuns vezes cada ponto da varredura
%   - Plota as CURVAS MÉDIAS dos algoritmos
%   - Guarda também os DESVIOS-PADRÃO para análise estatística
%==========================================================================

clear; close all; clc;
addpath(genpath(pwd)); % Garante acesso às subpastas

%% 0. Configuração de número de execuções
nRuns = 32;     % <-- você pode ajustar aqui o nº de execuções por ponto

fprintf('Definindo parâmetros BASE da simulação...\n');

base_params = struct();

% --- Topologia e Usuários (BASE) ---
base_params.nCUE       = 70;     % nº de CUEs base
base_params.nD2D       = 35;     % nº de pares D2D base
base_params.num_cells  = 9;      % Grid 3x3 (FU-2024)

base_params.reuse_factor = 3;    % Reuse-3
base_params.nRBs         = 273 * base_params.reuse_factor;
% => 3 conjuntos de 273 PRBs; cada MBS usa exatamente 273 PRBs

base_params.eNBradius  = 500;
base_params.maxD2DDist = 50;

% --- Potência e Interferência ---
base_params.txPowerCUE = 23;   % 23 dBm
base_params.txPowerD2D = 20;   % 0-20 dBm
base_params.gammaCUE   = 10;
base_params.gammaD2D   = 10;
base_params.InterfMargin_dB      = 3.0;
base_params.Max_D2D_per_PRB      = 2;
base_params.Interf_CCI_media_dBm = -110.0;

% --- PERFIL 1: CUEs eMBB (Prioridade Máxima) ---
base_params.eMBB_ratio         = 0.5;     % 50% dos CUEs são eMBB
base_params.eMBB_Req_Rate      = 50e6;    % 50 Mbps
base_params.eMBB_Req_BER       = 1e-4;
base_params.eMBB_Req_Lat       = 50e-3;   % 50 ms

% Tráfego eMBB
base_params.eMBB_Packet_size   = 1000;
base_params.eMBB_lambda        = 40000;   % 40 Mbps de carga
base_params.eMBB_Prop_latency  = 1e-6;
base_params.eMBB_Latency_SB    = 30e-3;   % 30 ms de backhaul

% --- PERFIL 2: D2Ds e CUEs Comuns (Best Effort) ---
base_params.Common_Req_Rate    = 0.5e6;   % 500 kbps
base_params.Common_Req_BER     = 1e-3;
base_params.Common_Req_Lat     = 100e-3;  % 100 ms
base_params.Common_lambda      = 400;     % ~400 kbps

% --- PERFIL 3: D2D URLLC (Agentes Críticos) ---
base_params.URLLC_Req_Rate     = 1e6;     % 1 Mbps
base_params.URLLC_Req_BER      = 1e-6;    % 10^-6
base_params.URLLC_Req_Lat      = 20e-3;   % 20 ms
base_params.URLLC_Packet_size  = 256*8;   % 256 bytes -> 2048 bits
base_params.URLLC_lambda       = 250;     % 250 pkts/s -> ~0.5 Mbps

% --- PARÂMETROS DAS MHAs ---
base_params.maxIter    = 15;  % Iterações do Loop Externo
base_params.colonySize = 15;
base_params.limit      = 5;
base_params.NLsize     = 5;

base_params.colonySize_PC = 10;
base_params.maxIter_PC    = 4;   % GWO rápido

% --- Parâmetros Físicos ---
phys = struct();
phys.freq           = 3.5e9;
phys.noiseDensity   = -174;
phys.plCell.a       = 128.1; phys.plCell.b = 37.6;
phys.plD2D.a        = 148;   phys.plD2D.b = 40;
phys.shadowStdCell  = 12;
phys.shadowStdD2D   = 10;
phys.gainTx_eNB     = 14;
phys.gainRx_UE      = 0;
phys.bwRB           = 360e3;

base_params.physParams = phys;

% Algoritmos a comparar
algorithms = {'ABC', 'OABC', 'WOA'};
nAlgo = numel(algorithms);

%% 1. Vetores de Varredura
vec_nCUE       = [40 50 60 70 80 90];        % 1) nº total de CUEs
vec_eMBB_ratio = [0.2 0.3 0.4 0.5 0.6 0.7]; % 2) fração eMBB
vec_nD2D       = [10 20 30 40 50 60];       % 3) nº de pares D2D

NnCUE   = numel(vec_nCUE);
NeMBB   = numel(vec_eMBB_ratio);
NnD2D   = numel(vec_nD2D);

%% Pré-alocação dos resultados (MÉDIAS e STD)

% 1) Throughput x nº CUE
thr_runs_nCUE     = zeros(nAlgo, NnCUE, nRuns); % bruto (para média/STD)

% 2) eMBB atendidos x nº eMBB
eMBB_total_runs   = zeros(NeMBB, nRuns);        % nº eMBB sorteados
eMBB_served_runs  = zeros(nAlgo, NeMBB, nRuns);

% 3) Throughput x nº D2D
thr_runs_nD2D     = zeros(nAlgo, NnD2D, nRuns);

% 4) URLLC atendidos x nº D2D
URLLC_total_runs  = zeros(NnD2D, nRuns);
URLLC_served_runs = zeros(nAlgo, NnD2D, nRuns);

%% 2. VARREDURA 1: Throughput x nº de CUEs
fprintf('\n=== VARREDURA 1: Throughput x nº de CUEs (média de %d execuções) ===\n', nRuns);
tic;
for idx = 1:NnCUE
    for r = 1:nRuns
        params      = base_params;
        params.nCUE = vec_nCUE(idx);
        % nD2D fixo = base_params.nD2D

        fprintf('  -> Run %2d/%2d | nCUE = %d, nD2D = %d...\n', ...
            r, nRuns, params.nCUE, params.nD2D);

        summary = run_multicell_scenario(params, algorithms);

        for a = 1:nAlgo
            thr_runs_nCUE(a, idx, r) = summary.throughput_Mbps(a);
        end
    end
end
toc;

%% 3. VARREDURA 2: nº eMBB atendidos x nº de eMBB
fprintf('\n=== VARREDURA 2: eMBB atendidos x nº de eMBB (média de %d execuções) ===\n', nRuns);
tic;
for idx = 1:NeMBB
    for r = 1:nRuns
        params             = base_params;
        params.nCUE        = base_params.nCUE;    % fixa nº CUEs
        params.eMBB_ratio  = vec_eMBB_ratio(idx); % varia fração de eMBB

        fprintf('  -> Run %2d/%2d | nCUE = %d, eMBB_ratio = %.2f...\n', ...
            r, nRuns, params.nCUE, params.eMBB_ratio);

        summary = run_multicell_scenario(params, algorithms);

        eMBB_total_runs(idx, r) = summary.eMBB_total;

        for a = 1:nAlgo
            eMBB_served_runs(a, idx, r) = summary.eMBB_served(a);
        end
    end
end
toc;

%% 4. VARREDURA 3: Throughput x nº de D2D
fprintf('\n=== VARREDURA 3: Throughput x nº de Pares D2D (média de %d execuções) ===\n', nRuns);
tic;
for idx = 1:NnD2D
    for r = 1:nRuns
        params      = base_params;
        params.nCUE = base_params.nCUE;      % fixa nº CUEs
        params.nD2D = vec_nD2D(idx);        % varia nº D2D

        fprintf('  -> Run %2d/%2d | nCUE = %d, nD2D = %d...\n', ...
            r, nRuns, params.nCUE, params.nD2D);

        summary = run_multicell_scenario(params, algorithms);

        URLLC_total_runs(idx, r) = summary.URLLC_total;

        for a = 1:nAlgo
            thr_runs_nD2D(a, idx, r)     = summary.throughput_Mbps(a);
            URLLC_served_runs(a, idx, r) = summary.URLLC_served(a);
        end
    end
end
toc;

%% 5. Cálculo das MÉDIAS e DESVIOS-PADRÃO

% --- (1) Throughput Total x nº de CUEs ---
thr_vs_nCUE      = mean(thr_runs_nCUE, 3);               % média
thr_std_vs_nCUE  = std( thr_runs_nCUE, 0, 3);            % std

% --- (2) CUEs eMBB atendidos x nº de CUEs eMBB ---
n_eMBB_total_vec = mean(eMBB_total_runs, 2).';           % média no eixo x
n_eMBB_satis_vs_nCUEe     = mean(eMBB_served_runs, 3);   % média atendidos
n_eMBB_satis_std_vs_nCUEe = std( eMBB_served_runs, 0, 3);

% --- (3) Throughput Total x nº de Pares D2D ---
thr_vs_nD2D      = mean(thr_runs_nD2D, 3);
thr_std_vs_nD2D  = std( thr_runs_nD2D, 0, 3);

% --- (4) D2D URLLC atendidos x nº de Pares D2D ---
n_URLLC_total_vec         = mean(URLLC_total_runs, 2).';
n_URLLC_satis_vs_nD2D     = mean(URLLC_served_runs, 3);
n_URLLC_satis_std_vs_nD2D = std( URLLC_served_runs, 0, 3);

%% 6. PLOTS (usando as MÉDIAS)
% Se quiser, depois você pode trocar 'plot' por 'errorbar' usando *_std.

% --- (1) Throughput Total x nº de CUEs ---
figure; hold on; grid on;
for a = 1:nAlgo
    plot(vec_nCUE, thr_vs_nCUE(a,:), '-o', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
    % Exemplo com barras de erro (comente/descomente se quiser):
    % errorbar(vec_nCUE, thr_vs_nCUE(a,:), thr_std_vs_nCUE(a,:), ...
    %     '-o','LineWidth',1.0,'HandleVisibility','off');
end
xlabel('Número de CUEs');
ylabel('Throughput Total Médio (Mbps)');
title(sprintf('Throughput Total x Número de CUEs (média de %d execuções)', nRuns));
legend('Location','best');
hold off;

% --- (2) CUEs eMBB atendidos x nº de CUEs eMBB ---
figure; hold on; grid on;
x_eMBB = n_eMBB_total_vec; % eixo x em nº de eMBB (médio)

for a = 1:nAlgo
    plot(x_eMBB, n_eMBB_satis_vs_nCUEe(a,:), '-s', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
    % errorbar(x_eMBB, n_eMBB_satis_vs_nCUEe(a,:), ...
    %          n_eMBB_satis_std_vs_nCUEe(a,:), '-s', ...
    %          'LineWidth',1.0,'HandleVisibility','off');
end
xlabel('Número de CUEs eMBB');
ylabel('Número Médio de CUEs eMBB Atendidos');
title(sprintf('Atendimento de QoS eMBB x Número de CUEs eMBB (média de %d execuções)', nRuns));
legend('Location','best');
hold off;

% --- (3) Throughput Total x nº de Pares D2D ---
figure; hold on; grid on;
for a = 1:nAlgo
    plot(vec_nD2D, thr_vs_nD2D(a,:), '-^', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
    % errorbar(vec_nD2D, thr_vs_nD2D(a,:), thr_std_vs_nD2D(a,:), ...
    %          '-^','LineWidth',1.0,'HandleVisibility','off');
end
xlabel('Número de Pares D2D');
ylabel('Throughput Total Médio (Mbps)');
title(sprintf('Throughput Total x Número de Pares D2D (média de %d execuções)', nRuns));
legend('Location','best');
hold off;

% --- (4) D2D URLLC atendidos x nº de Pares D2D ---
figure; hold on; grid on;
x_URLLC = n_URLLC_total_vec; % normalmente igual a vec_nD2D (médio)

for a = 1:nAlgo
    plot(x_URLLC, n_URLLC_satis_vs_nD2D(a,:), '-d', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
    % errorbar(x_URLLC, n_URLLC_satis_vs_nD2D(a,:), ...
    %          n_URLLC_satis_std_vs_nD2D(a,:), '-d', ...
    %          'LineWidth',1.0,'HandleVisibility','off');
end
xlabel('Número de Pares D2D (URLLC)');
ylabel('Número Médio de D2D URLLC Atendidos');
title(sprintf('Atendimento URLLC x Número de Pares D2D (média de %d execuções)', nRuns));
legend('Location','best');
hold off;

fprintf('\nVarreduras concluídas. Curvas médias e desvios-padrão calculados.\n');

%% ===================================================================== %%
%% FUNÇÃO AUXILIAR: roda UM cenário multicell completo e devolve métrica %%
%% ===================================================================== %%

function summary = run_multicell_scenario(params, algorithms)
% Roda setup + UA + CPS + RA+PC para um conjunto de parâmetros.
% Retorna métricas agregadas por algoritmo:
%   - throughput_Mbps(a)
%   - eMBB_total, eMBB_served(a)
%   - URLLC_total, URLLC_served(a)

    % ETAPA 0: Setup
    problem = setup_Multicell_Scenario(params);

    % ETAPA 1: Associação (UA)
    problem = perform_User_Association(problem, params);

    % ETAPA 2: Filtro de Candidatos (CPS)
    V_viabilidade = filter_candidates_cps(problem, params);

    nAlgo = numel(algorithms);
    summary.throughput_Mbps = zeros(1, nAlgo);
    summary.eMBB_served     = zeros(1, nAlgo);
    summary.URLLC_served    = zeros(1, nAlgo);

    summary.eMBB_total  = sum(problem.is_eMBB_CUE);
    if isfield(problem,'is_URLLC_D2D')
        summary.URLLC_total = sum(problem.is_URLLC_D2D);
    else
        summary.URLLC_total = 0;
    end

    % LOOP sobre algoritmos de RA
    for a = 1:nAlgo
        name = algorithms{a};
        fprintf('     -> Algoritmo: %s\n', name);

        rng('shuffle'); % semente aleatória para cada exec

        switch name
            case 'ABC'
                [bestAssign, ~, ~, ~] = ...
                    RRA_ABC_Multicell(problem, V_viabilidade, params);
            case 'OABC'
                [bestAssign, ~, ~, ~] = ...
                    RRA_OABC_Multicell(problem, V_viabilidade, params);
            case 'WOA'
                [bestAssign, ~, ~, ~] = ...
                    RRA_WOA_Multicell(problem, V_viabilidade, params);
            otherwise
                error('Algoritmo desconhecido: %s', name);
        end

        % Otimiza potência para a melhor alocação encontrada
        P_best = optimal_power_control(bestAssign, problem, params);

        % Calcula QoS
        [R_total, R_vec, BER_vec, L_vec] = ...
            calculate_qos_metrics(bestAssign, P_best, problem, params);

        % Throughput físico (Mbps)
        summary.throughput_Mbps(a) = R_total / 1e6;

        % --- Contagem de CUEs eMBB atendidos ---
        idx_eMBB = find(problem.is_eMBB_CUE);
        count_ok_eMBB = 0;
        for k = idx_eMBB'
            if (R_vec(k) >= problem.eMBB_Req_Rate) && ...
               (BER_vec(k) <= problem.eMBB_Req_BER)   && ...
               (L_vec(k)  <= problem.eMBB_Req_Lat)
                count_ok_eMBB = count_ok_eMBB + 1;
            end
        end
        summary.eMBB_served(a) = count_ok_eMBB;

        % --- Contagem de D2D URLLC atendidos ---
        if summary.URLLC_total > 0 && isfield(problem,'is_URLLC_D2D')
            idx_URLLC = find(problem.is_URLLC_D2D);  % índices locais (1..nD2D)
            count_ok_URLLC = 0;
            for d = idx_URLLC'
                kD2D = problem.nCUE + d;  % índice global do user D2D
                if (R_vec(kD2D) >= problem.URLLC_Req_Rate) && ...
                   (BER_vec(kD2D) <= problem.URLLC_Req_BER)   && ...
                   (L_vec(kD2D)  <= problem.URLLC_Req_Lat)
                    count_ok_URLLC = count_ok_URLLC + 1;
                end
            end
            summary.URLLC_served(a) = count_ok_URLLC;
        end
    end
end
