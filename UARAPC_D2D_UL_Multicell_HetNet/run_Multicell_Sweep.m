%==========================================================================
% run_Multicell_Sweep.m
%   Varre cenário multicell alterando nCUE, n_eMBB e nD2D e gera gráficos:
%   1) Throughput total x nº de CUEs
%   2) CUEs eMBB atendidos x nº de CUEs eMBB
%   3) Throughput total x nº de pares D2D
%   4) D2D URLLC atendidos x nº de pares D2D
%==========================================================================

clear; close all; clc;
addpath(genpath(pwd)); % Garante acesso às subpastas

%% 0. Parâmetros BASE de simulação (ajuste se quiser)
fprintf('Definindo parâmetros BASE da simulação...\n');

base_params = struct();

% --- Topologia e Usuários (BASE) ---
base_params.nCUE       = 70;     % nº de CUEs base
base_params.nD2D       = 35;     % nº de pares D2D base
base_params.num_cells  = 9;      % Grid 3x3 (FU-2024)
base_params.nRBs       = 273*7;  % 1911 PRBs (pool global para reuso 3)
base_params.reuse_factor = 3;

base_params.eNBradius  = 500; 
base_params.maxD2DDist = 50;  

% --- Potência e Interferência ---
base_params.txPowerCUE = 23;   % 23 dBm
base_params.txPowerD2D = 19;   % 19 dBm
base_params.gammaCUE   = 10; 
base_params.gammaD2D   = 10; 
base_params.InterfMargin_dB  = 3.0; 
base_params.Max_D2D_per_PRB  = 2;    
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
base_params.maxIter    = 10;  % Iterações do Loop Externo
base_params.colonySize = 20;  
base_params.limit      = 5;   
base_params.NLsize     = 10;  

base_params.colonySize_PC = 20;   
base_params.maxIter_PC    = 5;   % GWO rápido

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

%% 1. Vetores de Varredura (ajuste conforme desejado)
% 1) Variação do nº total de CUEs
vec_nCUE = [40 50 60 70 80 90];   % exemplo

% 2) Variação do nº de CUEs eMBB (via razão eMBB_ratio)
%    Mantemos nCUE fixo (base_params.nCUE) e variamos a fração de eMBB
vec_eMBB_ratio = [0.2 0.3 0.4 0.5 0.6 0.7];  % => nº eMBB = ratio * nCUE

% 3) Variação do nº de Pares D2D
vec_nD2D = [10 20 30 40 50 60];  % exemplo

% Pré-alocação dos resultados
thr_vs_nCUE        = zeros(nAlgo, numel(vec_nCUE));        % (1)
n_eMBB_satis_vs_nCUEe = zeros(nAlgo, numel(vec_eMBB_ratio)); % (2)
n_eMBB_total_vec   = zeros(1, numel(vec_eMBB_ratio));      % eixo x real

thr_vs_nD2D        = zeros(nAlgo, numel(vec_nD2D));        % (3)
n_URLLC_satis_vs_nD2D = zeros(nAlgo, numel(vec_nD2D));     % (4)
n_URLLC_total_vec  = zeros(1, numel(vec_nD2D));            % eixo x real

%% 2. VARREDURA 1: Throughput x nº de CUEs
fprintf('\n=== VARREDURA 1: Throughput x nº de CUEs ===\n');
tic;
for idx = 1:numel(vec_nCUE)
    params = base_params;
    params.nCUE = vec_nCUE(idx);
    
    fprintf('  -> Simulando cenário com nCUE = %d, nD2D = %d...\n', ...
        params.nCUE, params.nD2D);
    
    summary = run_multicell_scenario(params, algorithms);
    
    for a = 1:nAlgo
        thr_vs_nCUE(a, idx) = summary.throughput_Mbps(a);
    end
end

%% 3. VARREDURA 2: nº eMBB atendidos x nº de eMBB
fprintf('\n=== VARREDURA 2: eMBB atendidos x nº de eMBB ===\n');

for idx = 1:numel(vec_eMBB_ratio)
    params = base_params;
    params.nCUE        = base_params.nCUE;      % fixa nº CUEs
    params.eMBB_ratio  = vec_eMBB_ratio(idx);   % varia fração de eMBB
    
    fprintf('  -> Simulando cenário com nCUE = %d, eMBB_ratio = %.2f...\n', ...
        params.nCUE, params.eMBB_ratio);
    
    summary = run_multicell_scenario(params, algorithms);
    
    n_eMBB_total_vec(idx) = summary.eMBB_total;  % nº de eMBB realmente sorteados
    
    for a = 1:nAlgo
        n_eMBB_satis_vs_nCUEe(a, idx) = summary.eMBB_served(a);
    end
end

%% 4. VARREDURA 3: Throughput x nº de D2D
fprintf('\n=== VARREDURA 3: Throughput x nº de Pares D2D ===\n');

for idx = 1:numel(vec_nD2D)
    params = base_params;
    params.nCUE = base_params.nCUE;   % fixa nº CUEs
    params.nD2D = vec_nD2D(idx);      % varia nº D2D
    
    fprintf('  -> Simulando cenário com nCUE = %d, nD2D = %d...\n', ...
        params.nCUE, params.nD2D);
    
    summary = run_multicell_scenario(params, algorithms);
    
    n_URLLC_total_vec(idx) = summary.URLLC_total;  % normalmente = nD2D
    
    for a = 1:nAlgo
        thr_vs_nD2D(a, idx) = summary.throughput_Mbps(a);
        n_URLLC_satis_vs_nD2D(a, idx) = summary.URLLC_served(a);
    end
end

%% 5. PLOTS

% --- (1) Throughput Total x nº de CUEs ---
figure; hold on; grid on;
for a = 1:nAlgo
    plot(vec_nCUE, thr_vs_nCUE(a,:), '-o', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
end
xlabel('Número de CUEs');
ylabel('Throughput Total (Mbps)');
title('Throughput Total x Número de CUEs');
legend('Location','best');
hold off;

% --- (2) CUEs eMBB atendidos x nº de CUEs eMBB ---
figure; hold on; grid on;
x_eMBB = n_eMBB_total_vec; % eixo x em nº de eMBB

for a = 1:nAlgo
    plot(x_eMBB, n_eMBB_satis_vs_nCUEe(a,:), '-s', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
end
xlabel('Número de CUEs eMBB');
ylabel('Número de CUEs eMBB Atendidos');
title('Atendimento de QoS eMBB x Número de CUEs eMBB');
legend('Location','best');
hold off;

% --- (3) Throughput Total x nº de Pares D2D ---
figure; hold on; grid on;
for a = 1:nAlgo
    plot(vec_nD2D, thr_vs_nD2D(a,:), '-^', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
end
xlabel('Número de Pares D2D');
ylabel('Throughput Total (Mbps)');
title('Throughput Total x Número de Pares D2D');
legend('Location','best');
hold off;

% --- (4) D2D URLLC atendidos x nº de Pares D2D ---
figure; hold on; grid on;
x_URLLC = n_URLLC_total_vec; % normalmente igual a vec_nD2D

for a = 1:nAlgo
    plot(x_URLLC, n_URLLC_satis_vs_nD2D(a,:), '-d', 'LineWidth',1.5, ...
        'DisplayName', algorithms{a});
end
xlabel('Número de Pares D2D (URLLC)');
ylabel('Número de D2D URLLC Atendidos');
title('Atendimento URLLC x Número de Pares D2D');
legend('Location','best');
hold off;

fprintf('\nVarreduras concluídas e gráficos gerados.\n');

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
                [bestAssign, bestOF, convFit, convObj] = ...
                    RRA_ABC_Multicell(problem, V_viabilidade, params);
            case 'OABC'
                [bestAssign, bestOF, convFit, convObj] = ...
                    RRA_OABC_Multicell(problem, V_viabilidade, params);
            case 'WOA'
                [bestAssign, bestOF, convFit, convObj] = ...
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
toc;