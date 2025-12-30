% runTest_ABC_OABC_WOA_HO_GWO_AOA.m
clear; close all; clc;
tic

%% 1. Parâmetros de simulação (Artigo A1)
params.nCUE       = 5;
params.nD2D       = 5;
params.txPowerCUE = 23;   % CUE_Tx = 10dBm (A1)
params.txPowerD2D = 20;   % DUE_Tx = -10dBm (A1)
params.gammaCUE   = 10;
params.gammaD2D   = 10;
params.nRuns      = 10;
params.nRBs       = 10;
params.eNBradius  = 500;
params.maxD2DDist = 50;
params.maxIter    = 250;
params.colonySize = 30;
params.limit      = 5;
params.NLsize     = 10;

% Parâmetros específicos do AOA
params.AOA_C3 = 0.5; % ajustar conforme necessário
params.AOA_C4 = 0.5; % ajustar conforme necessário

%% 2. Parâmetros físicos
phys = struct();
phys.freq           = 1e9;
phys.noiseDensity   = -174;
phys.plCell.a       = 128.1; phys.plCell.b = 37.6;
phys.plD2D.a        = 148;   phys.plD2D.b = 40;
phys.shadowStdCell  = 10;    phys.shadowStdD2D = 12;
phys.gainTx_eNB     = 18;    phys.gainRx_UE = 0;
phys.bwSubch        = 1.44e6; % espaçamento de 60kHz → 1.44MHz por subch
phys.bwRB           = 1.44e6;
params.physParams   = phys;

%% 3. Inicialização do problema
problem = setup_Multicell_Scenario(params);

%% NOVA ETAPA: ASSOCIAÇÃO DE USUÁRIOS
problem = perform_User_Association(problem, params);

%% 4. Diagnóstico de soluções aleatórias
nSamples = 500;
objs = zeros(nSamples,1);
for k = 1:nSamples
    samp = randi([0, params.nCUE], 1, params.nD2D);
    objs(k) = evaluateSolution(samp, problem, params);
end
figure;
histogram(objs/1e6,30);
xlabel('Throughput [Mbps]');
ylabel('Frequência');
title('Fitness de Soluções Aleatórias');

%% 5. ABC puro
rng('shuffle');
[~, ~, convABC] = RRA_ABC_FullCSI_Llerena(problem, params);
allConvABC = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RRA_ABC_FullCSI_Llerena(problem, params);
    allConvABC(i,:) = conv;
end
meanConvABC = mean(allConvABC, 1);

rng(1234,'twister');
resABC = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RRA_ABC_FullCSI_Llerena(problem, params);
    resABC(i) = of;
end
bestABC = max(resABC)/1e6;
meanABC = mean(resABC)/1e6;
stdABC  = std(resABC)/1e6;

%% 6. OABC (ABC + OBL)
rng('shuffle');
[~, ~, convOABC] = RRA_OABC_FullCSI(problem, params);
allConvOABC = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RRA_OABC_FullCSI(problem, params);
    allConvOABC(i,:) = conv;
end
meanConvOABC = mean(allConvOABC, 1);

rng(1234,'twister');
resOABC = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RRA_OABC_FullCSI(problem, params);
    resOABC(i) = of;
end
bestOABC = max(resOABC)/1e6;
meanOABC = mean(resOABC)/1e6;
stdOABC  = std(resOABC)/1e6;

%% 7. WOA
rng('shuffle');
params.woaAgents  = params.colonySize;
params.woaMaxIter = params.maxIter;
[~, ~, convWOA] = RRA_WOA_FullCSI(problem, params);
allConvWOA = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RRA_WOA_FullCSI(problem, params);
    allConvWOA(i,:) = conv;
end
meanConvWOA = mean(allConvWOA, 1);

rng(1234,'twister');
resWOA = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RRA_WOA_FullCSI(problem, params);
    resWOA(i) = of;
end
bestWOA = max(resWOA)/1e6;
meanWOA = mean(resWOA)/1e6;
stdWOA  = std(resWOA)/1e6;

%% 8. HO
rng('shuffle');
params.hoAgents  = params.colonySize;
params.hoMaxIter = params.maxIter;
[~, ~, convHO] = RRA_HO_FullCSI(problem, params);
allConvHO = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RRA_HO_FullCSI(problem, params);
    allConvHO(i,:) = conv;
end
meanConvHO = mean(allConvHO, 1);

rng(1234,'twister');
resHO = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RRA_HO_FullCSI(problem, params);
    resHO(i) = of;
end
bestHO = max(resHO)/1e6;
meanHO = mean(resHO)/1e6;
stdHO  = std(resHO)/1e6;

%% 9. GWO (IGWO)
rng('shuffle');
params.gwoAgents  = params.colonySize;
params.gwoMaxIter = params.maxIter;
[~, ~, convGWO] = RRA_GWO_FullCSI(problem, params);
allConvGWO = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RRA_GWO_FullCSI(problem, params);
    allConvGWO(i,:) = conv;
end
meanConvGWO = mean(allConvGWO, 1);

rng(1234,'twister');
resGWO = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RRA_GWO_FullCSI(problem, params);
    resGWO(i) = of;
end
bestGWO = max(resGWO)/1e6;
meanGWO = mean(resGWO)/1e6;
stdGWO  = std(resGWO)/1e6;

%% 10. AOA (Archimedes Optimization Algorithm)
rng('shuffle');
[~, ~, convAOA] = RA_AOA_FullCSI(problem, params);
allConvAOA = zeros(params.nRuns, params.maxIter);
for i = 1:params.nRuns
    [~, ~, conv] = RA_AOA_FullCSI(problem, params);
    allConvAOA(i,:) = conv;
end
meanConvAOA = mean(allConvAOA, 1);

rng(1234,'twister');
resAOA = zeros(params.nRuns,1);
for i = 1:params.nRuns
    [~, of] = RA_AOA_FullCSI(problem, params);
    resAOA(i) = of;
end
bestAOA = max(resAOA)/1e6;
meanAOA = mean(resAOA)/1e6;
stdAOA  = std(resAOA)/1e6;

%% 11. Impressão de resultados
time = toc;
fprintf('--- Melhores Throughputs, Média e Desvios (%d runs) ---\n', params.nRuns);
fprintf('ABC  : %.2f  Média: %.2f  (±%.2f)\n', bestABC,  meanABC,  stdABC);
fprintf('OABC : %.2f  Média: %.2f  (±%.2f)\n', bestOABC, meanOABC, stdOABC);
fprintf('WOA  : %.2f  Média: %.2f  (±%.2f)\n', bestWOA, meanWOA, stdWOA);
fprintf('HO   : %.2f  Média: %.2f  (±%.2f)\n', bestHO,  meanHO,  stdHO);
fprintf('GWO  : %.2f  Média: %.2f  (±%.2f)\n', bestGWO, meanGWO, stdGWO);
fprintf('AOA  : %.2f  Média: %.2f  (±%.2f)\n', bestAOA, meanAOA, stdAOA);
fprintf('Tempo total: %.1f s\n', time);

%% 12. Plot de convergência
meanConvWOA_thru = -meanConvWOA;  % inverter se necessário
figure; hold on;
plot(meanConvABC/1e6,     'b-','LineWidth',1.5);
plot(meanConvOABC/1e6,    'm-','LineWidth',1.5);
plot(meanConvWOA_thru/1e6,'g-.','LineWidth',1.5);
plot(meanConvHO/1e6,      'c:','LineWidth',1.5);
plot(meanConvGWO/1e6,     'k--','LineWidth',1.5);
plot(meanConvAOA/1e6,     'r-.','LineWidth',1.5);
hold off;
legend('ABC','OABC','WOA','HO','GWO','AOA','Location','best');
xlabel('Iteração'); ylabel('Throughput [Mbps]');
title('Convergência comparativa');
grid on;

%% 13. Overhead de sinalização (em bits)
maxD2DPairs = 20;
d2d_list    = 1:maxD2DPairs;
bitsPerPair = ceil(log2(params.nCUE + 1));
overhead_D2D = d2d_list * bitsPerPair;

maxCUE   = 10;
cue_list = 1:maxCUE;
overhead_CUE = params.nD2D * ceil(log2(cue_list + 1));

figure;
plot(d2d_list, overhead_D2D, '-o','LineWidth',1.5); hold on;
plot(cue_list, overhead_CUE, '-s','LineWidth',1.5);
hold off;
xlabel('Nº de pares D2D / Nº de CUEs');
ylabel('Overhead (bits)');
title('Overhead vs Nº de pares D2D e Nº de CUEs');
legend('Overhead D2D','Overhead CUE','Location','northwest');
grid on;
