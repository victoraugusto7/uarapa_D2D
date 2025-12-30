function problem = setup_Multicell_Scenario(params)
% setup_Multicell_Scenario: Configura o cenário 3x3 Hexagonal (Estilo FU-2024)
%   Cria 9 células em disposição "Honeycomb" e distribui usuários.

    %% 0. Reprodutibilidade
    if isfield(params,'seedTopology')
        rng(params.seedTopology);
    end

    %% 1. Copia parâmetros básicos
    problem.nCUE      = params.nCUE;
    problem.nD2D      = params.nD2D;
    problem.eNBradius = params.eNBradius;
    problem.maxD2DDist= params.maxD2DDist;
    problem.txPowerCUE = params.txPowerCUE;
    problem.txPowerD2D = params.txPowerD2D;
    problem.gammaCUE   = params.gammaCUE;
    problem.gammaD2D   = params.gammaD2D;
    problem.bwRB      = params.physParams.bwRB;
    problem.nRBs      = params.nRBs;
    problem.physParams = params.physParams;
    
    % Campo de reuso
    problem.reuse_factor = params.reuse_factor; 

    %% 2. Geração de Posições (Topologia Hexagonal 3x3 - FU-2024)
    
    R = problem.eNBradius;
    % Distância Horizontal entre centros (width of hexagon)
    D = R * sqrt(3); 
    % Distância Vertical entre linhas (3/2 * R)
    V = 1.5 * R;
    
    if isfield(params, 'num_cells')
        nCells = params.num_cells;
    else
        nCells = 9; 
    end
    
    pos_MBSs = zeros(nCells, 2);
    
    % Construção da Grade Hexagonal 3x3 (Honeycomb)
    % Linha 1 (Topo): 3 células
    % Linha 2 (Meio): 3 células (Deslocada para a direita por D/2)
    % Linha 3 (Baixo): 3 células
    
    % Coordenadas manuais para garantir o visual "Puzzle" perfeito
    
    % Linha Superior (y = +V)
    pos_MBSs(1,:) = [-D, V];
    pos_MBSs(2,:) = [ 0, V];
    pos_MBSs(3,:) = [ D, V];
    
    % Linha do Meio (y = 0) -> Deslocada +D/2
    pos_MBSs(4,:) = [-D + D/2, 0];
    pos_MBSs(5,:) = [ 0 + D/2, 0];
    pos_MBSs(6,:) = [ D + D/2, 0];
    
    % Linha Inferior (y = -V)
    pos_MBSs(7,:) = [-D, -V];
    pos_MBSs(8,:) = [ 0, -V];
    pos_MBSs(9,:) = [ D, -V];
    
    % Se nCells for diferente de 9, trunca ou expande (mas focamos em 9)
    problem.pos_MBSs = pos_MBSs(1:nCells, :);
    problem.nCells = nCells;
    
    %% 3. Lógica de Frequência (Reuse-3 Padrão)
    
    problem.PRBs_available_local = cell(problem.nCells, 1);
    R_set_size = floor(problem.nRBs / 3);
    
    % Padrão de cores para evitar vizinhos iguais no Grid Hexagonal
    % Linha 1: 1, 2, 3
    % Linha 2: 3, 1, 2 (Shifted)
    % Linha 3: 2, 3, 1 (Shifted)
    MBS_set_map = [1, 2, 3, 3, 1, 2, 2, 3, 1]; 
    
    if problem.reuse_factor == 1
        % Reuse 1: Todos usam tudo
        for b = 1:problem.nCells
             problem.PRBs_available_local{b} = 1:problem.nRBs;
        end
        problem.MBS_PRB_Set_ID = ones(1, problem.nCells);
    else
        % Reuse 3
        problem.MBS_PRB_Set_ID = MBS_set_map(1:nCells);
        for b = 1:problem.nCells
            set_idx = problem.MBS_PRB_Set_ID(b);
            start_idx = (set_idx - 1) * R_set_size + 1;
            end_idx = set_idx * R_set_size;
            problem.PRBs_available_local{b} = start_idx:end_idx;
        end
    end
    
    %% 4. Distribuição de Usuários (Ajustada para o Grid)
    
    % Raio de Simulação que cobre todo o grid 3x3
    % Centro aproximado em (0,0), raio deve cobrir até as bordas (~2*D)
    simulationRadius = 2.5 * R; 
    
    % --- Posições dos CUEs ---
    % Gera pontos, mas filtra para garantir que caiam DENTRO de alguma célula
    problem.pos_CUEs = zeros(problem.nCUE, 2);
    count = 0;
    while count < problem.nCUE
        pt = getPointsInCircle(simulationRadius, 1);
        % Verifica se está dentro de algum hexágono (distância ao centro < R)
        % Simplificação: Verifica distância ao MBS mais próximo
        dists = vecnorm(problem.pos_MBSs - pt, 2, 2);
        if min(dists) <= R
            count = count + 1;
            problem.pos_CUEs(count, :) = pt;
        end
    end
    
    % --- Posições dos D2Ds ---
    problem.pos_D2DTx = zeros(problem.nD2D, 2);
    problem.pos_D2DRx = zeros(problem.nD2D, 2);
    
    count = 0;
    while count < problem.nD2D
        pt = getPointsInCircle(simulationRadius, 1);
        dists = vecnorm(problem.pos_MBSs - pt, 2, 2);
        if min(dists) <= R
            count = count + 1;
            problem.pos_D2DRx(count, :) = pt;
            shift = getPointsInCircle(problem.maxD2DDist, 1);
            problem.pos_D2DTx(count, :) = pt + shift;
        end
    end

    % --- 2e. Etiquetar Usuários (Cenário Híbrido: eMBB + URLLC) ---
    fprintf('Etiquetando usuários (eMBB e URLLC)...\n');
    
    % CUEs: params.eMBB_ratio são eMBB (Comandantes)
    num_eMBB_CUEs = round(params.eMBB_ratio * problem.nCUE);
    problem.is_eMBB_CUE = zeros(problem.nCUE, 1);
    indices_CUE = randperm(problem.nCUE, num_eMBB_CUEs);
    problem.is_eMBB_CUE(indices_CUE) = 1; 
    fprintf('  - CUEs: %d eMBB (Comandantes) / %d Comuns\n', num_eMBB_CUEs, problem.nCUE - num_eMBB_CUEs);

    % D2Ds: Todos são URLLC (Agentes)
    problem.is_URLLC_D2D = ones(problem.nD2D, 1); 
    problem.is_eMBB_D2D  = zeros(problem.nD2D, 1); 
    fprintf('  - D2Ds: %d URLLC (Agentes de Segurança)\n', problem.nD2D);

    % --- 2f. Parâmetros de QoS ---
    if isfield(params, 'eMBB_Req_Rate_Test')
        problem.eMBB_Req_Rate = params.eMBB_Req_Rate_Test; 
    else
        problem.eMBB_Req_Rate = params.eMBB_Req_Rate;  
    end
    problem.eMBB_Req_BER  = params.eMBB_Req_BER;  
    problem.eMBB_Req_Lat  = params.eMBB_Req_Lat; 
    
    if isfield(params, 'URLLC_Req_Rate')
        problem.URLLC_Req_Rate = params.URLLC_Req_Rate;
        problem.URLLC_Req_BER  = params.URLLC_Req_BER;
        problem.URLLC_Req_Lat  = params.URLLC_Req_Lat;
    else
        problem.URLLC_Req_Rate = 1e6; 
        problem.URLLC_Req_BER  = 1e-6;
        problem.URLLC_Req_Lat  = 20e-3;
    end
    
    % SINR Mínimo
    problem.eMBB_Req_SINR_dB   = 10.0; 
    problem.Common_Req_SINR_dB = 7.0; 
    problem.URLLC_Req_SINR_dB  = 7.0; 

    %% 5. Cálculo de Ganhos de Canal (Globais)
    
    % CUE -> MBSs
    problem.gain_CUE_to_MBSs = zeros(problem.nCUE, problem.nCells);
    for i = 1:problem.nCUE
        for b = 1:problem.nCells
            distKm = norm(problem.pos_CUEs(i,:) - problem.pos_MBSs(b,:))/1e3;
            problem.gain_CUE_to_MBSs(i, b) = ...
                computeGain(distKm, 'cellular', problem.physParams);
        end
    end

    % D2D Tx -> D2D Rx
    problem.gain_D2D_pair = zeros(problem.nD2D, 1);
    for i = 1:problem.nD2D
        distKm = norm(problem.pos_D2DTx(i,:) - problem.pos_D2DRx(i,:))/1e3;
        problem.gain_D2D_pair(i) = ...
            computeGain(distKm, 'D2D', problem.physParams);
    end

    % Interferências
    problem.gain_CUE_to_D2DRx = zeros(problem.nD2D, problem.nCUE);
    for i = 1:problem.nD2D
        for j = 1:problem.nCUE
            distKm = norm(problem.pos_CUEs(j,:) - problem.pos_D2DRx(i,:))/1e3;
            problem.gain_CUE_to_D2DRx(i, j) = ...
                computeGain(distKm, 'cellular', problem.physParams);
        end
    end
    
    problem.gain_D2DTx_to_MBSs = zeros(problem.nD2D, problem.nCells);
    for i = 1:problem.nD2D 
        for b = 1:problem.nCells
            distKm = norm(problem.pos_D2DTx(i,:) - problem.pos_MBSs(b,:))/1e3;
            problem.gain_D2DTx_to_MBSs(i, b) = ...
                computeGain(distKm, 'cellular', problem.physParams);
        end
    end
    
    problem.gain_D2DTx_to_D2DRx = zeros(problem.nD2D, problem.nD2D);
    for i = 1:problem.nD2D 
        for k = 1:problem.nD2D 
            if i == k
                problem.gain_D2DTx_to_D2DRx(i, k) = 0; 
                continue; 
            end
            distKm = norm(problem.pos_D2DTx(k,:) - problem.pos_D2DRx(i,:))/1e3;
            problem.gain_D2DTx_to_D2DRx(i, k) = ...
                computeGain(distKm, 'D2D', problem.physParams);
        end
    end
    
    fprintf('Cenário Multicelular (%d Células, Reuso %d, Layout FU-2024) criado.\n', ...
        problem.nCells, problem.reuse_factor);
end