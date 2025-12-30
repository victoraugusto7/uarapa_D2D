function P_star = optimal_power_control(X, problem, params)
% OPTIMAL_POWER_CONTROL: Orquestra o I-GWO para encontrar o vetor de potência P_star.
%   X: Vetor Solução de Alocação de PRBs (fixo).
%   P_star: Vetor [N_total_users x 1] de potências ótimas (W).
%
% Requer: I-GWO.m e initialization_g.m no diretório.

    % N_total_users (Dimensão do vetor P: 75x1)
    nTotalUsers = problem.nCUE + problem.nD2D; 
    dim = nTotalUsers;

   % fprintf('  [Loop Interno: Executando GWO para Controle de Potência]...\n');

    %% 1. Definir Limites do Espaço de Busca [lb, ub] (em Watts)
    
    P_max_c_watts = 10^((params.txPowerCUE - 30) / 10);
    P_max_d_watts = 10^((params.txPowerD2D - 30) / 10);
    P_min_epsilon = 1e-10; % Um valor muito pequeno, mas diferente de zero
    
    lb = zeros(dim, 1) + P_min_epsilon; % <-- CORRIGIDO
    ub = zeros(dim, 1); % P_max
    
    for k = 1:nTotalUsers
        if k <= problem.nCUE % CUE
            ub(k) = P_max_c_watts;
        else % D2D
            ub(k) = P_max_d_watts;
        end
    end
    
    % Transpõe lb e ub para o formato de linha exigido pelo solver IGWO
    lb_row = lb';
    ub_row = ub';
    
    %% 2. Definir a Função Objetivo (fobj)
    % Minimizar o negativo da taxa é equivalente a maximizar a taxa.
    fobj = @(P_vector) -power_control_fitness_multicell(P_vector, X, problem, params);

    %% 3. Parâmetros e Chamada do I-GWO
    
    gwoAgents = params.colonySize_PC; 
    gwoMaxIter = params.maxIter_PC;   
    
    % --- CHAMADA REAL DO I-GWO ---
    % O I-GWO é um minimizador, e fobj minimiza (-taxa).
    [~, Alpha_pos, ~] = IGWO(dim, gwoAgents, gwoMaxIter, lb_row, ub_row, fobj);
    
    % P_star deve ser o vetor P otimizado (em coluna)
    P_star = Alpha_pos'; 
    
   % fprintf('  [Loop Interno: PC concluído. Potência ótima P* encontrada pelo I-GWO]\n');
end