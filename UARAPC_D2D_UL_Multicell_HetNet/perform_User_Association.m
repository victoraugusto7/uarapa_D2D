function problem = perform_User_Association(problem, ~)
% perform_User_Association: Etapa 1 (UA) - FINAL
%   Associa CUEs e D2Ds às MBSs com base na Distância Mínima e aplica um
%   Filtro de Cobertura para remover usuários fora de alcance.

    fprintf('Executando Etapa 1: Associação de Usuários (UA) por Proximidade...\n');

    % --- 1. Inicialização ---
    problem.assoc_CUE = zeros(problem.nCUE, 1);
    problem.assoc_D2D = zeros(problem.nD2D, 1);
    
    nCells = problem.nCells;
    pos_MBSs = problem.pos_MBSs;
    R_cell_max = problem.eNBradius * 1.5; % Raio de Cobertura Máximo (750m)

    % --- 2. Associação de CUEs (CUE i -> MBS mais próxima) ---
    for i = 1:problem.nCUE
        dist_min = inf;
        b_best = 0; % Inicializa com 0 (não associado)
        
        for b = 1:nCells
            dist_ib = norm(problem.pos_CUEs(i,:) - pos_MBSs(b,:));
            
            if dist_ib < dist_min
                dist_min = dist_ib;
                b_best = b;
            end
        end
        
        % Filtro de Cobertura: Só associa se a MBS mais próxima estiver dentro do R_max
        if dist_min <= R_cell_max
            problem.assoc_CUE(i) = b_best;
        end
    end

    % --- 3. Associação de D2Ds (DUE-Tx j -> MBS mais próxima) ---
    for j = 1:problem.nD2D
        dist_min = inf;
        b_best = 0; % Inicializa com 0 (não associado)
        
        for b = 1:nCells
            dist_jb = norm(problem.pos_D2DTx(j,:) - pos_MBSs(b,:));
            
            if dist_jb < dist_min
                dist_min = dist_jb;
                b_best = b;
            end
        end
        
        % Filtro de Cobertura: Só associa se a MBS mais próxima estiver dentro do R_max
        if dist_min <= R_cell_max
            problem.assoc_D2D(j) = b_best;
        end
    end
    
    % Diagnóstico:
    unassoc_cue = sum(problem.assoc_CUE == 0);
    unassoc_d2d = sum(problem.assoc_D2D == 0);
    fprintf('Associação de Usuários (UA) concluída.\n');
    fprintf('  - CUEs rejeitados (fora de cobertura): %d\n', unassoc_cue);
    fprintf('  - D2Ds rejeitados (fora de cobertura): %d\n', unassoc_d2d);
end