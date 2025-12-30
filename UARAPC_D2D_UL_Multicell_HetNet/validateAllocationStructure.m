function is_valid = validateAllocationStructure(X, V_viabilidade, problem, ~)
% VALIDATEALLOCATIONSTRUCTURE: Valida a Matriz de Solução X (3 x N_RBs).
%   Verifica regras rígidas de rádio e estrutura Underlay.
%
%   Entrada:
%       X: Matriz [3 x N_RBs]
%          Linha 1: CUEs
%          Linhas 2-3: D2Ds
%
%   Retorna:
%       is_valid: true se a solução respeita todas as regras estruturais.

    is_valid = true; 
    nCUEs = problem.nCUE;
    nTotalUsers = nCUEs + problem.nD2D;
    nRBs = problem.nRBs;

    % --- 1. Validação Por PRB (Colunas) ---
    for r = 1:nRBs
        % A. Verifica a Linha 1 (CUEs)
        u_cue = X(1, r);
        if u_cue > 0
            % Deve ser um CUE válido
            if u_cue > nCUEs
                is_valid = false; return; % Erro: D2D na linha de CUE
            end
            
            % Valida Regras de Rádio para o CUE
            if ~check_radio_rules(u_cue, r, V_viabilidade, problem)
                is_valid = false; return;
            end
        end
        
        % B. Verifica as Linhas 2 e 3 (D2Ds - Underlay)
        for row = 2:3
            u_d2d = X(row, r);
            if u_d2d > 0
                % Deve ser um D2D válido
                if u_d2d <= nCUEs
                    is_valid = false; return; % Erro: CUE na linha de D2D
                end
                
                % Regra de Underlay: D2D só pode usar PRB se houver um CUE (Linha 1)
                if u_cue == 0
                    is_valid = false; return; % Violação de Underlay (D2D sozinho no PRB)
                end
                
                % Valida Regras de Rádio para o D2D
                if ~check_radio_rules(u_d2d, r, V_viabilidade, problem)
                    is_valid = false; return;
                end
            end
        end
        
        % C. Verifica Duplicatas no Mesmo PRB
        % (Não podemos ter o mesmo D2D duas vezes no mesmo PRB, linhas 2 e 3)
        users_in_prb = X(X(:,r) > 0, r);
        if length(unique(users_in_prb)) < length(users_in_prb)
            is_valid = false; return; % Usuário duplicado no mesmo PRB
        end
    end
    
    % --- 2. Validação de Quantidade (Por Usuário) ---
    % Conta quantas vezes cada usuário aparece na matriz inteira
    % (Isso pode ser pesado, então fazemos de forma vetorizada se possível)
    
    % Achata a matriz para vetor para contagem rápida
    X_flat = X(:);
    X_flat = X_flat(X_flat > 0); % Apenas usuários ativos
    
    % Conta ocorrências
    counts = histcounts(X_flat, 1:(nTotalUsers+1)); 
    
    for k = 1:nTotalUsers
        count_k = counts(k);
        
        is_cue = (k <= nCUEs);
        is_embb_cue = (is_cue && problem.is_eMBB_CUE(k));
        
        % Regra para Comuns (CUEs e D2Ds): Máximo 1 PRB
        if ~is_embb_cue
            if count_k > 1
                is_valid = false; return;
            end
        end
        
        % Regra para eMBB: Nenhuma restrição de máximo (Otimização Pura)
        % Opcional: Pode-se exigir count_k > 0 para evitar inanição total,
        % mas a função de fitness (QoS) já penaliza isso com taxa zero.
    end

end

%% Função Auxiliar Local
function ok = check_radio_rules(u, r, V, problem)
    ok = true;
    
    % 1. Viabilidade (CPS Fase 2)
    if V(u, r) == 0
        ok = false; return;
    end
    
    % 2. Reuso 3 (Pertence à MBS correta?)
    if u <= problem.nCUE
        b = problem.assoc_CUE(u);
    else
        b = problem.assoc_D2D(u - problem.nCUE);
    end
    
    if b == 0
        ok = false; return; % Usuário não associado
    end
    
    PRB_set = problem.PRBs_available_local{b};
    if ~ismember(r, PRB_set)
        ok = false; return;
    end
end