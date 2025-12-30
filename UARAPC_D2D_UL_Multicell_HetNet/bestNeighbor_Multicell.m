function newSol = bestNeighbor_Multicell(sol, NL, V_viabilidade, problem, params)
% BESTNEIGHBOR_MULTICELL: Operador de Vizinhança (Versão Matriz 3xR).
%   Modifica aleatoriamente uma posição da matriz X (3 x N_RBs).
%   - Linha 1: Altera CUEs.
%   - Linhas 2 e 3: Altera D2Ds (Underlay).

    nCUEs = problem.nCUE;
    nD2Ds = problem.nD2D;
    nTotalUsers = nCUEs + nD2Ds;
    nRBs = problem.nRBs;
    
    % Domínios separados
    dom_CUE = [0, 1:nCUEs];
    dom_D2D = [0, (nCUEs+1):nTotalUsers];
    
    % Melhor solução local começa como a atual
    bestSol = sol;
    
    % Avalia a solução atual (Fitness base)
    bestFit = evaluateSolution_Multicell(sol, V_viabilidade, problem, params);
    
    % Loop de Vizinhança (Tenta NL mutações)
    for t = 1:NL
        cand = sol;
        
        % 1. Escolhe uma posição aleatória na Matriz (Linha, Coluna)
        r_idx = randi(nRBs);    % Coluna (PRB)
        row_idx = randi(3);     % Linha (1=CUE, 2/3=D2D)
        
        k_old = cand(row_idx, r_idx);
        
        % 2. Escolhe um novo valor (k_new) baseado na Linha
        if row_idx == 1
            % Linha 1: Apenas CUEs
            k_new = dom_CUE(randi(length(dom_CUE)));
        else
            % Linhas 2 e 3: Apenas D2Ds
            k_new = dom_D2D(randi(length(dom_D2D)));
        end
        
        if k_new == k_old
            continue; 
        end

        % --- VALIDAÇÃO DA MUTAÇÃO (Golden Rules) ---
        is_valid_mutation = true;
        
        if k_new > 0
            % a) Identifica a MBS e Tipo do novo utilizador
            if k_new <= nCUEs 
                b = problem.assoc_CUE(k_new);
                is_embb = problem.is_eMBB_CUE(k_new);
            else 
                d2d_idx = k_new - nCUEs;
                b = problem.assoc_D2D(d2d_idx);
                is_embb = false; % D2Ds são sempre Comuns/URLLC aqui
            end
            
            % REGRA 1: Utilizador deve estar associado
            if b == 0
                is_valid_mutation = false;
            end
            
            % REGRA 2: Reuso 3 (O PRB deve pertencer à MBS do utilizador)
            if is_valid_mutation
                PRB_set_b = problem.PRBs_available_local{b};
                if ~ismember(r_idx, PRB_set_b)
                    is_valid_mutation = false;
                end
            end

            % REGRA 3: Viabilidade (Fase 2 - CPS)
            if is_valid_mutation
                if V_viabilidade(k_new, r_idx) == 0
                    is_valid_mutation = false; 
                end
            end
            
            % REGRA 4: Quantidade para Comuns (Underlay Rígido)
            % Se for Comum, não pode ter mais de 1 PRB em TODA a matriz.
            if is_valid_mutation && ~is_embb
                % Conta quantas vezes k_new aparece na matriz inteira
                count = sum(cand(:) == k_new);
                if count >= 1
                     is_valid_mutation = false;
                end
            end
            
            % REGRA 5: Evitar duplicata no mesmo PRB (CUE vs D2D)
            % Não podemos ter o mesmo usuário na Linha 1 e na Linha 2 do mesmo PRB
            if is_valid_mutation
                users_in_prb = cand(:, r_idx); % Quem já está neste PRB
                if ismember(k_new, users_in_prb)
                    is_valid_mutation = false;
                end
            end
            
            % REGRA 6 (Opcional): Underlay "Puro"
            % D2D (Linha 2/3) só pode entrar se houver CUE (Linha 1)?
            % Se quiser forçar isso, descomente abaixo:
            % if row_idx > 1 && cand(1, r_idx) == 0
            %    is_valid_mutation = false;
            % end
        end
        
        % --- Aplica a Mudança e Avalia ---
        if is_valid_mutation
            cand(row_idx, r_idx) = k_new;
            
            % Avalia a nova solução candidata
            fit = evaluateSolution_Multicell(cand, V_viabilidade, problem, params);
            
            if fit > bestFit
                bestFit = fit;
                bestSol = cand;
            end
        end
    end
    
    newSol = bestSol;
end