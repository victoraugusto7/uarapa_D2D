function total_rate = power_control_fitness_multicell(P_vector, X, problem, params)
% POWER_CONTROL_FITNESS_MULTICELL: Objetivo do GWO - Maximizar Taxa.
%   Adaptado para MATRIZ X (3 x N_RBs).
%   P_vector: Vetor [N_total_users x 1] de potências (W).
%   X: Matriz de Alocação [3 x N_RBs].

    % 1. Inicialização e Parâmetros
    total_rate = 0;
    
    % Converte SINR Mínimo para escala linear
    xi_min_c_lin  = 10^(problem.Common_Req_SINR_dB / 10); % 7 dB
    xi_min_embb_lin = 10^(problem.eMBB_Req_SINR_dB / 10); % 10 dB

    % 2. Loop sobre cada PRB
    for r = 1:problem.nRBs
        
        % Recupera os utilizadores neste PRB (Coluna r, todas as 3 linhas)
        users_in_prb = X(:, r); % Vetor 3x1
        
        % Se o PRB estiver vazio (todas as linhas = 0), pula
        if all(users_in_prb == 0)
            continue;
        end
        
        % --- CÁLCULO PARA CADA USUÁRIO NO PRB 'r' ---
        % (CUE na linha 1, D2Ds nas linhas 2 e 3)
        
        for row = 1:3
            k = users_in_prb(row);
            
            if k == 0
                continue; % Slot vazio
            end
            
            % Calcula SINR e Taxa para o par (k, r)
            % (A função calculate_SINR... já foi atualizada para ler a matriz X)
            [SINR_kr, R_kr] = calculate_SINR_Rate_for_PRB(k, r, P_vector, X, problem, params);

            % 3. Imposição da Restrição de SINR Mínimo (Penalidade GWO)
            is_cue = (k <= problem.nCUE);
            
            if is_cue && problem.is_eMBB_CUE(k)
                SINR_min_k = xi_min_embb_lin;
            else
                SINR_min_k = xi_min_c_lin;
            end

            % Penalidade para o GWO se violar o SINR mínimo
            if SINR_kr < SINR_min_k
                R_kr = -1e9; % Penalidade negativa severa (GWO foge disso)
            end

            total_rate = total_rate + R_kr;
        end
    end
end