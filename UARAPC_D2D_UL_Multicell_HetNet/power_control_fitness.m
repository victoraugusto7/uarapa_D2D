function total_rate = power_control_fitness(P_vector, X, problem, params)
% POWER_CONTROL_FITNESS: Função de Fitness para o Loop Interno (MHA-PC)
%   Objetivo: Maximizar a Taxa Total (Sum Rate) para um dado vetor P (potência).
%   P_vector: Vetor [N_total_users x 1] de potências em Watts (W).
%   X: Vetor Solução de Alocação de PRBs (fixo).

    % 1. Inicialização e Parâmetros
    total_rate = 0;
    
    % Converte o vetor de potências contínuas para a escala permitida [Pmin, Pmax]
    P_max_d_watts = 10^((params.txPowerD2D - 30) / 10);
    P_max_c_watts = 10^((params.txPowerCUE - 30) / 10);

    % Restrições do Domínio: Aplica a restrição de Potência Máxima
    P_vector(P_vector > P_max_c_watts) = P_max_c_watts; % Para CUEs (índices 1 a N)
    P_vector(P_vector > P_max_d_watts) = P_max_d_watts; % Para D2Ds (índices N+1 a N+M)
    P_vector(P_vector < 0) = 0;

    % Parâmetros de QoS (em linear)
    xi_min_c_lin = 10^(problem.Common_Req_SINR_dB / 10); % 7 dB
    xi_min_embb_lin = 10^(problem.eMBB_Req_SINR_dB / 10); % 10 dB

    % 2. Mapear Alocação (X)
    [~, ~, ~] = map_allocation_from_X(X, problem);
    
    % 3. Loop sobre cada PRB (o menor bloco de recurso)
    for r = 1:problem.nRBs
        
        % Usuários que utilizam o PRB 'r'
        users_in_prb = find(X == r);
        if isempty(users_in_prb)
            continue;
        end
        
        % --- CÁLCULO SINR E TAXA PARA CADA USUÁRIO NO PRB 'r' ---
        
        % 3a. Cálculo da Interferência (A Parte Mais Crítica)
        % Esta função (auxiliar) será desenvolvida para calcular I_total para o PRB 'r'
        % dado o vetor de potências P_vector.
        I_map = calculate_interference_matrix(X, P_vector, problem);

        for k_idx = 1:length(users_in_prb)
            k = users_in_prb(k_idx); % Índice global do usuário
            
            % Calcula o SINR
            [SINR_kr, R_kr] = calculate_SINR_Rate_for_PRB(k, r, P_vector, I_map, problem);

            % 4. Imposição da Restrição de SINR Mínimo
            % (Penaliza a taxa se a potência P_k for insuficiente para o SINR mínimo)
            
            is_embb_cue = (k <= problem.nCUE && problem.is_eMBB_CUE(k));
            
            if is_embb_cue
                SINR_min_k = xi_min_embb_lin;
            else
                SINR_min_k = xi_min_c_lin;
            end

            if SINR_kr < SINR_min_k
                R_kr = 0; % Penalidade (Taxa zero para este PRB)
            end

            total_rate = total_rate + R_kr;
        end
    end
end