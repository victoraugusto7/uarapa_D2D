function [SINR_kr, R_kr] = calculate_SINR_Rate_for_PRB(k, r, P_vector, X, problem, params)
% CALCULATE_SINR_RATE_FOR_PRB: Calcula SINR e Taxa (Shannon) para o par (k, r).
%   Suporta a estrutura de Matriz X (3 linhas) para Underlay.
%
%   k: Índice global do utilizador
%   r: Índice do PRB
%   P_vector: Vetor de potências (Watts)
%   X: Matriz de Alocação [3 x N_RBs]

    %% 1. Parâmetros e Ruído
    % Ruído térmico (Watts)
    noise_density_dbW_Hz = params.physParams.noiseDensity - 30;
    noise_watts = 10^(noise_density_dbW_Hz / 10) * problem.bwRB;

    % Interferência Intercelular (CCI) Média (Baseline para Reuso 3)
    % Como discutido, usamos um valor médio para representar a interferência das células distantes
    I_CCI_media_watts = 10^((params.Interf_CCI_media_dBm - 30) / 10);
    
    % Ganho útil (Sinal)
    % Identifica se é CUE->MBS ou D2D->D2D
    [G_k_util, b_servidora] = get_channel_gain_util_and_mbs(k, problem);
    
    % Potência do utilizador k
    P_k = P_vector(k);
    
    %% 2. Sinal Útil
    S_kr = P_k * G_k_util;

    %% 3. Cálculo da Interferência Intracelular (I_intra)
    % Quem mais está neste PRB r? (Verifica as 3 linhas da coluna r)
    users_in_prb = X(:, r); % Vetor 3x1 [CUE; D2D1; D2D2]
    
    I_intra = 0;
    
    for row = 1:3
        u_interferer = users_in_prb(row);
        
        % Se for 0 (vazio) ou for o próprio utilizador k, ignora
        if u_interferer == 0 || u_interferer == k
            continue;
        end
        
        % Temos um interferente 'u' no mesmo PRB!
        P_u = P_vector(u_interferer);
        
        % Calcula o ganho de interferência de u para o receptor de k
        G_u_to_k_receiver = get_interference_gain(u_interferer, k, b_servidora, problem);
        
        I_intra = I_intra + (P_u * G_u_to_k_receiver);
    end
    
    %% 4. SINR e Taxa
    I_total = noise_watts + I_CCI_media_watts + I_intra;
    
    SINR_kr = S_kr / I_total;
    
    % Taxa de Shannon (bps)
    R_kr = problem.bwRB * log2(1 + SINR_kr);
end

%% --- Funções Auxiliares Locais ---

function [G_util, b_serv] = get_channel_gain_util_and_mbs(k, problem)
    if k <= problem.nCUE
        % É um CUE
        b_serv = problem.assoc_CUE(k);
        G_util = problem.gain_CUE_to_MBSs(k, b_serv);
    else
        % É um D2D
        d2d_idx = k - problem.nCUE;
        b_serv = problem.assoc_D2D(d2d_idx); % Apenas para referência de localização
        G_util = problem.gain_D2D_pair(d2d_idx);
    end
end

function G_int = get_interference_gain(u_tx, k_victim, b_victim, problem)
% Calcula o ganho do canal entre o Transmissor u e o Receptor do utilizador k
    
    is_victim_CUE = (k_victim <= problem.nCUE);
    is_interferer_CUE = (u_tx <= problem.nCUE);
    
    if is_victim_CUE
        % Vítima é CUE -> Receptor é a MBS (b_victim)
        if is_interferer_CUE
             % CUE interferindo em MBS (Raro intra-celula no mesmo PRB, mas possível se for Overlay)
             % Neste modelo Underlay, CUEs ocupam a linha 1, então não deve haver outro CUE aqui.
             G_int = 0; 
        else
             % D2D interferindo em MBS
             d2d_idx = u_tx - problem.nCUE;
             G_int = problem.gain_D2DTx_to_MBSs(d2d_idx, b_victim);
        end
    else
        % Vítima é D2D -> Receptor é o D2D-Rx de k
        d2d_victim_idx = k_victim - problem.nCUE;
        
        if is_interferer_CUE
            % CUE interferindo em D2D-Rx
            G_int = problem.gain_CUE_to_D2DRx(d2d_victim_idx, u_tx);
        else
            % D2D interferindo em D2D-Rx
            d2d_tx_idx = u_tx - problem.nCUE;
            G_int = problem.gain_D2DTx_to_D2DRx(d2d_victim_idx, d2d_tx_idx);
        end
    end
end