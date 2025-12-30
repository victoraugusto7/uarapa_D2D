function N_req = calculate_N_PRB_required(problem)
% CALCULATE_N_PRB_REQUIRED: Calcula o número de PRBs (N_k^PRB) para CUEs eMBB.
%   O cálculo é feito para todos os eMBB CUEs, com base no requisito de 50 Mbps e 10 dB.
    
    % Parâmetros eMBB
    R_target_bps = problem.eMBB_Req_Rate; % 50 Mbps
    W_RB = problem.bwRB;                 % 360 kHz
    xi_min_dB = problem.eMBB_Req_SINR_dB; % 10 dB
    
    % Conversão para Linear
    xi_min_lin = 10^(xi_min_dB / 10); % ~10.0
    
    % Taxa Máxima Teórica por PRB (com 10 dB)
    R_per_PRB_max = W_RB * log2(1 + xi_min_lin); 
    
    % N_PRBs necessários para atingir 50 Mbps
    N_req_float = R_target_bps / R_per_PRB_max;
    
    % Arredondamento para cima (garantir a taxa)
    N_req = ceil(N_req_float); % Esperado ser aproximadamente 41
    
    % Adiciona o número de PRBs necessários ao problema (para uso futuro)
    problem.N_PRBs_eMBB_req = N_req; 
end