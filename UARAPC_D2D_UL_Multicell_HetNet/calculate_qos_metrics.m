function [R_total, R_k_vector, BER_k_vector, L_k_vector] = calculate_qos_metrics(X, P_star, problem, params)
% CALCULATE_QOS_METRICS: Versão com D2D URLLC
%   - CUEs (Uplink): inclui Backhaul (L_sb) + fila.
%   - D2Ds (Direto): NÃO inclui Backhaul. Apenas fila + ar + processamento.
%   - D2Ds URLLC: usam lambda e tamanho de pacote dedicados.

    %% 1. Inicialização
    nCUEs        = problem.nCUE;
    nTotalUsers  = nCUEs + problem.nD2D;
    
    R_k_vector   = zeros(nTotalUsers, 1);
    BER_k_vector = zeros(nTotalUsers, 1);
    L_k_vector   = zeros(nTotalUsers, 1);
    
    % Parâmetros de Latência Fixos
    L_prop = params.eMBB_Prop_latency;   % propagação
    L_sb   = params.eMBB_Latency_SB;     % backhaul (apenas CUEs)
    L_proc = 1e-3;                       % processamento (estimado p/ D2D)

    %% 2. Mapeamento X -> PRBs por usuário
    [User_to_PRBs_Map, ~] = map_allocation_from_X(X, problem);
    
    %% 3. Loop Principal
    for k = 1:nTotalUsers
        
        PRBs_k = User_to_PRBs_Map{k};
        
        % Se não tem PRB, usuário morto
        if isempty(PRBs_k)
            R_k_vector(k)   = 0;
            BER_k_vector(k) = 1;
            L_k_vector(k)   = 1e12; 
            continue;
        end
        
        % 4. Taxa total e pior SINR (para BER conservadora)
        R_k_total_sum = 0;
        min_SINR_k    = inf;
        
        for r_idx = 1:length(PRBs_k)
            r = PRBs_k(r_idx);
            [SINR_kr, R_kr] = calculate_SINR_Rate_for_PRB(k, r, P_star, X, problem, params);
            R_k_total_sum = R_k_total_sum + R_kr;
            if SINR_kr < min_SINR_k
                min_SINR_k = SINR_kr;
            end
        end
        
        R_k_vector(k) = R_k_total_sum;
        
        % 5. BER (BPSK/QPSK approx)
        if min_SINR_k > 0
             BER_k_vector(k) = 0.5 * erfc(sqrt(min_SINR_k)); 
        else
             BER_k_vector(k) = 1.0;
        end
        
        % 6. Latência (M/G/1) diferenciada por tipo de usuário
        if k <= nCUEs
            % --------- CUEs ----------
            is_D2D = false;
            if problem.is_eMBB_CUE(k)
                % eMBB
                lambda_k = params.eMBB_lambda; 
                pkt_size = params.eMBB_Packet_size;
            else
                % CUE Comum (Best-effort)
                if isfield(params, 'Common_lambda')
                    lambda_k = params.Common_lambda;
                else
                    lambda_k = 400;
                end
                pkt_size = 1000;
            end
        else
            % --------- D2D ----------
            is_D2D   = true;
            d2d_idx  = k - nCUEs;
            
            if isfield(problem,'is_URLLC_D2D') && problem.is_URLLC_D2D(d2d_idx)
                % D2D URLLC: parâmetros dedicados
                if isfield(params,'URLLC_lambda')
                    lambda_k = params.URLLC_lambda;
                else
                    lambda_k = 250;
                end
                
                if isfield(params,'URLLC_Packet_size')
                    pkt_size = params.URLLC_Packet_size;
                else
                    pkt_size = 256*8; % fallback
                end
            else
                % D2D comum / best-effort
                if isfield(params, 'Common_lambda')
                    lambda_k = params.Common_lambda;
                else
                    lambda_k = 400;
                end
                pkt_size = 1000;
            end
        end
        
        % Cálculo M/G/1
        R_k_bps        = R_k_total_sum; 
        Stability_demand = lambda_k * pkt_size;

        if R_k_bps > Stability_demand
             % tempo médio de transmissão de um pacote
             L_tran  = pkt_size / R_k_bps; 
             
             % tempo médio de fila (M/G/1 approx)
             L_queue = (lambda_k * pkt_size^2) / ...
                       (2 * R_k_bps * (R_k_bps - Stability_demand));
             
             if is_D2D
                 % D2D: link direto (sem backhaul)
                 L_k_vector(k) = L_prop + L_tran + L_queue + L_proc;
             else
                 % CUE: via rede (com backhaul)
                 L_k_vector(k) = L_sb + L_prop + L_tran + L_queue;
             end
        else
             % Fila instável
             L_k_vector(k) = 1e12;
        end
    end
    
    R_total = sum(R_k_vector);
end
