function [fitness, throughput_real] = evaluateSolution_Multicell_old(X, V_viabilidade, problem, params)
% Retorna [Fitness_com_Penalidade, Throughput_Real_Limpo]
% Agora com penalidade também para D2D URLLC.

    % 1. Validação Estrutural
    if ~validateAllocationStructure(X, V_viabilidade, problem, params)
        fitness         = -1e15; 
        throughput_real = 0;
        return; 
    end
    
    % 2. Otimização de Potência (GWO)
    P_star = optimal_power_control(X, problem, params);

    % 3. Cálculo de Métricas
    [R_total, R_k_vec, BER_k_vec, L_k_vec] = calculate_qos_metrics(X, P_star, problem, params);
    
    % Este é o valor que queremos no gráfico (Throughput Físico)
    throughput_real = R_total; 

    % 4. Cálculo de Penalidades (Para o Fitness)
    penalty_score = 0;
    P_base        = 1e9; 
    
    %% 4.a) CUEs eMBB (Prioridade Máxima)
    idx_eMBB = find(problem.is_eMBB_CUE);
    for k = idx_eMBB'
        % (i) Taxa
        if R_k_vec(k) < problem.eMBB_Req_Rate
            deficit       = (problem.eMBB_Req_Rate - R_k_vec(k));
            penalty_score = penalty_score + P_base + (deficit * 100); 
        end
        % (ii) BER
        if BER_k_vec(k) > problem.eMBB_Req_BER
            penalty_score = penalty_score + P_base;
        end
        % (iii) Latência
        if L_k_vec(k) > problem.eMBB_Req_Lat
            penalty_score = penalty_score + P_base;
        end
    end
    
    %% 4.b) D2D URLLC (Pares Agentes Críticos)
    if isfield(problem, 'is_URLLC_D2D')
        nCUEs          = problem.nCUE;
        idx_D2D_URLLC  = find(problem.is_URLLC_D2D);  % índices 1..nD2D
        
        for d = idx_D2D_URLLC'
            k = nCUEs + d;  % índice global do usuário D2D
            
            % (i) Taxa mínima URLLC
            if R_k_vec(k) < problem.URLLC_Req_Rate
                deficit       = (problem.URLLC_Req_Rate - R_k_vec(k));
                penalty_score = penalty_score + P_base + (deficit * 100);
            end
            
            % (ii) BER máxima URLLC
            if BER_k_vec(k) > problem.URLLC_Req_BER
                penalty_score = penalty_score + P_base;
            end
            
            % (iii) Latência máxima URLLC
            if L_k_vec(k) > problem.URLLC_Req_Lat
                penalty_score = penalty_score + P_base;
            end
        end
    end
    
    % 5. Fitness Final (Matemático)
    fitness = R_total - penalty_score;
end
