function V = filter_candidates_cps(problem, params)
% filter_candidates_cps: Etapa 2 (CPS) - Modo Relaxado + Respeito ao Reuse-3
%   Objetivo: Aprovar o máximo de usuários viáveis (focando em SNR),
%   mas restringindo a viabilidade aos PRBs do conjunto local de cada célula.

    fprintf('Executando Etapa 2: Seleção de Candidatos (CPS) [Relaxado + Reuse-3]...\n');
    is_urllc = false;

    %% 1. Parâmetros
    % SINR alvo relaxado (portão de admissão)
    XI_MIN_COMMON_DB = problem.Common_Req_SINR_dB - 3; 
    XI_MIN_EMBB_DB   = problem.eMBB_Req_SINR_dB   - 3;   

    % URLLC (D2D) - alvo relaxado também
    XI_MIN_URLLC_DB  = problem.URLLC_Req_SINR_dB - 3;

    
    % Potências máximas em Watts
    P_max_c_watts = 10^((params.txPowerCUE - 30) / 10);
    P_max_d_watts = 10^((params.txPowerD2D - 30) / 10);
    
    % Ruído térmico por PRB
    noise_density_dbW_Hz = params.physParams.noiseDensity - 30;
    noise_watts = 10^(noise_density_dbW_Hz / 10) * problem.bwRB;

    % Limiares lineares de SINR
    xi_min_common_lin = 10^(XI_MIN_COMMON_DB / 10);
    xi_min_embb_lin   = 10^(XI_MIN_EMBB_DB   / 10);

    xi_min_common_lin = 10^(XI_MIN_COMMON_DB / 10);
    xi_min_embb_lin   = 10^(XI_MIN_EMBB_DB   / 10);
    xi_min_urllc_lin  = 10^(XI_MIN_URLLC_DB  / 10);

    
    % Interferência média inter-cell (CCI)
    if isfield(params, 'Interf_CCI_media_dBm')
        I_CCI_watts = 10^((params.Interf_CCI_media_dBm - 30) / 10);
    else
        I_CCI_watts = 0;
    end

    % Interferência total de teste (sem intra-cell)
    I_total_teste = noise_watts + I_CCI_watts;

    %% 2. Inicialização
    nTotalUsers = problem.nCUE + problem.nD2D;
    V = zeros(nTotalUsers, problem.nRBs); 

    %% 3. Loop de Filtragem (por usuário)
    for k = 1:nTotalUsers
        
        % --- A. Identificação do usuário, célula servidora e ganho útil ---
        if k <= problem.nCUE
            i = k;
            b_serv = problem.assoc_CUE(i);
            if b_serv == 0, continue; end   % CUE fora de cobertura

            G_k    = problem.gain_CUE_to_MBSs(i, b_serv);
            P_max_k = P_max_c_watts;
            
            if problem.is_eMBB_CUE(i)
                xi_min_k = xi_min_embb_lin;
                is_embb = true;
            else
                xi_min_k = xi_min_common_lin;
                is_embb = false;
            end
        else
    
            j = k - problem.nCUE;
            b_serv = problem.assoc_D2D(j);
            if b_serv == 0, continue; end
            
            G_k    = problem.gain_D2D_pair(j);
            P_max_k = P_max_d_watts;
        
            % D2D pode ser URLLC ou "Comum"
            if isfield(problem,'is_URLLC_D2D') && problem.is_URLLC_D2D(j)
                xi_min_k = xi_min_urllc_lin;
                is_embb  = false;      % não é eMBB, mas é crítico (tratado no fitness)
                is_urllc = true;
            else
                xi_min_k = xi_min_common_lin;
                is_embb  = false;
                is_urllc = false;
            end
        end

        
        if G_k == 0, continue; end
        
        % --- B. Potência requerida para atingir SINR alvo relaxado ---
        P_req = (xi_min_k * I_total_teste) / G_k;
        
        % --- C. Decisão de admissão ---
        % PRBs locais dessa célula (reuse-3)
        rb_local = problem.PRBs_available_local{b_serv};
        
        if is_embb
            % eMBB: imunidade quase total (até ~3 dB acima do P_max)
            if P_req <= (P_max_k * 2)
         V(k, rb_local) = 1;
        end
        elseif exist('is_urllc','var') && is_urllc
            % D2D URLLC: um pouco mais rígido que Comum (sem margem extra)
        if P_req <= P_max_k
            V(k, rb_local) = 1;
        end 
    else
        % Comuns/D2D não-URLLC: filtro normal
        if P_req <= P_max_k
        V(k, rb_local) = 1;
        end
end


    end
    
    % Diagnóstico
    aprovados = sum(any(V, 2));
    fprintf('  => CPS Relaxado: %d/%d usuários aprovados para a Fase 3.\n', aprovados, nTotalUsers);
end
