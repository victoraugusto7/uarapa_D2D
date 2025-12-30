function V = filter_candidates_cps_v2(problem, params)
% filter_candidates_cps: Etapa 2 (CPS) - VERSÃO RELAXADA (High Approval)
%   Objetivo: Aprovar o máximo de usuários viáveis (focando em SNR).
%   Deixa a gestão da interferência (SINR) para a Fase 3 (Otimização).

    fprintf('Executando Etapa 2: Seleção de Candidatos (CPS) [Modo Relaxado]...\n');

    %% 1. Parâmetros
    % Usamos um SINR alvo mais baixo para admissão (apenas para passar no portão)
    % Se o alvo real é 10dB, aqui aceitamos quem tem potencial para 3dB.
    XI_MIN_COMMON_DB = problem.Common_Req_SINR_dB - 3; 
    XI_MIN_EMBB_DB   = problem.eMBB_Req_SINR_dB - 3;   
    
    % Potências em Watts
    P_max_c_watts = 10^((params.txPowerCUE - 30) / 10);
    P_max_d_watts = 10^((params.txPowerD2D - 30) / 10);
    
    % Ruído Térmico
    noise_density_dbW_Hz = params.physParams.noiseDensity - 30;
    noise_watts = 10^(noise_density_dbW_Hz / 10) * problem.bwRB;

    % Limiares Lineares
    xi_min_common_lin = 10^(XI_MIN_COMMON_DB / 10);
    xi_min_embb_lin   = 10^(XI_MIN_EMBB_DB / 10);
    
    % Interferência de Borda (Inter-cell) - Valor fixo médio/baixo
    if isfield(params, 'Interf_CCI_media_dBm')
        I_CCI_watts = 10^((params.Interf_CCI_media_dBm - 30) / 10);
    else
        I_CCI_watts = 0; % Ou muito baixo
    end

    %% 2. Inicialização
    nTotalUsers = problem.nCUE + problem.nD2D;
    V = zeros(nTotalUsers, problem.nRBs); 

    %% 3. Loop de Filtragem
    for k = 1:nTotalUsers
        
        % --- A. Identificação ---
        if k <= problem.nCUE
            i = k;
            b_serv = problem.assoc_CUE(i);
            if b_serv == 0, continue; end
            
            G_k = problem.gain_CUE_to_MBSs(i, b_serv);
            P_max_k = P_max_c_watts;
            
            % eMBB tem prioridade total (Imunidade) ou check leve?
            % Vamos fazer um check leve de SNR para não aprovar quem está "morto"
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
            
            G_k = problem.gain_D2D_pair(j);
            P_max_k = P_max_d_watts;
            xi_min_k = xi_min_common_lin; 
            is_embb = false;
        end
        
        if G_k == 0, continue; end
        
        % --- B. CÁLCULO DE INTERFERÊNCIA (AQUI ESTÁ A MUDANÇA) ---
        
        % ANTES (Severo): I_total = Ruído + CCI + I_Intra_Pior_Caso
        % AGORA (Relaxado): I_total = Ruído + CCI (Ignora Intra-cell)
        % Assumimos que o ABC vai alocar PRBs ortogonais ou controlar a potência
        % para mitigar a I_intra.
        
        I_total_teste = noise_watts + I_CCI_watts;
        
        % --- C. DECISÃO (Check de Potência) ---
        
        % P_req = (SINR_alvo * I_total) / Ganho_Canal
        P_req = (xi_min_k * I_total_teste) / G_k;
        
        % Se for eMBB, aplicamos "Imunidade Quase Total" 
        % (Aceita mesmo se P_req for um pouco alto, confiando no ganho de processamento)
        if is_embb
             if P_req <= (P_max_k * 2) % Margem de 3dB extra para eMBB
                 V(k, :) = 1; 
             end
        else
             % Para comuns/D2D, check normal de SNR
             if P_req <= P_max_k
                 V(k, :) = 1;
             end
        end
    end
    
    % Diagnóstico Rápido
    aprovados = sum(any(V, 2));
    fprintf('  => CPS Relaxado: %d/%d usuários aprovados para a Fase 3.\n', aprovados, nTotalUsers);
end