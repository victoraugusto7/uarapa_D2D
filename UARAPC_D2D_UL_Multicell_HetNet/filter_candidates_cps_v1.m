function V = filter_candidates_cps_v1(problem, params)
% filter_candidates_cps: Etapa 2 (CPS) - VERSÃO CORRIGIDA (Proteção Index 0)
%   Filtra (k, r) usando o teste P_min vs P_max.
%   Ignora usuários não associados (b_servidora = 0).

    fprintf('Executando Etapa 2: Seleção de Candidatos (CPS) [Finalizada]...\n');

    %% 1. Definição dos Limiares e Parâmetros
    XI_MIN_COMMON_DB = problem.Common_Req_SINR_dB; % 7.0 dB
    XI_MIN_EMBB_DB   = problem.eMBB_Req_SINR_dB;   % 10.0 dB
    K_max = params.Max_D2D_per_PRB;

    %% 2. Conversão para Escala Linear (Watts)
    P_max_c_watts = 10^((params.txPowerCUE - 30) / 10);
    P_max_d_watts = 10^((params.txPowerD2D - 30) / 10);
    
    noise_density_dbW_Hz = params.physParams.noiseDensity - 30;
    noise_watts = 10^(noise_density_dbW_Hz / 10) * problem.bwRB;

    xi_min_common_lin = 10^(XI_MIN_COMMON_DB / 10);
    xi_min_embb_lin   = 10^(XI_MIN_EMBB_DB / 10);
    
    % --- Interferência Média (Se quiser usar depois) ---
    % I_CCI_media_watts = 10^((params.Interf_CCI_media_dBm - 30) / 10); 

    %% 3. Inicialização
    nCUEs = problem.nCUE;
    nD2Ds = problem.nD2D;
    nTotalUsers = nCUEs + nD2Ds;
    V = zeros(nTotalUsers, problem.nRBs); % Matriz de Viabilidade (Saída)

    %% 4. Loop Principal: Testar (k)
    
    for k = 1:nTotalUsers % Iterar por todos os usuários
        
        % --- A. Identifica o usuário 'k' e Pega Parâmetros Úteis ---
        if k <= nCUEs
            % Usuário k é um CUE (i)
            i = k;
            b_servidora = problem.assoc_CUE(i);
            
            % --- PROTEÇÃO CONTRA USER NÃO ASSOCIADO ---
            if b_servidora == 0
                continue; % Pula este usuário (Viabilidade = 0)
            end
            % ------------------------------------------
            % --- IMUNIDADE: Força passagem de todos eMBB para a Fase 3 ---
            if problem.is_eMBB_CUE(i)
                V(k, :) = 1; % Aprova em todos os PRBs
                continue;    % Pula o teste de P_min para este usuário
            end
            G_k = problem.gain_CUE_to_MBSs(i, b_servidora); % Ganho útil
            P_max_k = P_max_c_watts;
            
            % Define SINR min
            if problem.is_eMBB_CUE(i)
                xi_min_k = xi_min_embb_lin;
            else
                xi_min_k = xi_min_common_lin;
            end
            isCUE = true;
            
        else
            % Usuário k é um D2D (j)
            j = k - nCUEs;
            b_servidora = problem.assoc_D2D(j); % Célula de ancoragem
            
            % --- PROTEÇÃO CONTRA USER NÃO ASSOCIADO ---
            if b_servidora == 0
                continue; % Pula este usuário
            end
            % ------------------------------------------

            G_k = problem.gain_D2D_pair(j); % Ganho útil
            P_max_k = P_max_d_watts;
            xi_min_k = xi_min_common_lin; % D2Ds são Comuns/URLLC (7dB)
            isCUE = false;
        end
        
        % 1. Filtro Rápido (Se Ganho útil for zero, reprova)
        if G_k == 0
            continue; 
        end
        
        % --- B. CÁLCULO I_INTRA (Pior Caso Local) ---
        I_intra = 0;
        CUEs_locais = find(problem.assoc_CUE == b_servidora);
        D2Ds_locais = find(problem.assoc_D2D == b_servidora);

        if isCUE % Vítima é CUE i (Interferência de D2Ds locais)
            if ~isempty(D2Ds_locais)
                gains_d2d_mbs = problem.gain_D2DTx_to_MBSs(D2Ds_locais, b_servidora);
                [~, idxs] = maxk(gains_d2d_mbs, min(K_max, length(D2Ds_locais)));
                I_intra = sum(P_max_d_watts * gains_d2d_mbs(idxs));
            end
        else % Vítima é D2D j (Interferência de CUEs + Outros D2Ds locais)
            if ~isempty(CUEs_locais)
                gains_cue_d2d = problem.gain_CUE_to_D2DRx(j, CUEs_locais);
                [max_gain, ~] = max(gains_cue_d2d);
                I_intra = P_max_c_watts * max_gain;
            end
            D2Ds_locais_outros = D2Ds_locais(D2Ds_locais ~= j);
            if ~isempty(D2Ds_locais_outros)
                gains_d2d_d2d = problem.gain_D2DTx_to_D2DRx(j, D2Ds_locais_outros);
                [~, idxs] = maxk(gains_d2d_d2d, min(K_max - 1, length(D2Ds_locais_outros)));
                I_intra = I_intra + sum(P_max_d_watts * gains_d2d_d2d(idxs));
            end
        end
        % --- FIM DO CÁLCULO I_intra ---

        % --- C. DECISÃO DE ADMISSÃO (P_min vs P_max) ---
        
        % I_total_teste = Ruído + I_Intra Pior Caso (CCI é zero no Reuso 3)
        I_total_teste = noise_watts + I_intra;
            
        % P_min = (SINR_min_linear * (Ruído + I_Total)) / Ganho_útil
        P_min_k = (xi_min_k * I_total_teste) / G_k;
        
        % A Decisão de Admissão
        if P_min_k <= P_max_k
            V(k, :) = 1; % APROVADO
        else
            V(k, :) = 0; % REPROVADO
        end
    end
    
    fprintf('Filtragem CPS (Versão Final Protegida) concluída. Matriz V criada.\n');
end