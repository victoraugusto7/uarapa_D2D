function V = filter_candidates_cps_threshold(problem, params)
% filter_candidates_cps_threshold: Etapa 2 (CPS) - Versão com Limiares (Thresholds)
%   Cria uma matriz de viabilidade V(k, r) (Usuário x PRB)
%   baseada em limiares de Sinal Útil e Interferência Gerada.
%
%   V(k,r) = 1 se o usuário k PODE usar o PRB r (Aprovado)
%   V(k,r) = 0 se o usuário k NÃO PODE usar o PRB r (Reprovado)
%
%   Nota: Este filtro aprova um *usuário*, não um par (usuário, PRB).
%   Um usuário aprovado (V(k, :) = 1) pode ser testado pela MHA em qualquer PRB.

    fprintf('Executando Etapa 2: Seleção de Candidatos (CPS) por Limiar...\n');

    %% 1. Definição dos Limiares (Valores de exemplo, em dBm)
    % (Estes valores virão de 'params')
    
    % --- Para Usuários eMBB (Rigoroso) ---
    LIMIAR_SINAL_UTIL_MIN_EMBB_DBM = params.eMBB_Sinal_min_dBm;   % Ex: -90 dBm
    LIMIAR_INTERF_MAX_EMBB_DBM   = params.eMBB_Interf_max_dBm; % Ex: -115 dBm

    % --- Para Usuários Comuns (Flexível) ---
    LIMIAR_SINAL_UTIL_MIN_COMUM_DBM = params.Common_Sinal_min_dBm;   % Ex: -100 dBm
    LIMIAR_INTERF_MAX_COMUM_DBM   = params.Common_Interf_max_dBm; % Ex: -105 dBm

    %% 2. Ganhos de Canal (em dB)
    % Convertendo os ganhos lineares (calculados na Etapa 0) para dB
    % G_dB = 10 * log10(G_linear)
    
    gain_CUE_MBSs_dB = 10 * log10(problem.gain_CUE_to_MBSs);
    gain_D2D_pair_dB = 10 * log10(problem.gain_D2D_pair);
    gain_CUE_D2DRx_dB = 10 * log10(problem.gain_CUE_to_D2DRx);
    gain_D2DTx_MBSs_dB = 10 * log10(problem.gain_D2DTx_to_MBSs);
    gain_D2DTx_D2DRx_dB = 10 * log10(problem.gain_D2DTx_to_D2DRx);
    
    % Define ganhos de -inf (links de 0) para um valor muito pequeno
    gain_D2DTx_D2DRx_dB(isinf(gain_D2DTx_D2DRx_dB)) = -300; 

    %% 3. Inicialização
    nCUEs = problem.nCUE;
    nD2Ds = problem.nD2D;
    nTotalUsers = nCUEs + nD2Ds;
    nRBs = problem.nRBs;
    V = zeros(nTotalUsers, nRBs); % Matriz de Viabilidade (Saída)

    %% 4. Filtragem dos CUEs
    for i = 1:nCUEs
        b = problem.assoc_CUE(i); % MBS servidora do CUE i
        
        % Define os limiares corretos (eMBB ou Comum)
        if problem.is_eMBB_CUE(i)
            limiar_sinal = LIMIAR_SINAL_UTIL_MIN_EMBB_DBM;
            limiar_interf = LIMIAR_INTERF_MAX_EMBB_DBM;
        else
            limiar_sinal = LIMIAR_SINAL_UTIL_MIN_COMUM_DBM;
            limiar_interf = LIMIAR_INTERF_MAX_COMUM_DBM;
        end

        % --- Teste 1: Sinal Útil (do CUE i para sua MBS b) ---
        S_util_dBm = params.txPowerCUE + gain_CUE_MBSs_dB(i, b);
        
        if S_util_dBm < limiar_sinal
            continue; % REPROVADO (Sinal útil muito fraco)
        end
        
        % --- Teste 2: Interferência Gerada (do CUE i para D2Ds) ---
        % Encontra a interferência MÁXIMA que este CUE causa em QUALQUER D2D-Rx
        interf_gerada_dBm = params.txPowerCUE + max(gain_CUE_D2DRx_dB(:, i));
        
        if interf_gerada_dBm > limiar_interf
            continue; % REPROVADO (Interfere demais nos D2Ds)
        end
        
        % --- APROVADO ---
        V(i, :) = 1; % Aprovado para todos os PRBs
    end
    
    %% 5. Filtragem dos D2Ds
    idx_offset = problem.nCUE; % Offset para o índice global
    
    for j = 1:nD2Ds
        % (Como D2Ds são 100% Comuns, usamos limiares Comuns)
        limiar_sinal = LIMIAR_SINAL_UTIL_MIN_COMUM_DBM;
        limiar_interf = LIMIAR_INTERF_MAX_COMUM_DBM;
        
        % --- Teste 1: Sinal Útil (interno do par D2D j) ---
        S_util_dBm = params.txPowerD2D + gain_D2D_pair_dB(j);
        
        if S_util_dBm < limiar_sinal
            continue; % REPROVADO (Sinal D2D interno muito fraco)
        end
        
        % --- Teste 2: Interferência Gerada (D2D j -> Vítimas) ---
        
        % Vítima 1: Interferência máxima que o D2D-Tx j causa em QUALQUER MBS
        interf_max_para_MBSs = params.txPowerD2D + max(gain_D2DTx_MBSs_dB(j, :));
        
        % Vítima 2: Interferência máxima que o D2D-Tx j causa em QUALQUER D2D-Rx
        interf_max_para_D2Ds = params.txPowerD2D + max(gain_D2DTx_D2DRx_dB(:, j));
        
        interf_gerada_max_dBm = max(interf_max_para_MBSs, interf_max_para_D2Ds);
        
        if interf_gerada_max_dBm > limiar_interf
            continue; % REPROVADO (Interfere demais)
        end
        
        % --- APROVADO ---
        V(idx_offset + j, :) = 1;
    end
    
    fprintf('Filtragem CPS (Limiares) concluída. Matriz V criada.\n');
end