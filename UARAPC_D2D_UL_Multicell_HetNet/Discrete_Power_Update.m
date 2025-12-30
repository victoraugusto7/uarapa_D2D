function P_new = Discrete_Power_Update(X, P_old, problem, params)
% DISCRETE_POWER_UPDATE
%   Atualiza o vetor de potências em NÍVEIS DISCRETOS, dado X.
%   Estratégia:
%   - 4 níveis para CUEs:  {0.25, 0.5, 0.75, 1} * Pmax_CUE
%   - 4 níveis para D2Ds:  {0.25, 0.5, 0.75, 1} * Pmax_D2D
%   - Aumenta / reduz nível com base nas métricas (R, BER, Lat).
%
% Inputs:
%   X       : matriz [3 x nRBs]
%   P_old   : vetor [nUsers x 1] (Watts) (pode ser vazio -> inicializa)
%   problem : struct com nCUE, nD2D, flags de eMBB/URLLC, etc.
%   params  : struct com txPowerCUE, txPowerD2D, Common_Req_*, etc.
%
% Output:
%   P_new   : vetor [nUsers x 1] atualizado (Watts)

    nCUE   = problem.nCUE;
    nD2D   = problem.nD2D;
    nUsers = nCUE + nD2D;

    % Se P_old vier vazio ou com dimensão errada, inicializa
    if nargin < 2 || numel(P_old) ~= nUsers
        P_old = init_power_vector(problem, params);
    end
    P_old = P_old(:); % garante coluna

    % dBm -> Watts (pico)
    Pmax_CUE  = 10^((params.txPowerCUE  - 30)/10);
    Pmax_D2D  = 10^((params.txPowerD2D  - 30)/10);

    levels_CUE = Pmax_CUE * [0.25 0.5 0.75 1.0];
    levels_D2D = Pmax_D2D * [0.25 0.5 0.75 1.0];

    % Mapa de alocação (pra saber quem usa PRB)
    [User_to_PRBs_Map, ~] = map_allocation_from_X(X, problem);

    % Métricas atuais com P_old
    [~, R_vec, BER_vec, L_vec] = ...
        calculate_qos_metrics(X, P_old, problem, params);

    P_new = P_old;

    % Requisitos básicos
    eMBB_R_req  = problem.eMBB_Req_Rate;
    eMBB_BERreq = problem.eMBB_Req_BER;
    eMBB_L_req  = problem.eMBB_Req_Lat;

    % Para Common (se existir)
    if isfield(params, 'Common_Req_Rate')
        Common_R_req  = params.Common_Req_Rate;
        Common_BERreq = params.Common_Req_BER;
        Common_L_req  = params.Common_Req_Lat;
    else
        Common_R_req  = 0.5e6;
        Common_BERreq = 1e-3;
        Common_L_req  = 100e-3;
    end

    % Para URLLC
    if isfield(problem, 'URLLC_Req_Rate')
        URLLC_R_req  = problem.URLLC_Req_Rate;
        URLLC_BERreq = problem.URLLC_Req_BER;
        URLLC_L_req  = problem.URLLC_Req_Lat;
    else
        URLLC_R_req  = 1e6;
        URLLC_BERreq = 1e-6;
        URLLC_L_req  = 20e-3;
    end

    % Flags de tipo de usuário
    is_eMBB_CUE   = problem.is_eMBB_CUE(:);
    is_URLLC_D2D  = false(nD2D,1);
    if isfield(problem,'is_URLLC_D2D')
        is_URLLC_D2D = problem.is_URLLC_D2D(:);
    end

    % Loop em todos usuários
    for k = 1:nUsers

        PRBs_k = User_to_PRBs_Map{k};

        % Se não tem PRB alocado: zera potência
        if isempty(PRBs_k)
            P_new(k) = 0;
            continue;
        end

        % Decide se é CUE ou D2D
        if k <= nCUE
            % ---- CUE ----
            isCritical = is_eMBB_CUE(k); % eMBB -> crítico

            R_req  = eMBB_R_req;
            BERreq = eMBB_BERreq;
            L_req  = eMBB_L_req;

            % Se não for eMBB, usa requisitos Common (menos rígidos)
            if ~isCritical
                R_req  = Common_R_req;
                BERreq = Common_BERreq;
                L_req  = Common_L_req;
            end

            levels = levels_CUE;

        else
            % ---- D2D ----
            d_idx = k - nCUE;
            isCritical = is_URLLC_D2D(d_idx);

            if isCritical
                R_req  = URLLC_R_req;
                BERreq = URLLC_BERreq;
                L_req  = URLLC_L_req;
            else
                % D2D best-effort (se existir no futuro)
                R_req  = Common_R_req;
                BERreq = Common_BERreq;
                L_req  = Common_L_req;
            end

            levels = levels_D2D;
        end

        % Se por algum motivo R_req=0, não faz muito sentido ajustar
        if R_req <= 0
            continue;
        end

        Rk   = R_vec(k);
        BERk = BER_vec(k);
        Lk   = L_vec(k);

        % Encontra índice de nível discreto mais próximo
        [~, idxLevel] = min(abs(levels - P_old(k)));

        violating = (Rk < R_req) || (BERk > BERreq) || (Lk > L_req);

        if isCritical
            % Usuários críticos (eMBB CUE ou D2D URLLC)

            if violating
                % tenta AUMENTAR potência um nível (se ainda não está no máx)
                if idxLevel < numel(levels)
                    idxLevel = idxLevel + 1;
                end
            else
                % margem de folga? Se taxa bem acima e latência bem abaixo,
                % podemos TENTAR reduzir um nível
                highRate   = (Rk >= 1.4 * R_req);
                lowLatency = (Lk <= 0.6 * L_req);

                if highRate && lowLatency && idxLevel > 1
                    idxLevel = idxLevel - 1;
                end
            end

        else
            % Usuários best-effort (CUE comum / D2D não-URLLC)

            % Se taxa muito abaixo do requisito mínimo
            if Rk < 0.7 * R_req
                if idxLevel < numel(levels)
                    idxLevel = idxLevel + 1;
                end
            elseif Rk > 2.0 * R_req
                % Muita sobra -> reduz
                if idxLevel > 1
                    idxLevel = idxLevel - 1;
                end
            end
        end

        % Atualiza potência no nível escolhido
        P_new(k) = levels(idxLevel);
    end
end
