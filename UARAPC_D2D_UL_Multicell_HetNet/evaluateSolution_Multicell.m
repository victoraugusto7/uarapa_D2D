function [U, R_total, P_lin, metrics] = evaluateSolution_Multicell(X, V_viabilidade, problem, params)
% =========================================================================
% Avaliação unificada para RA (X) + PC nested (otimiza P dado X).
% Retorna:
%   U       : fitness/utility (MAIOR = melhor)
%   R_total : throughput total (bps)
%   P_lin   : potências (W), dimensão (nUsers x 1)
%   metrics : struct com vetores QoS, PRBs por usuário e contagens atendidos
% =========================================================================

    % -------- 0) Normalizações defensivas
    if isempty(X)
        U = -Inf; R_total = 0; P_lin = []; metrics = struct(); return;
    end

    % Garantir dimensão esperada: (1+K) x nRBs
    K = params.Max_D2D_per_PRB;
    nRBs = problem.nRBs;

    if size(X,1) ~= (1+K) && size(X,2) == (1+K)
        X = X.'; % tenta corrigir transposto
    end

    if size(X,1) ~= (1+K) || size(X,2) ~= nRBs
        % Penaliza forte, mas não explode execução:
        U = -1e12; R_total = 0; P_lin = []; metrics = struct('flag','BAD_DIM');
        return;
    end

    % -------- 1) (Opcional) Checagem de viabilidade CPS
    % A ideia: V_viabilidade filtra candidatos na fase 2.
    % Se X já é gerado respeitando V, você pode passar [] sem problema.
    if ~isempty(V_viabilidade)
        if ~is_solution_compatible_with_V(X, V_viabilidade, problem, params)
            U = -1e12; R_total = 0; P_lin = []; metrics = struct('flag','VIOLATES_V');
            return;
        end
    end

    % -------- 2) PC nested (otimiza P dado X)
    P_lin = optimal_power_control(X, problem, params);

    % -------- 3) QoS por usuário
    [R_total, R_vec, BER_vec, L_vec] = calculate_qos_metrics(X, P_lin, problem, params);

    % -------- 4) PRBs por usuário (CUE na linha 1, D2D nas linhas 2..)
    prbCount = compute_prb_counts_from_X(X, problem, params);

    % -------- 5) Utility / fitness (soft)
    [U, sat] = compute_utility_soft(R_total, R_vec, BER_vec, L_vec, problem, params);

    % -------- 6) Estrutura de métricas (para relatório)
    metrics = struct();
    metrics.R_total_bps = R_total;
    metrics.R_vec_bps   = R_vec;
    metrics.BER_vec     = BER_vec;
    metrics.L_vec_s     = L_vec;
    metrics.P_lin_W     = P_lin;
    metrics.P_dBm       = 10*log10(max(P_lin, eps)) + 30;
    metrics.prbCount    = prbCount;

    metrics.sat_eMBB_frac  = sat.sat_e;
    metrics.sat_URLLC_frac = sat.sat_u;
    metrics.n_eMBB_ok       = sat.n_e_ok;
    metrics.n_URLLC_ok      = sat.n_u_ok;
end

% ========================= helpers =========================

function ok = is_solution_compatible_with_V(X, V, problem, ~)
    % Regra simples: para cada RB, o CUE escolhido precisa ser candidato e
    % os D2D escolhidos (se != 0) também.
    ok = true;
    nCUE = problem.nCUE;
    nD2D = problem.nD2D;
    nRBs = problem.nRBs;

    for r = 1:nRBs
        cue = X(1,r);
        if cue < 1 || cue > nCUE || ~V(cue,r)
            ok = false; return;
        end
        for k = 2:size(X,1)
            d = X(k,r);
            if d == 0, continue; end
            if d < 1 || d > nD2D || ~V(nCUE + d, r)
                ok = false; return;
            end
        end
    end
end

function prbCount = compute_prb_counts_from_X(X, problem, ~)
    nCUE = problem.nCUE;
    nD2D = problem.nD2D;

    prbCount = zeros(nCUE + nD2D, 1);

    % CUEs (linha 1): índices globais 1..nCUE
    cueIdx = X(1,:);
    for i = 1:nCUE
        prbCount(i) = sum(cueIdx == i);
    end

    % D2D (linhas 2..): índices locais 1..nD2D, armazenados como 0 ou j
    d2dMat = X(2:end,:);
    for j = 1:nD2D
        prbCount(nCUE + j) = sum(d2dMat(:) == j);
    end
end

function [U, sat] = compute_utility_soft(R_total, R_vec, BER_vec, L_vec, problem, params)
    nCUE = problem.nCUE;

    idx_e = find(problem.is_eMBB_CUE(:));
    idx_u = nCUE + find(problem.is_URLLC_D2D(:));

    % Satisfação (máscaras)
    satMask_e = (R_vec(idx_e) >= params.eMBB_Req_Rate) & ...
                (BER_vec(idx_e) <= params.eMBB_Req_BER) & ...
                (L_vec(idx_e) <= params.eMBB_Req_Lat);

    satMask_u = (R_vec(idx_u) >= params.URLLC_Req_Rate) & ...
                (BER_vec(idx_u) <= params.URLLC_Req_BER) & ...
                (L_vec(idx_u) <= params.URLLC_Req_Lat);

    sat.sat_e   = mean(satMask_e);
    sat.sat_u   = mean(satMask_u);
    sat.n_e_ok  = sum(satMask_e);
    sat.n_u_ok  = sum(satMask_u);

    % Penalidades normalizadas (estáveis)
    rel_def    = @(x, req) max(0, (req - x) ./ max(req, eps));
    rel_excess = @(x, req) max(0, (x - req) ./ max(req, eps));
    log_excess = @(x, req) max(0, log10(max(x, eps) ./ max(req, eps)));

    pen_e = mean(rel_def(R_vec(idx_e),   params.eMBB_Req_Rate)) + ...
            mean(log_excess(BER_vec(idx_e), params.eMBB_Req_BER)) + ...
            mean(rel_excess(L_vec(idx_e), params.eMBB_Req_Lat));

    pen_u = mean(rel_def(R_vec(idx_u),   params.URLLC_Req_Rate)) + ...
            mean(log_excess(BER_vec(idx_u), params.URLLC_Req_BER)) + ...
            mean(rel_excess(L_vec(idx_u), params.URLLC_Req_Lat));

    % Pesos (preferencialmente parametrizáveis)
    if ~isfield(params,'w_thr'), params.w_thr = 1.0; end           % por Mbps
    if ~isfield(params,'w_sat_e'), params.w_sat_e = 4000; end
    if ~isfield(params,'w_sat_u'), params.w_sat_u = 4000; end
    if ~isfield(params,'w_pen'), params.w_pen = 1500; end

    R_total_Mbps = R_total / 1e6;

    U = params.w_thr * R_total_Mbps + ...
        params.w_sat_e * sat.sat_e + ...
        params.w_sat_u * sat.sat_u - ...
        params.w_pen * (0.5*(pen_e + pen_u));
end
