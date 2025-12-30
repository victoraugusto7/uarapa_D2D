function X = generate_valid_solution(problem, V_viabilidade, ~)

    nCUEs = problem.nCUE;
    nD2Ds = problem.nD2D;
    nRBs  = problem.nRBs;
    X = zeros(3, nRBs);

    % --- N_PRB aproximado para eMBB (a partir do modelo) ---
    gamma_lin   = 10^(problem.eMBB_Req_SINR_dB / 10);
    R_per_PRB   = problem.bwRB * log2(1 + gamma_lin);   % [bps]
    Nreq_eMBB   = ceil(problem.eMBB_Req_Rate / R_per_PRB);

    %% 1. CUEs
    cue_list = randperm(nCUEs);

    for u = cue_list
        b = problem.assoc_CUE(u);
        if b == 0, continue; end

        PRBs_licenciados = problem.PRBs_available_local{b};
        viability_row    = V_viabilidade(u, :);

        candidates = find(viability_row == 1 & X(1,:) == 0);
        candidates = intersect(candidates, PRBs_licenciados);
        if isempty(candidates), continue; end

        is_embb = problem.is_eMBB_CUE(u);
        if is_embb
            % tenta alocar Nreq_eMBB PRBs, limitado pelos candidatos
            num_to_alloc = min(length(candidates), Nreq_eMBB);
        else
            % CUE comum: 1 PRB (best effort)
            num_to_alloc = 1;
        end

        prbs_selected = candidates(randperm(length(candidates), num_to_alloc));
        X(1, prbs_selected) = u;
    end

    %% 2. D2Ds (igual ao seu código atual)
    d2d_list = randperm(nD2Ds);

    for k = 1:length(d2d_list)
        j_idx = d2d_list(k);
        u_d2d = nCUEs + j_idx;

        b_d2d = problem.assoc_D2D(j_idx);
        if b_d2d == 0, continue; end

        PRBs_licenciados = problem.PRBs_available_local{b_d2d};
        viability_d2d    = V_viabilidade(u_d2d, :);

        needed = 1;  % 1 PRB por D2D URLLC (por enquanto)

        potential_prbs = find(X(1,:) > 0 & viability_d2d == 1);
        potential_prbs = intersect(potential_prbs, PRBs_licenciados);
        if isempty(potential_prbs), continue; end

        potential_prbs = potential_prbs(randperm(length(potential_prbs)));

        count = 0;
        for r = potential_prbs
            if X(2,r) == 0
                X(2,r) = u_d2d;
                count = count + 1;
            elseif X(3,r) == 0
                X(3,r) = u_d2d;
                count = count + 1;
            end

            if count >= needed
                break;
            end
        end
    end
end
