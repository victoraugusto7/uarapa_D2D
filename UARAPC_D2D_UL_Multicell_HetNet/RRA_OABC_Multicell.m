function [bestAssign, bestOF, convFit, convObj] = RRA_OABC_Multicell(problem, V_viabilidade, params)
% Retorna convFit e convObj

    %% 1) Extrai Parâmetros
    nRBs       = problem.nRBs; 
    maxIter    = params.maxIter;
    colonySize = params.colonySize;
    limit      = params.limit;
    NL0        = params.NLsize;

    %% 2) Inicialização + OBL
    Solutions = zeros(colonySize, 3, nRBs);
    Fitness   = zeros(colonySize, 1);
    RealRate  = zeros(colonySize, 1);
    
    % População Normal
    for i = 1:colonySize
        X = generate_valid_solution(problem, V_viabilidade, params);
        Solutions(i, :, :) = X;
        [Fitness(i), RealRate(i)] = evaluateSolution_Multicell(X, V_viabilidade, problem, params);
    end

    % População Oposta
    Opp = zeros(colonySize, 3, nRBs);
    oppFit = zeros(colonySize, 1);
    oppReal = zeros(colonySize, 1);
    
    for i = 1:colonySize
        X_opp = generate_valid_solution(problem, V_viabilidade, params);
        Opp(i, :, :) = X_opp;
        [oppFit(i), oppReal(i)] = evaluateSolution_Multicell(X_opp, V_viabilidade, problem, params);
    end

    % Seleção
    allSol  = cat(1, Solutions, Opp); 
    allFit  = [Fitness; oppFit];
    allReal = [RealRate; oppReal];
    
    [~, ord] = sort(allFit, 'descend');
    keep     = ord(1:colonySize);
    
    Solutions = allSol(keep, :, :);
    Fitness   = allFit(keep);
    RealRate  = allReal(keep); 
    
    Trial     = zeros(colonySize, 1);
    convFit   = zeros(maxIter, 1);
    convObj   = zeros(maxIter, 1);

    [bestOF, idxMax] = max(Fitness);
    bestAssign       = squeeze(Solutions(idxMax, :, :));
    bestReal         = RealRate(idxMax);

    %% 3) Loop Principal OABC
    for iter = 1:maxIter
        NL_curr = max(1, round(NL0 * (1 - iter/maxIter)));

        % 3.1 Empregadas
        for i = 1:colonySize
            currSol = squeeze(Solutions(i, :, :));
            newSol = bestNeighbor_Multicell(currSol, NL_curr, V_viabilidade, problem, params);
            [newFit, newReal] = evaluateSolution_Multicell(newSol, V_viabilidade, problem, params);
            
            if newFit > Fitness(i)
                Solutions(i, :, :) = newSol;
                Fitness(i)         = newFit;
                RealRate(i)        = newReal;
                Trial(i)           = 0;
            else
                Trial(i) = Trial(i) + 1;
            end
        end

        % 3.2 Observadoras
        [~, rank] = sort(Fitness);
        prob = zeros(size(Fitness));
        prob(rank) = (1:colonySize) / sum(1:colonySize);
        
        i = 1; t = 0;
        while t < colonySize
            if rand < prob(i)
                currSol = squeeze(Solutions(i, :, :));
                newSol = bestNeighbor_Multicell(currSol, NL_curr, V_viabilidade, problem, params);
                [newFit, newReal] = evaluateSolution_Multicell(newSol, V_viabilidade, problem, params);
                
                if newFit > Fitness(i)
                    Solutions(i, :, :) = newSol;
                    Fitness(i)         = newFit;
                    RealRate(i)        = newReal;
                    Trial(i)           = 0;
                else
                    Trial(i) = Trial(i) + 1;
                end
                t = t + 1;
            end
            i = mod(i, colonySize) + 1;
        end

        % 3.3 Exploradoras
        for i = 1:colonySize
            if Trial(i) > limit
                X = generate_valid_solution(problem, V_viabilidade, params);
                Solutions(i, :, :) = X;
                [Fitness(i), RealRate(i)] = evaluateSolution_Multicell(X, V_viabilidade, problem, params);
                Trial(i) = 0;
            end
        end

        % 3.4 Atualiza Global
        [fmax, idxMax] = max(Fitness);
        if fmax > bestOF
            bestOF     = fmax;
            bestAssign = squeeze(Solutions(idxMax, :, :));
            bestReal   = RealRate(idxMax);
        end
        
        % 3.5 Intensificação Leve
        for k = 1:2
            cand = bestNeighbor_Multicell(bestAssign, 1, V_viabilidade, problem, params);
            [val, valReal] = evaluateSolution_Multicell(cand, V_viabilidade, problem, params);
            if val > bestOF
                bestOF     = val;
                bestReal   = valReal;
                bestAssign = cand;
            end
        end

        convFit(iter) = bestOF;
        convObj(iter) = bestReal;
    end
end