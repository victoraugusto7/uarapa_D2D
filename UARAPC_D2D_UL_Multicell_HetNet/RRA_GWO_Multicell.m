function [bestAssign, bestOF, convFit, convObj] = RRA_GWO_Multicell(problem, V_viabilidade, params)
% Retorna convFit e convObj

    gwoAgents   = params.colonySize; 
    maxIter     = params.maxIter;
    nRBs        = problem.nRBs;
    
    Solutions = zeros(gwoAgents, 3, nRBs);
    Fitness   = zeros(gwoAgents, 1);
    RealRate  = zeros(gwoAgents, 1);
    
    for i = 1:gwoAgents
        X = generate_valid_solution(problem, V_viabilidade, params);
        Solutions(i, :, :) = X;
        [Fitness(i), RealRate(i)] = evaluateSolution_Multicell(X, V_viabilidade, problem, params);
    end

    [bestOF, idxMax] = max(Fitness);
    bestAssign       = squeeze(Solutions(idxMax, :, :));
    bestReal         = RealRate(idxMax);
    
    convFit = zeros(maxIter, 1);
    convObj = zeros(maxIter, 1);
    
    convFit(1) = bestOF;
    convObj(1) = bestReal; 
    
    for iter = 2:maxIter
        % Simulação de busca local para benchmark
        if rand < 0.3
            idx_rand = randi(gwoAgents);
            currSol = squeeze(Solutions(idx_rand, :, :));
            newSol = bestNeighbor_Multicell(currSol, 1, V_viabilidade, problem, params);
            [newFit, newReal] = evaluateSolution_Multicell(newSol, V_viabilidade, problem, params);
            
            if newFit > Fitness(idx_rand)
                Solutions(idx_rand, :, :) = newSol;
                Fitness(idx_rand) = newFit;
                RealRate(idx_rand) = newReal;
            end
        end

        [fmax, idxMax] = max(Fitness);
        if fmax > bestOF
            bestOF     = fmax;
            bestAssign = squeeze(Solutions(idxMax, :, :));
            bestReal   = RealRate(idxMax);
        end
        
        convFit(iter) = bestOF;
        convObj(iter) = bestReal;
    end
end