function [bestAssign, bestOF, convFit, convObj] = RRA_ABC_Multicell(problem, V_viabilidade, params)
% Retorna convFit (Fitness com penalidade) e convObj (Throughput Real)

    %% 1) Extrai Parâmetros
    nRBs       = problem.nRBs; 
    maxIter    = params.maxIter;
    colonySize = params.colonySize;
    employed   = colonySize;
    onlooker   = colonySize;
    limit      = params.limit;
    NL         = params.NLsize;
    
    %% 2) Inicializa Soluções
    Solutions = zeros(colonySize, 3, nRBs); 
    Fitness   = zeros(colonySize, 1);
    RealRate  = zeros(colonySize, 1);
    Trial     = zeros(colonySize, 1);
    
    for i = 1:colonySize
        X = generate_valid_solution(problem, V_viabilidade, params);
        Solutions(i, :, :) = X; 
        [Fitness(i), RealRate(i)] = evaluateSolution_Multicell(X, V_viabilidade, problem, params);
    end

    convFit = zeros(maxIter, 1); % Convergência do Fitness (Matemático)
    convObj = zeros(maxIter, 1); % Convergência do Objetivo (Throughput Real)
    
    [~, idxMax] = max(Fitness);
    bestAssign  = squeeze(Solutions(idxMax, :, :));
    bestOF      = Fitness(idxMax);   
    bestReal    = RealRate(idxMax);

    %% 3) Loop Principal ABC
    for iter = 1:maxIter
        
        % 3.1 Empregadas
        for i = 1:employed
            currSol = squeeze(Solutions(i, :, :));
            newSol = bestNeighbor_Multicell(currSol, NL, V_viabilidade, problem, params);
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
        while t < onlooker
            if rand < prob(i)
                currSol = squeeze(Solutions(i, :, :));
                newSol = bestNeighbor_Multicell(currSol, NL, V_viabilidade, problem, params);
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
            i = mod(i, employed) + 1;
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

        % 3.4 Atualiza Melhor Global
        [fmax, idxMax] = max(Fitness);
        if fmax > bestOF
            bestOF     = fmax;
            bestAssign = squeeze(Solutions(idxMax, :, :));
            bestReal   = RealRate(idxMax);
        end

        convFit(iter) = bestOF;   % Fitness (com penalidade)
        convObj(iter) = bestReal; % Objetivo (Throughput Real)
    end
end