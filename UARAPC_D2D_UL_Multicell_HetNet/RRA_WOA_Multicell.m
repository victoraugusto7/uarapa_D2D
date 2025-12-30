function [bestAssign, bestOF, convFit, convObj] = RRA_WOA_Multicell(problem, V_viabilidade, params)
% RRA_WOA_Multicell: Algoritmo de Otimização das Baleias (Discreto/Híbrido)
% Adaptado para Matriz de Alocação 3xN e Controle de Potência.

    %% 1. Parâmetros
    nRBs       = problem.nRBs; 
    maxIter    = params.maxIter;
    popSize    = params.colonySize;
    
    %% 2. Inicialização
    Solutions = zeros(popSize, 3, nRBs); 
    Fitness   = zeros(popSize, 1);
    RealRate  = zeros(popSize, 1);
    
    % Gera população inicial válida
    for i = 1:popSize
        X = generate_valid_solution(problem, V_viabilidade, params);
        Solutions(i, :, :) = X; 
        [Fitness(i), RealRate(i)] = evaluateSolution_Multicell(X, V_viabilidade, problem, params);
    end

    % Encontra a melhor baleia (Líder)
    [bestOF, idxMax] = max(Fitness);
    bestAssign       = squeeze(Solutions(idxMax, :, :));
    bestReal         = RealRate(idxMax);

    convFit = zeros(maxIter, 1);
    convObj = zeros(maxIter, 1);

    %% 3. Loop Principal do WOA
    for t = 1:maxIter
        
        % a: decresce linearmente de 2 para 0
        a = 2 - t * (2 / maxIter); 
        
        % a2: decresce de -1 para -2 (para a espiral)
        a2 = -1 + t * ((-1)/maxIter);

        for i = 1:popSize
            
            r1 = rand(); 
            r2 = rand();
            
            A = 2 * a * r1 - a;  % Vetor A
            C = 2 * r2;          % Vetor C
            
            b = 1;               % Parâmetro da espiral logarítmica
            l = (a2-1)*rand + 1; % Número aleatório em [-1, 1]
            p = rand();          % Probabilidade de escolha de comportamento

            currentSol = squeeze(Solutions(i, :, :));
            newSol     = currentSol; % Candidato
            
            if p < 0.5
                if abs(A) < 1
                    % --- [1] ENCIRCLING PREY (Explotação) ---
                    % Aproxima-se da Melhor Solução (Crossover Discreto)
                    % Copia uma porcentagem de colunas da Melhor Solução para o Agente
                    mix_rate = 0.5; % Taxa de mistura
                    mask = rand(1, nRBs) < mix_rate;
                    newSol(:, mask) = bestAssign(:, mask);
                    
                else
                    % --- [2] SEARCH FOR PREY (Exploração) ---
                    % Aproxima-se de uma Baleia Aleatória
                    randIdx = randi(popSize);
                    randSol = squeeze(Solutions(randIdx, :, :));
                    
                    mix_rate = 0.3; % Mistura leve
                    mask = rand(1, nRBs) < mix_rate;
                    newSol(:, mask) = randSol(:, mask);
                end
            else
                % --- [3] SPIRAL BUBBLE-NET (Ataque Refinado) ---
                % Simula o movimento espiral fazendo uma busca local na Melhor Solução
                % Quanto mais perto do fim (l), menor o raio da busca (NL)
                
                NL_curr = max(1, round(params.NLsize * abs(l)));
                % A baleia segue a líder, mas tenta melhorar (Busca Local)
                newSol = bestNeighbor_Multicell(bestAssign, NL_curr, V_viabilidade, problem, params);
            end
            
            % --- Avaliação e Seleção ---
            % O WOA aceita se melhorar (Greedy)
            [fNew, rNew] = evaluateSolution_Multicell(newSol, V_viabilidade, problem, params);
            
            if fNew > Fitness(i)
                Solutions(i, :, :) = newSol;
                Fitness(i)         = fNew;
                RealRate(i)        = rNew;
            end
            
            % Atualiza o Líder Global se necessário
            if fNew > bestOF
                bestOF     = fNew;
                bestAssign = newSol;
                bestReal   = rNew;
            end
        end
        
        % Registra convergência
        convFit(t) = bestOF;   % Fitness (com penalidade)
        convObj(t) = bestReal; % Throughput Real (Mbps)
        
    end
end