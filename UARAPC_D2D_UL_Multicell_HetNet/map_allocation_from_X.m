function [User_to_PRBs_Map, PRB_to_Users_Map] = map_allocation_from_X(X, problem)
% MAP_ALLOCATION_FROM_X: Converte a MATRIZ X (3 x N_RBs) em mapas.
%   Agora suporta a estrutura Underlay (3 linhas).
%
%   X: Matriz [3 x N_RBs].
%      Linha 1: CUEs
%      Linhas 2-3: D2Ds
%
% Saída:
%   User_to_PRBs_Map: Cell Array {1 x N_users} -> Lista de PRBs de cada um.
%   PRB_to_Users_Map: Cell Array {1 x N_RBs} -> Quem está no PRB r.

    nTotalUsers = problem.nCUE + problem.nD2D;
    nRBs = problem.nRBs;
    
    % 1. Inicializa Mapas
    User_to_PRBs_Map = cell(1, nTotalUsers);
    PRB_to_Users_Map = cell(1, nRBs);
    
    % 2. Varredura da Matriz X (3 linhas x nRBs colunas)
    for r = 1:nRBs
        
        % Quem está neste PRB? (Pega os valores das 3 linhas)
        users_in_prb = [];
        
        for row = 1:3
            u = X(row, r);
            
            if u > 0
                % Adiciona este PRB à lista do utilizador u
                User_to_PRBs_Map{u} = [User_to_PRBs_Map{u}, r];
                
                % Adiciona o utilizador u à lista deste PRB
                users_in_prb = [users_in_prb, u];
            end
        end
        
        % Guarda a lista completa de ocupantes deste PRB (CUE + D2Ds)
        PRB_to_Users_Map{r} = users_in_prb;
    end
end