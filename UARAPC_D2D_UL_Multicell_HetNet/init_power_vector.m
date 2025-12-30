function P0 = init_power_vector(problem, params)
% INIT_POWER_VECTOR: inicializa vetor de potências (em Watts)
%   - CUEs: params.txPowerCUE (dBm)
%   - D2Ds: params.txPowerD2D (dBm)
%
% Saída:
%   P0: [nUsers x 1] em Watts

    nCUE  = problem.nCUE;
    nD2D  = problem.nD2D;
    nUser = nCUE + nD2D;

    % dBm -> Watts
    Pmax_CUE_W  = 10^((params.txPowerCUE - 30)/10);
    Pmax_D2D_W  = 10^((params.txPowerD2D - 30)/10);

    P0 = zeros(nUser,1);

    % CUEs
    P0(1:nCUE) = Pmax_CUE_W;

    % D2Ds
    if nD2D > 0
        P0(nCUE+1 : nUser) = Pmax_D2D_W;
    end
end
