function gain = computeGain(distKm, linkType, physParams)
% computeGain: calcula o ganho de canal incluindo path-loss, shadowing, fading e ganhos de antena
%  distKm    : vetor ou matriz de distâncias em km
%  linkType  : 'cellular' ou 'D2D'
%  physParams: struct com campos:
%      .plCell.a, .plCell.b       : coeficientes path-loss célula (a=128.1, b=37.6)
%      .plD2D.a, .plD2D.b         : coeficientes path-loss D2D (a=148, b=40)
%      .shadowStdCell             : σ de shadowing célula em dB
%      .shadowStdD2D              : σ de shadowing D2D em dB
%      .gainTx_eNB, .gainRx_UE    : ganhos de antena em dBi

% 1. Path-loss e shadowing (log-normal)
if strcmpi(linkType,'cellular')
    a = physParams.plCell.a;
    b = physParams.plCell.b;
    sigma = physParams.shadowStdCell;
else
    a = physParams.plD2D.a;
    b = physParams.plD2D.b;
    sigma = physParams.shadowStdD2D;
end
% PL em dB
PL_dB = a + b * log10(distKm);
% shadowing em dB
shadow_dB = sigma .* randn(size(distKm));
% ganho combinado em dB (negativo)
combinedPL_dB = -(PL_dB + shadow_dB);
% converte para escala linear
lossLinear = 10.^(combinedPL_dB/10);

% 2. Small-scale fading (Rayleigh)
% h ~ CN(0,1) -> |h|^2
h = (randn(size(distKm)) + 1j*randn(size(distKm))) / sqrt(2);
fading = abs(h).^2;

% 3. Ganhos de antena
% Para cellular: Tx=CUE(0dBi), Rx=eNB(14dBi); para D2D ambos 0dBi (UE)
gUE = 10^(physParams.gainRx_UE/10);      % UE gain (0 dBi)
gENB = 10^(physParams.gainTx_eNB/10);    % eNB gain (14 dBi)
% usar produto de ambos para qualquer link
gAnt = gUE * gENB;

% 4. Ganho total do canal
gain = lossLinear .* fading * gAnt;
end

%% Função de taxa de Shannon (opcional)
function R = shannonRate(sinr, bw)
% sinr: relação sinal-ruído linear
% bw  : largura de banda [Hz]
R = bw * log2(1 + sinr);
end

%% Função de SINR FullCSI (opcional)
function sinr = computeSINR_FullCSI(assign, problem, params)
% assign: vetor [1 x nD2D] com índice do CUE reutilizado (0=nenhum)
% problem: struct retornado por RRAProblem_FullCSI
% params : struct de parâmetros (contém txPowerCUE, txPowerD2D)

% Notação:
% cada CUE gera interferências se D2Ds reusearem seu RB
% loop para cada D2D i: sinal = pD2D * gain_D2D_pair(i)
% interferência de CUE de reuse j: pCUE * gain_D2D_to_CUE(i,j)

nD2D = problem.nD2D;
sinr = zeros(nD2D,1);
for i = 1:nD2D
    % sinal útil
    signal = db2pow(params.txPowerD2D) * problem.gain_D2D_pair(i);
    % ruído térmico
    noise = db2pow(params.physParams.noiseDensity) * problem.bwRB;
    % interferência CUE
    if assign(i)>0
        j = assign(i);
        interfCUE = db2pow(params.txPowerCUE) * problem.gain_D2D_to_CUE(i,j);
    else
        interfCUE = 0;
    end
    % interferência de outros D2Ds (full CSI assume todos ativos)
    interfD2D = 0;
    for k = 1:nD2D
        if k~=i && assign(k)==assign(i)
            interfD2D = interfD2D + db2pow(params.txPowerD2D) * problem.gain_D2D_pair(k);
        end
    end
    sinr(i) = signal / (noise + interfCUE + interfD2D);
end
end
