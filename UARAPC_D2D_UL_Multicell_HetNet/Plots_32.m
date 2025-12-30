%% 6. VISUALIZAÇÃO FINAL (Gráficos)
fprintf('\n--- Gerando Gráficos... ---\n');

% FIG 1 e FIG 3 ficam iguais...

% ---------------- FIG 4: Convergência MÉDIA (Fitness) + STD em 5,10,15 -----
figure(4); clf; hold on; grid on;
it_marks = [5 10 15];                          % iterações onde terá STD
it_marks = it_marks(it_marks <= params.maxIter);  % garante que não extrapola

for a = 1:nAlgos
    name = algorithms{a};

    mean_convFit = mean(convFit_all.(name), 1);   % média em cada iteração
    std_fit      = std_convFit.(name);           % std já calculado antes

    % plota curva média
    h = plot(1:params.maxIter, mean_convFit, 'LineWidth',1.5, ...
        'DisplayName', name);

    % mesma cor da curva para as barrinhas
    col = get(h, 'Color');

    % adiciona barras de erro apenas nas iterações 5,10,15
    errorbar(it_marks, mean_convFit(it_marks), std_fit(it_marks), ...
        'LineStyle','none', ...        % só as barras/markers, sem linha
        'Marker','o', 'MarkerSize',6, ...
        'Color', col, ...
        'HandleVisibility','off');     % não repetir na legenda
end
xlabel('Iteração');
ylabel('Fitness (Com Penalidades)');
title(sprintf('Fig 4: Convergência MÉDIA do Fitness (%d execuções)', nRuns));
legend('show','Location','best');

% ------------- FIG 5: Convergência MÉDIA (Throughput) + STD em 5,10,15 ----
figure(5); clf; hold on; grid on;
for a = 1:nAlgos
    name = algorithms{a};

    mean_convObj_Mbps = mean(convObj_all.(name), 1) / 1e6; % média em Mbps
    std_obj_Mbps      = std_convObj.(name) / 1e6;          % std em Mbps

    % curva média
    h = plot(1:params.maxIter, mean_convObj_Mbps, ...
        'LineWidth',1.5, 'DisplayName', name);

    col = get(h,'Color');

    % barras de erro só em 5,10,15
    errorbar(it_marks, mean_convObj_Mbps(it_marks), std_obj_Mbps(it_marks), ...
        'LineStyle','none', ...
        'Marker','s','MarkerSize',6, ...
        'Color',col, ...
        'HandleVisibility','off');
end
xlabel('Iteração');
ylabel('Throughput Total MÉDIO (Mbps)');
title(sprintf('Fig 5: Função Objetivo - Throughput REAL MÉDIO do Sistema (%d execuções)', nRuns));
legend('show','Location','best');

fprintf('Visualização Concluída.\n');
toc;
