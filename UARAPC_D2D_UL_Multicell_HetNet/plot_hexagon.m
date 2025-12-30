function h = plot_hexagon(center_x, center_y, R)
% plot_hexagon: Desenha um hexágono "deitado" (corrigido)
%   e retorna o handle do gráfico.

    % Ângulos rotacionados em 30 graus para "deitar" o hexágono
    angles = 30:60:390; 
    
    % Calcula as coordenadas (x, y) dos vértices
    vx = center_x + R * cosd(angles);
    vy = center_y + R * sind(angles);
    
    % Plota o hexágono (preto, pontilhado, fino)
    h = plot(vx, vy, 'k:', 'LineWidth', 0.5);
end