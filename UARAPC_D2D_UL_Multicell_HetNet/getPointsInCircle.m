function pts = getPointsInCircle(radius, n)
% getPointsInCircle   Gera n pontos uniformes dentro de um disco de dado raio r 
%   pts = getPointsInCircle(radius, n) retorna um n×2 com coordenadas [x,y]
%   distribuídas uniformemente no disco centrado na origem.

  % raios distribuídos segundo a raiz quadrada de rand para uniformidade
  r = radius * sqrt(rand(n,1));
  % ângulos uniformes de 0 a 2π
  theta = 2*pi*rand(n,1);
  % coordenadas cartesianas
  x = r .* cos(theta);
  y = r .* sin(theta);
  pts = [x, y];
end
