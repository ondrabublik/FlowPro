function data = gamm(n)
% funkce generuje geometrii GAMM kanalu

x = linspace(1,2,n);
h = 0.1;
ys = (h-(0.5^2)/h)/2;
xs = 1.5;
R = (ys^2 + 0.5^2)^(1/2);
y = sqrt(R^2-(x-xs).^2) + ys;

x = [x,3,3,0,0];
y = [y,0,1,1,0];

data = [x',y'];
