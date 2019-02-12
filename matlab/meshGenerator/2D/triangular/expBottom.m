function data = expBottom
x = linspace(0,3,30)';
y = 0.15*exp(-(x-1.5).^2/0.1);
x = [x;3;0];
y = [y;1;1];
data = [x,y];