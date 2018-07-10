function data = lavalNozzle2(n)

bx = [0 0 1 2 3];
by = [1 0 0 0.8 1];
t = linspace(0,1,n);
xL = zeros(length(t),1);
yL = zeros(length(t),1);
for i = 1:length(t)
    xL(i) = bezier(bx,t(i));
    yL(i) = bezier(by,t(i));
end

a = 0.1;
x0 = -1.7;
t = linspace(pi,0,round(n/2));
xtop = x0 - (yL(1)+a)*sin(t);
ytop = - (yL(1)+a)*cos(t);


xi = -1.5;
x = [xi, xL', xL(end), xtop, xL(end), xL(end:-1:1)',xi];
y = [yL(1), yL', yL(end)+a, ytop, -(yL(end)+a), -yL(end:-1:1)',-yL(1)];

data = [y',-x'];

% plot(x,y)
% axis equal


function z = bezier(x,t)
n = length(x);
y = x;
for k = 1:n-1
    for i = 2:n+1-k
        y(i-1) = (1-t)*x(i-1) + t*x(i);
    end
    x = y;
end
z = y(1);