function data = lavalNozzle(n)

x0 = -1.3;
t = linspace(-1,2,n);
s = (1-atan(10*t)/pi*2)/2;
xL = t;
yL = s*1 + (1-s).*sqrt(abs(t)+0.1)/1.5;

x = [xL, xL(end), x0, x0, xL(end), xL(end:-1:1)];
y = [yL, yL(end)+0.1, yL(end)+0.1, -(yL(end)+0.1), -(yL(end)+0.1), -yL(end:-1:1)];

data = [y',-x'];

plot(x,y)