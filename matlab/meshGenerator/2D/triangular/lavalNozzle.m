function data = lavalNozzle(n)

x0 = -2.3;
t = linspace(-1,2,n);
s = (1-atan(10*t)/pi*2)/2;
xL = 2*t;
ymax = 1.8;
yL = s*ymax + (1-s).*sqrt(abs(t)+0.1)/1.5;

a = 1;
x = [xL, xL(end), x0, x0, xL(end), xL(end:-1:1)];
y = [yL, yL(end)+a, yL(end)+a, -(yL(end)+a), -(yL(end)+a), -yL(end:-1:1)];

data = [x',y'];

plot(x,y)