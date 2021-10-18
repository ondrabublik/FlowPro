function data = halfCylinder(n,r)

t1 = linspace(0,pi,n);
t2 = linspace(r,-r,round(n/2));

x1 = r*cos(t1);
y1 = r*sin(t1);
x2 = linspace(-r,r,round(n/2));
y2 = zeros(size(x2));

x = [x1(2:end-1),x2]';
y = [y1(2:end-1),y2]';

data = [x,y];
