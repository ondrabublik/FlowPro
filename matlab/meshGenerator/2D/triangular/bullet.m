function data = bullet(n)

r = 1;
a = 3;
h = 5;
t = linspace(0,pi,n);
xtop = r*cos(t);
ytop = a*sin(t);

x = [r,xtop,-r];
y = [-h, ytop,-h];

data = [x',y'];

plot(x,y)
axis equal