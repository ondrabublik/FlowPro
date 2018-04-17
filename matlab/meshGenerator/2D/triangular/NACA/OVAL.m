function data = OVAL(n)
n = n+1;
t = linspace(0,1,n)';
t = t;
a = 0.5;
b = 0.05;
x = -a*cos(2*pi*t);
y = b*sin(2*pi*t);

figure;
plot(x,y,'.');
axis equal
data = [x(1:n-1),y(1:n-1)];