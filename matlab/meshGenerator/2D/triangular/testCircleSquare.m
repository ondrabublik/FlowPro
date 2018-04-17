function data = testCircleSquare

alfa = 5;
n = 30;
t1 = linspace(0,pi/2-alfa/180*pi,n);
t2 = linspace(pi/2+alfa/180*pi,pi,n);

data = [-1 1 cos(t1) cos(t1(end)) cos(t2(1)) cos(t2); ...
        -1 -1 sin(t1) 0 0 sin(t2)]';