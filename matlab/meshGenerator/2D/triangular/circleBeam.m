function data = circleBeam(n,beta)
beta = beta/180*pi;
h = 0.01;
L = 1;
r = 0.1;
alfa = asin(h/(2*r));
t = linspace(alfa,2*pi-alfa,n)';

data = [r*sin(t(end)),L+r];
data = [data;[r*sin(t(1)),L+r]];
data = [data;[r*sin(t),r*cos(t)]];
data(:,2) = data(:,2)-r;

data = [cos(beta)*data(:,1) - sin(beta)*data(:,2), sin(beta)*data(:,1) + cos(beta)*data(:,2)];
