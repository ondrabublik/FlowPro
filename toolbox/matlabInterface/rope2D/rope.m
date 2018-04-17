function rope
global n h k m b G
L = 1;
n = 10;
h = L/(n-1);
X = [linspace(0,L,n)',zeros(n,1)];
V = zeros(size(X));
G = zeros(size(X));
G(:,2) = -9.81;

m = 0.5/n;
k = 1000;
b = 0.1;

t = 0;
dt = 0.001;

figure
tisk = plot(X(:,1),X(:,2),'marker','o');
axis equal
axis([-1 1 -2 0.1]);

iter = 10000;
for op = 1:iter
   
    [K1x,K1v] = fun(X,V);
    [K2x,K2v] = fun(X+dt/2*K1x,V+dt/2*K1v);
    X = X + dt*K2x;
    V = V + dt*K2v;
   
    set(tisk,'Xdata',X(:,1));
    set(tisk,'Ydata',X(:,2));
    drawnow;
%     pause(0.05);
   
    t = t + dt;
end

function [Kx,Kv] = fun(X,V)
global n h k m b G

   F = zeros(n,2);
   for i = 2:n
       rm = X(i-1,:) - X(i,:);
       rn = norm(rm);
       rm = rm/rn*(rn-h);
       F(i,:) = k*rm;
       if(i < n)
           r = (X(i-1,:) + X(i+1,:))/2 - X(i,:);
           rp = X(i+1,:) - X(i,:);
           rn = norm(rp);
           rp = rp/rn*(rn-h);
           F(i,:) = F(i,:) + k*r + k*rp;
       end
   end

   Kx = V;
   Kv = (F - b*V + G)/m;
   Kx(1,:) = 0;
   Kv(1,:) = 0;
   
   
   
   
   