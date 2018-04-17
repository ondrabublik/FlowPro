function ropeBDF
global n h k kt m b G
L = 1;
n = 10;
h = L/(n-1);
x = linspace(0,L,n)';
y = zeros(n,1);

u = zeros(4*n,1);
for i = 1:n
    u(4*(i-1) + 1) = x(i);
    u(4*(i-1) + 2) = y(i);
end

m = 0.1/n;
k = 10000;
kt = 1000;
b = 0;
G = [0, -9.81];

Ix = 1:4:4*n;
Iy = 2:4:4*n;

t = 0;
dt = 0.001;

figure
tisk = plot(u(Ix),u(Iy),'marker','o');
axis equal
axis([-1 1 -2 0.1]);

iter = 10000;
dR = sparse(4*n,4*n);
E = speye(4*n,4*n);
hd = 1e-6;
for op = 1:iter
    
%     K1 = fun(u);
%     K2 = fun(u + dt/2*K1);
%     u = u + dt*K2;

    r = fun(u);
    for i = 1:4*n
        u(i) = u(i) + hd;
        dR(:,i) = (fun(u) - r)/hd;
        u(i) = u(i) - hd;
    end

    u = u + (E/dt - dR)\r;
    
    set(tisk,'Xdata',u(Ix));
    set(tisk,'Ydata',u(Iy));
    drawnow;
%     pause(0.05);
   
    t = t + dt;
end

function R = fun(u)
global n h k kt m b G

   R = zeros(4*n,1);
   I = 1:2;
   for i = 2:n
       rm = u(4*(i-2)+I) - u(4*(i-1)+I);
       rn = norm(rm);
       rm = rm/rn*(rn-h);
       F = k*rm;
       if(i < n)
           r = (u(4*(i-2)+I) + u(4*i+I))/2 - u(4*(i-1)+I);
           rp = u(4*i+I) - u(4*(i-1)+I);
           rn = norm(rp);
           rp = rp/rn*(rn-h);
           F = F + kt*r + k*rp;
       end
       R(4*(i-1) + 1) = u(4*(i-1) + 3);
       R(4*(i-1) + 2) = u(4*(i-1) + 4);
       R(4*(i-1) + 3) = (F(1) - b*u(4*(i-1) + 3) + G(1))/m;
       R(4*(i-1) + 4) = (F(2) - b*u(4*(i-1) + 4) + G(2))/m;
   end
   R(1) = 0;
   R(2) = 0;
   
   
   