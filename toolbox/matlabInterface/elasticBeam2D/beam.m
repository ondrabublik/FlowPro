function beam
L = 0.5;
n = 10;
h = L/(n-1);
E = 0.1;
Is = 1;
nu = 2*n;
M = zeros(nu,nu);
K = zeros(nu,nu);

rho = 100;
Ap = 0.1;
y = linspace(0,L,n);
         
u = zeros(nu-2,1);
uOld = zeros(nu-2,1);
uOld2 = zeros(nu-2,1);

a = 4818;
g = 729 * h;
c = 1482;
d = -321 * h;
e = 172 * h * h;
f = -73 * h * h;
Me = rho*Ap*h/12600*[a, g, c, d; g, e, -d, f; c, -d, a, -g; d, f, -g, e];
Ke = E*Is/h^3*[12 6*h -12 6*h; 6*h 4*h^2 -6*h 2*h^2; -12 -6*h 12 -6*h; 6*h 2*h^2 -6*h 4*h^2];
Fe = h/2*[1,h/6,1,-h/6]';
k = 1:4;
l = 1:4;
for i = 1:n-1
    M((i-1)*2 + k, (i-1)*2 + l) = M((i-1)*2 + k, (i-1)*2 + l) + Me(k,l);
    K((i-1)*2 + k, (i-1)*2 + l) = K((i-1)*2 + k, (i-1)*2 + l) + Ke(k,l);
end
I = 3:nu;
M = M(I,I);
K = K(I,I);
alfa = 0;
beta = 0;
B = alfa*M + beta*K;
dt = 0.01;
Ai = inv(M/dt^2 + B/dt + K);

t = 0;
figure
J = 1:2:nu-2;
tisk = plot(y,[0;u(J)]);
axis([0 0.5 -0.5 0.5]);
iter = fix(50/dt);
tt = zeros(iter,2);
for op = 1:iter
    F = zeros(nu,1);
    if(t < 1)
        q0 = 5;
    else
        q0 = 0;
    end
    for i = 1:n-1
        F((i-1)*2 + k) = F((i-1)*2 + k) + q0*Fe(k);
    end
    F = F(I);
    u = Ai*(M*(2*uOld - uOld2)/dt^2 + B*uOld/dt + F);
    uOld2 = uOld;
    uOld = u;
    t = t + dt;
    tt(op,:) = [t,u(end)];
    w = [0;u(J)];
    set(tisk,'Ydata',w);
    dh = zeros(n,1);
    for i = 2:n
        dh(i) = h - sqrt(max(h^2 - (w(i)-w(i-1))^2,0));
    end
    x = y-dh';
    set(tisk,'Xdata',x);
    
    sum = 0;
    for i = 2:n
        sum = sum + sqrt((y(i)-y(i-1))^2 + (x(i)-x(i-1))^2);
    end
    sum
    drawnow;
end

figure
plot(tt(:,1),tt(:,2));

% axis equal;
