function ebeam
mu = 1;
EI = 1;
L = 0.5;

n = 10;
u = zeros(n,2);
uold = u;

x = linspace(0,L,n)';
y = zeros(size(x));

dt = 0.001;
I = 2:n;

h2 = (L/(n-1))^2;
D = sparse(n,n);
D(2,1) = -1/h2;
D(2,2) = 1/h2;
for i = 3:n-1
    D(i,i-1) = 1/h2;
    D(i,i) = -2/h2;
    D(i,i+1) = 1/h2;
end
D2 = D*D;
D2 = D2(I,I);
A = [-1/mu*EI*D2, zeros(size(D2)); zeros(size(D2)), ones(size(D2))];
spy(A)
return

figure
tisk = plot(x,y+u(:,2));
t = 0;
for i = 1:10
    q = ones(n,1)*sin(t);
    for s = 1:10
        u(I,1) = uold(I,1) - dt/mu*(EI*D2*u(I,2) + q(I));
        u(I,2) = uold(I,2) + dt*u(I,1);
    end
    t = t + dt;
    
    set(tisk,'ydata',y+u(:,2));
    drawnow
    pause(0.1);
end







