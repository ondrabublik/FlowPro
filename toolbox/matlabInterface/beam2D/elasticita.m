function elasticita
global nq QP QO np PX PY Ptyp;
clc;

% nacteni geometrie
load 'geometrie';
GP = geometrie{1};
GQ = geometrie{2};
GE = geometrie{3};

% deklarace potrebnych matic
deklaraceMatic(GP,GQ,GE);

% ridici hodnoty
E = 2.1e5;
mu = 0;
C1 = E/(1-mu^2);
C2 = mu;
lam = 0.5*(1-C2);
t = 0.01;
P = [1e-5,0]; % zatizeni

% sestaveni matice tuhosti
display('Sestavovani matice tuhosti...');
K = sparse(2*np,2*np);
F = zeros(2*np,1);
for i = 1:nq % opakovani pres jednotlive trojuhelniky
    x1 = PX(QP(i,1));
    y1 = PY(QP(i,1));
    x2 = PX(QP(i,2));
    y2 = PY(QP(i,2));
    x3 = PX(QP(i,3));
    y3 = PY(QP(i,3));
    Ke = matice_tuhosti(x1,y1,x2,y2,x3,y3,QO(i),C1,C2,lam,t);
    I = [2*QP(i,1)-1,2*QP(i,1), 2*QP(i,2)-1,2*QP(i,2), 2*QP(i,3)-1,2*QP(i,3)];
    K(I,I) = K(I,I) + Ke;
    
    Fe = zeros(6,1);
    for k = 1:3
        if(Ptyp(QP(i,k)) == 1)
            Fe = Fe + sila(PX(QP(i,k)),PY(QP(i,k)),x1,y1,x2,y2,x3,y3,P);
        end
    end
    
    F(I) = F(I) + Fe;
end

I = [];
for i = 1:np
    if(Ptyp(i) ~= -1)
        I = [I;2*i-1;2*i];
    end
end

x = zeros(2*np,1);

display('Vypocet...');
tic
% M = cholinc(K(I,I),'0');
% x(I) = cgs(K(I,I),F(I),1e-14,1000,M);
x(I) = K(I,I)\F(I);
toc

u = zeros(np,1);
v = zeros(np,1);
for i = 1:np
    u(i) = x(2*i-1);
    v(i) = x(2*i);
end 

S = napeti(u,v,C1,C2);

display('Tisk vysledku...');

figure;
title('Posunuti v x');
trisurf(QP,PX,PY,u);
shading interp;
view(2);
colorbar
axis equal;

figure;
title('Posunuti v y');
trisurf(QP,PX,PY,v);
shading interp;
view(2);
colorbar
axis equal;

figure;
hold on;
triplot(QP,PX,PY,'g');
triplot(QP,PX+u,PY+v,'k');
axis equal;

figure;
title('Redukovane napeti');
trisurf(QP,PX+u,PY+v,S);
shading interp;
view(2);
colorbar
axis equal;



%__________________________________________________________________________
function Ke = matice_tuhosti(x1,y1,x2,y2,x3,y3,O,C1,C2,lam,t)
x12 = x1 - x2;
y12 = y1 - y2;
x13 = x1 - x3;
y13 = y1 - y3;
x23 = x2 - x3;
y23 = y2 - y3;

Ke = C1*t/(4*O)* ...
    [        y23^2 + lam*x23^2,      -C2*x23*y23 - lam*y23*x23,        -y13*y23 - lam*x23*x13,     C2*x13*y23 + lam*y13*x23,         y23*y12 + lam*x23*x12,    -C2*y23*x12 - lam*x23*y12; ...
     -C2*x23*y23 - lam*y23*x23,              x23^2 + lam*y23^2,      C2*x23*y13 + lam*x13*y23,       -x23*x13 - lam*y23*y13,     -C2*x23*y12 - lam*y23*x12,        x12*x23 + lam*y23*y12; ...
        -y13*y23 - lam*x23*x13,       C2*x23*y13 + lam*x13*y23,             y13^2 + lam*x13^2,     C2*x13*y13 - lam*x13*y13,        -y13*y12 - lam*x13*x12,     C2*x12*y13 + lam*x13*y12; ...
      C2*x13*y23 + lam*y13*x23,         -x23*x13 - lam*y23*y13,      C2*x13*y13 - lam*x13*y13,            x13^2 + lam*y13^2,      C2*x13*y12 + lam*y13*x12,       -x13*x12 - lam*y13*y12; ...
         y23*y12 + lam*x23*x12,      -C2*x23*y12 - lam*y23*x12,        -y13*y12 - lam*x13*x12,     C2*x13*y12 + lam*y13*x12,             y12^2 + lam*x12^2,    -C2*x12*y12 - lam*x12*y12; ...
     -C2*y23*x12 - lam*x23*y12,          x12*x23 + lam*y23*y12,      C2*x12*y13 + lam*x13*y12,       -x13*x12 - lam*y13*y12,     -C2*x12*y12 - lam*x12*y12,            x12^2 + lam*y12^2];
 
 
 
%__________________________________________________________________________
function Fe = sila(x,y,x1,y1,x2,y2,x3,y3,P)
U = [x,y,1,0,0,0; 0,0,0,x,y,1];
S = [x1, y1,  1,  0,  0,  0; ...
      0,  0,  0, x1, y1,  1; ...
     x2, y2,  1,  0,  0,  0; ...
      0,  0,  0, x2, y2,  1; ...
     x3, y3,  1,  0,  0,  0; ...
      0,  0,  0, x3, y3,  1];
  
V = U*inv(S);
Fe = (P*V)';
 


%__________________________________________________________________________
function Sred = napeti(u,v,C1,C2)
global nq np QP PX PY

Se = zeros(3*nq,1);

B = [1 0 0 0 0 0; 0 0 0 0 1 0; 0 1 0 1 0 0];
D = C1*[1, C2, 0; C2, 1, 0; 0, 0, 0.5*(1-C2)];

for i = 1:nq
    x1 = PX(QP(i,1));
    y1 = PY(QP(i,1));
    x2 = PX(QP(i,2));
    y2 = PY(QP(i,2));
    x3 = PX(QP(i,3));
    y3 = PY(QP(i,3));
    S = [x1, y1,  1,  0,  0,  0; ...
          0,  0,  0, x1, y1,  1; ...
         x2, y2,  1,  0,  0,  0; ...
          0,  0,  0, x2, y2,  1; ...
         x3, y3,  1,  0,  0,  0; ...
          0,  0,  0, x3, y3,  1];
     
    U = [u(QP(i,1)),v(QP(i,1)),u(QP(i,2)),v(QP(i,2)),u(QP(i,3)),v(QP(i,3))]';  
    I = [3*i-2,3*i-1,3*i];
    Se(I) = D*B*inv(S)*U;
end

Sred_e = zeros(nq,1); % redukovane napeti
for i = 1:nq
    s1 = (Se(3*i-2) + Se(3*i-1))/2 + sqrt(((Se(3*i-2) - Se(3*i-1))/2)^2 + Se(3*1)^2); % hlavni napeti
    s2 = (Se(3*i-2) + Se(3*i-1))/2 - sqrt(((Se(3*i-2) - Se(3*i-1))/2)^2 + Se(3*1)^2);
    Sred_e(i) = sqrt(s1^2 + s2^2 - s1*s2);
end

% prevod do vrcholu
Sred = zeros(np,1);
poc = zeros(np,1);
for i = 1:nq
    for k = 1:3
        Sred(QP(i,k)) = Sred(QP(i,k)) + Sred_e(i);
        poc(QP(i,k)) = poc(QP(i,k)) + 1;
    end
end
Sred = Sred./poc;