function Cmesh

ns = 15; % pocet elementu symetrie
n2 = 30;

% data = NACA(np);
data = load('NACA64A010.txt');
I = size(data,1):-1:2;
data = data(I,:);
np = size(data,1)/2;
n1 = 2*np + 2*ns;
X = zeros(n1+1,n2);
Y = zeros(n1+1,n2);
xs1 = linspace(1,0,ns+1)';
xs2 = linspace(0,1,ns+1)';

% spodni stena
X(:,1) = [xs1;data(2:2*np,1)-1;xs2];
Y(:,1) = [zeros(ns+1,1);data(2:2*np,2);zeros(ns+1,1)];

[X,Y] = Hyp_gen2(X,Y);

figure;
hold on;
for i = 1:n1
    plot(X(i,:),Y(i,:));
end
for j = 1:n2
    plot([X(:,j);X(1,j)],[Y(:,j);Y(1,j)]);
end
axis equal;

%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen2(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 30;
lambda = 6;

de = 1/n1;
dn = 1/n2;

a = 0.5;
for j = 1:n2-1
    for i = 2:n1-1
        dx = (X(i+1,j)-X(i-1,j))/(2*de);
        dy = (Y(i+1,j)-Y(i-1,j))/(2*de);
        
        g11 = dx.^2 + dy.^2;
        
        V = K*sqrt(g11)*exp(-lambda*(1-j/n2));

        X(i,j+1) =  a*X(i,j) + (1-a)/2*(X(i+1,j) + X(i-1,j)) - dn*(V/g11)*dy;
        Y(i,j+1) =  a*Y(i,j) + (1-a)/2*(Y(i+1,j) + Y(i-1,j)) + dn*(V/g11)*dx;
    end
    g11 = (((X(2,j)-X(1,j))/de).^2 + ((Y(2,j)-Y(1,j))/de).^2);
    V = K*sqrt(g11)*exp(-lambda*(1-j/n2));
    X(1,j+1) = X(1,j) - dn*(V/g11)*(Y(2,j) - Y(1,j))/de;
    Y(1,j+1) = Y(1,j) + dn*(V/g11)*(X(2,j) - X(1,j))/de;
        
    g11 = (((X(n1,j)-X(n1-1,j))/de).^2 + ((Y(n1,j)-Y(n1-1,j))/de).^2);
    V = K*sqrt(g11)*exp(-lambda*(1-j/n2));
    X(n1,j+1) = X(n1,j) - dn*(V/g11)*(Y(n1,j) - Y(n1-1,j))/de;
    Y(n1,j+1) = Y(n1,j) + dn*(V/g11)*(X(n1,j) - X(n1-1,j))/de;
end

Xp = X;
Yp = Y;

function [Xp,Yp] = AlgebGen(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
F = zeros(n1,1);
N = zeros(n2,1);

for i = 1:n1
    F(i) = (i-1)/(n1-1);
end

for j = 1:n2
    N(j) = (j-1)/(n2-1);
end

Xp = zeros(n1,n2);
Yp = zeros(n1,n2);

for i = 1:n1
    for j = 1:n2
        Xp(i,j) = (1-F(i))*X(1,j)+F(i)*X(n1,j)+(1-N(j))*X(i,1)+N(j)*X(i,n2)-(1-F(i))*(1-N(j))*X(1,1)-F(i)*(1-N(j))*X(n1,1)-(1-F(i))*N(j)*X(1,n2)-F(i)*N(j)*X(n1,n2);
        Yp(i,j) = (1-F(i))*Y(1,j)+F(i)*Y(n1,j)+(1-N(j))*Y(i,1)+N(j)*Y(i,n2)-(1-F(i))*(1-N(j))*Y(1,1)-F(i)*(1-N(j))*Y(n1,1)-(1-F(i))*N(j)*Y(1,n2)-F(i)*N(j)*Y(n1,n2);
    end
end

function [Xp,Yp] = ElipticGen(X,Y)
Xp = X;
Yp = Y;
n1 = length(X(:,1));
n2 = length(X(1,:));
g11 = zeros(n1,n2);
g22 = zeros(n1,n2);
g12 = zeros(n1,n2);
xtemp = zeros(n1,n2);
ytemp = zeros(n1,n2);

[X,Y] = AlgebGen(X,Y);
I = 2:n1-1;
J = 2:n2-1;
n = 300;
for op = 1:n
    g11(I,J) = ((X(I+1,J)-X(I-1,J)).^2 + (Y(I+1,J)-Y(I-1,J)).^2)/4;
    g22(I,J) = ((X(I,J+1)-X(I,J-1)).^2 + (Y(I,J+1)-Y(I,J-1)).^2)/4;
    g12(I,J) = ((X(I+1,J)-X(I-1,J)).*(X(I,J+1)-X(I,J-1)) + (Y(I+1,J)-Y(I-1,J)).*(Y(I,J+1)-Y(I,J-1)))/4;
    
    xtemp(I,J) = 1./(2*(g11(I,J)+g22(I,J))).*(g22(I,J).*X(I+1,J)-0.5*g12(I,J).*X(I+1,J+1)+0.5*g12(I,J).*X(I+1,J-1)+g11(I,J).*X(I,J+1)+g11(I,J).*X(I,J-1)+g22(I,J).*X(I-1,J)-0.5*g12(I,J).*X(I-1,J-1)+0.5*g12(I,J).*X(I-1,J+1));
    ytemp(I,J) = 1./(2*(g11(I,J)+g22(I,J))).*(g22(I,J).*Y(I+1,J)-0.5*g12(I,J).*Y(I+1,J+1)+0.5*g12(I,J).*Y(I+1,J-1)+g11(I,J).*Y(I,J+1)+g11(I,J).*Y(I,J-1)+g22(I,J).*Y(I-1,J)-0.5*g12(I,J).*Y(I-1,J-1)+0.5*g12(I,J).*Y(I-1,J+1));
    
    err = sum(sum((X(I,J)-xtemp(I,J)).^2 + (Y(I,J)-ytemp(I,J)).^2));
    
    X(I,J) = xtemp(I,J);
    Y(I,J) = ytemp(I,J);
end
Xp = X;
Yp = Y;
