function generuj_sit5
np = 50; % pocet elementu profilu
ns = 20; % pocet elementu symetrie
n1 = 2*np + 2*ns;
n2 = 40;

X = zeros(n1+1,n2);
Y = zeros(n1+1,n2);

data = NACA(np);
xs1 = linspace(2,0,ns+1)';
xs2 = linspace(0,2,ns+1)';

% spodni stena
X(:,1) = [xs1;data(2:2*np,1)-1;xs2];
Y(:,1) = [zeros(ns+1,1);data(2:2*np,2);zeros(ns+1,1)];

[X,Y] = Hyp_gen2(X,Y);
[X,Y] = Laplace(X,Y);

figure;
hold on;
for i = 1:n1+1
    plot(X(i,:),Y(i,:));
end
for j = 1:n2
    plot([X(:,j);X(1,j)],[Y(:,j);Y(1,j)]);
end
plot(X(:,n2),Y(:,n2),'.r');
axis equal;
size(X)

% tvorba trojuhelnikove site
[n1,n2] = size(X);
np = n1*n2;
PX = zeros(np,1);
PY = zeros(np,1);
ind = zeros(n1,n2);
k = 1;
for i = 1:n1
    for j = 1:n2
        ind(i,j) = k;
        PX(k) = X(i,j);
        PY(k) = Y(i,j);
        k = k + 1;
    end
end

nt = 2*(n1-1)*(n2-1);
TP = zeros(nt,3);
k = 1;
for i = 1:n1-1
    for j = 1:n2-1
        TP(k,:) = [ind(i,j),ind(i+1,j),ind(i+1,j+1)];
        TP(k+1,:) = [ind(i+1,j+1),ind(i,j+1),ind(i,j)];
        k = k + 2;
    end
end

figure
triplot(TP,PX,PY);
axis equal;


%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen2(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 5;
lambda = 5;

de = 1/n1;
dn = 2/n2;

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

%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen4(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 10;
lambda = 4;

de = 1/n1;
dn = 2/n2;

% explicitni
a = 0.5;
for j = 1:2
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

% implicitni
A = sparse(2*n1,2*n1);
b = zeros(2*n1,1);
koef = dn/(2*de);
I1 = 1:2:2*n1;
I2 = 2:2:2*n1;
for j = 2:n2-1
    for i = 1:n1
        if(i == n1)
            ip = 1;
        else
            ip = i+1;
        end
        if(i == 1)
            im = n1;
        else
            im = i-1;
        end
        
        dxde = (X(ip,j)-X(im,j))/(2*de);
        dyde = (Y(ip,j)-Y(im,j))/(2*de);
        
        g11 = dxde^2 + dyde^2;
        
        V0 = K*sqrt(g11)*exp(-lambda*(1-j/n2));
        
        dxdn = -dyde*V0/g11;
        dydn = dxde*V0/g11;
        
        ap = dxde*dxdn - dyde*dydn;
        bp = dxde*dydn + dxdn*dyde;
        cp = g11;
        
        C = 1/cp*[ap bp; bp -ap];
        tlum = -sqrt((ap^2 + bp^2)/cp^2);
        
        V = K*sqrt(g11)*exp(-lambda*(1-(j+1)/n2));
        
        S = (V+V0)/cp*[-dyde; dxde];
        
        I = 2*(i-1)+1:2*i;
        Ip = 2*(ip-1)+1:2*ip;
        Im = 2*(im-1)+1:2*im;
        
        A(I,I) = (1-2*tlum)*eye(2);
        A(I,Ip) = koef*C + tlum*eye(2);
        A(I,Im) = -koef*C + tlum*eye(2);
        b(I) = dn*S + [X(i,j);Y(i,j)];
    end
    R = A\b;

    X(:,j+1) = R(I1);
    Y(:,j+1) = R(I2);
end

Xp = X;
Yp = Y;

%__________________________________________________________________________
function [Xp,Yp] = Laplace(X,Y)
[n1,n2] = size(X);
I = 2:n1-1;
J = 2:n2-1;
for op = 1:10
    X(I,J) = 1/4*(X(I+1,J)+X(I,J+1)+X(I-1,J)+X(I,J-1));
    Y(I,J) = 1/4*(Y(I+1,J)+Y(I,J+1)+Y(I-1,J)+Y(I,J-1));
end
Xp = X;
Yp = Y;
