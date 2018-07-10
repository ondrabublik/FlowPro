function generuj_sit4
np = 80; % pocet elementu profilu
n1 = 2*np;
n2 = 60;

X = zeros(n1,n2);
Y = zeros(n1,n2);

data = NACA(np,0);
% data = OVAL(2*np);

I = 1:n1;
X(I,1) = data(I,1);
Y(I,1) = data(I,2);

[X,Y] = Hyp_gen4(X,Y);
% [X,Y] = Laplace([X;X(1,:)],[Y;Y(1,:)]);
% X = X(:,1:n2);
% Y = Y(:,1:n2);

% tisk ctvercove site
figure;
hold on;
for i = 1:n1
    plot(X(i,:),Y(i,:));
end
for j = 1:n2
    plot([X(:,j);X(1,j)],[Y(:,j);Y(1,j)]);
end
plot(data(:,1),data(:,2),'.r');
axis equal;

% tvorba trojuhelnikove site
[n1,n2] = size(X);
np = n1*n2;
PX = zeros(np,1);
PY = zeros(np,1);
ind = zeros(n1+1,n2);
k = 1;
for i = 1:n1+1
    if(i == n1+1)
        for j = 1:n2
            ind(i,j) = ind(1,j);
        end
    else
        for j = 1:n2
            ind(i,j) = k;
            PX(k) = X(i,j);
            PY(k) = Y(i,j);
            k = k + 1;
        end
    end
end

nt = 2*(n1-1)*(n2-1);
TP = zeros(nt,3);
k = 1;
for i = 1:n1
    for j = 1:n2-1
        TP(k,:) = [ind(i,j),ind(i+1,j),ind(i,j+1)];
        TP(k+1,:) = [ind(i+1,j),ind(i+1,j+1),ind(i,j+1)];
        k = k + 2;
    end
end

figure
triplot(TP,PX,PY);
axis equal;

P = [PX,PY];
save 'P' P;
T = TP;
save 'T' T;

mesh{1} = P;
mesh{2} = T;
save mesh mesh;

%__________________________________________________________________________
function [Xp,Yp] = Hyp_gen4(X,Y)
n1 = length(X(:,1));
n2 = length(X(1,:));
K = 12;
lambda = 8;

de = 1/n1;
dn = 3*1/n2;

% explicitni
jexpl = 2;
for j = 1:jexpl
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
        
        g11 = (((X(ip,j)-X(im,j))/(2*de)).^2 + ((Y(ip,j)-Y(im,j))/(2*de)).^2);
        
        V = K*sqrt(g11)*exp(-lambda*(1-j/n2));

        X(i,j+1) = X(i,j) - dn*(V/g11)*(Y(ip,j) - Y(im,j))/(2*de);
        Y(i,j+1) = Y(i,j) + dn*(V/g11)*(X(ip,j) - X(im,j))/(2*de);
    end
end

% implicitni
A = sparse(2*n1,2*n1);
b = zeros(2*n1,1);
koef = dn/(2*de);
I1 = 1:2:2*n1;
I2 = 2:2:2*n1;
for j = jexpl+1:n2-1
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
