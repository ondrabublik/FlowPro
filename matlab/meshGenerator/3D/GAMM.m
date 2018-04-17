function GAMM
n11q = 7;
n12q = 10;
n13q = 7;
n2q = 15;
n3q = 15;

n1 = n11q+n12q+n13q + 1; % pocet bodu site profilu
n2 = n2q + 1;
n3 = n3q + 1;

display(['Pocet bunek site je: ', num2str(n1*n2*n3)]);

X = zeros(n1,n2,n3);
Y = zeros(n1,n2,n3);
Z = zeros(n1,n2,n3);

z = linspace(0,1,n3);
I = 1:n1;
J = 1:n2;
for k = 1:n3
    [Xp,Yp] = GAMM2D(n11q,n12q,n13q,n2q,1.5+z(k)/2);
    X(I,J,k) = Xp;
    Y(I,J,k) = Yp;
    Z(I,J,k) = z(k);
end

ind = zeros(n1,n2,n3);
P = zeros(n1*n2*n3,3);
s = 1;
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            ind(i,j,k) = s;
            P(s,:) = [X(i,j,k),Y(i,j,k),Z(i,j,k)];
            s = s + 1;
        end
    end
end

nq = (n1-1)*(n2-1)*(n3-1);
QP = zeros(nq,8);
s = 1;
for i = 1:n1-1
    for j = 1:n2-1
        for k = 1:n3-1
            v = [ind(i,j,k),ind(i+1,j,k),ind(i+1,j+1,k),ind(i,j+1,k),ind(i,j,k+1),ind(i+1,j,k+1),ind(i+1,j+1,k+1),ind(i,j+1,k+1)];
            QP(s,:) = v;
            s = s + 1;
        end
    end
end

% boundary condition
eps = 1e-5;
bound = [];
for i = 1:length(QP(:,1))
    for k = 1:6
        index = face(k);
        xs = sum(P(QP(i,index),1))/4;
        if(abs(xs) < eps)
            bound = [bound; [-2, QP(i,index)-1]];
        end
        if(abs(xs-3.5) < eps)
            bound = [bound; [-3, QP(i,index)-1]];
        end 
    end
end

dlmwrite('geometry/vertices.txt', P, 'delimiter', ' ', 'precision', 16);
dlmwrite('geometry/elements.txt', QP-1, 'delimiter', ' ');
dlmwrite('geometry/boundaryType.txt', bound, 'delimiter', ' ');
dlmwrite('geometry/elementType.txt', 6*ones(size(QP,1),1), 'delimiter', ' ');

if(1 == 1)
    figure;
    hold on;
    for k = 1:n3;
        for i = 1:n1
            plot3(X(i,:,k),Y(i,:,k),Z(i,:,k));
        end
        for j = 1:n2
            plot3(X(:,j,k),Y(:,j,k),Z(:,j,k));
        end
    end
    axis equal;
end


function [X,Y] = GAMM2D(n11,n12,n13,n2,xs)
n1 = n11+n12+n13;
X = zeros(n1+1,n2+1);
Y = zeros(n1+1,n2+1);

% prvni cast
x = linspace(0,xs-0.5,n11+1);
y = linspace(0,1,n2+1);
for i = 1:n11
    for j = 1:n2
        X(i,j) = x(i);
        Y(i,j) = y(i);
    end
end

% prvni cast
x = linspace(0,xs-0.5,n11+1);
y = linspace(0,1,n2+1);
for i = 1:n11+1
    for j = 1:n2+1
        X(i,j) = x(i);
        Y(i,j) = y(j);
    end
end

% druha cast
x = linspace(xs-0.5,xs+0.5,n12+1);
for i = 1:n12+1
    y = linspace(bulb(x(i)-xs),1,n2+1);
    for j = 1:n2+1
        X(i+n11,j) = x(i);
        Y(i+n11,j) = y(j);
    end
end

% treti cast
x = linspace(xs+0.5,3.5,n13+1);
y = linspace(0,1,n2+1);
for i = 1:n13+1
    for j = 1:n2+1
        X(i+n11+n12,j) = x(i);
        Y(i+n11+n12,j) = y(j);
    end
end

% figure
% hold on
% for i = 1:n1
%     for j = 1:n2
%         plot([X(i,j),X(i+1,j),X(i+1,j+1),X(i,j+1),X(i,j)],[Y(i,j),Y(i+1,j),Y(i+1,j+1),Y(i,j+1),Y(i,j)]);
%     end
% end


function y = bulb(x)
R = 26/20;
alfa = asin(x/R);
y = R*cos(alfa) - (R-0.1);


function ind = face(k)
    switch k
        case(1)
            ind = [0, 1, 2, 3]+1;
        case(2)
            ind = [4, 7, 6, 5]+1;
        case(3)
            ind = [0, 4, 5, 1]+1;
        case(4)
            ind = [1, 5, 6, 2]+1;
        case(5)
            ind = [2, 6, 7, 3]+1;
        case(6)
            ind = [3, 7, 4, 0]+1;
    end











