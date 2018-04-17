function CUBE

n1 = 10;
n2 = 10;
n3 = 10;

display(['Pocet bunek site je: ', num2str(n1*n2*n3)]);

X = zeros(n1,n2,n3);
Y = zeros(n1,n2,n3);
Z = zeros(n1,n2,n3);

z = linspace(0,1,n3);
I = 1:n1;
J = 1:n2;
for k = 1:n3
    [Xp,Yp] = fun(n1,n2,z(k));
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
bound = [];

save 'P' P;
save 'QP' QP;

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


dlmwrite('geometry/vertices.txt', P, 'delimiter', ' ', 'precision', 16);
dlmwrite('geometry/meshIndexes.txt', QP-1, 'delimiter', ' ');
dlmwrite('geometry/boundaryType.txt', bound(:,[5,3,4]), 'delimiter', ' ');
dlmwrite('geometry/elementType.txt', 4*ones(size(TP,1),1), 'delimiter', ' ');


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


function [X,Y] = fun(n1,n2,z)
X = zeros(n1,n2);
Y = zeros(n1,n2);

% prvni cast
x = linspace(-0.5,0.5,n1);
y = linspace(-0.5,0.5,n2);
for i = 1:n1
    for j = 1:n2
        X(i,j) = x(i);
        Y(i,j) = y(j);
    end
end














