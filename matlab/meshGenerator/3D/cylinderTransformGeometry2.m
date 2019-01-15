function cylinderTransformGeometry2

n1 = 10; % pocet bodu site profilu
n2 = 15;
n3 = 5;
L = 2;
r1 = 0.5;
r2 = 1.5;

display(['Pocet bunek site je: ', num2str((n1-1)*(n2-1)*(n3-1))]);

X = zeros(n1,n2,n3);
Y = zeros(n1,n2,n3);
Z = zeros(n1,n2,n3);

x = linspace(0,L,n1);
y = linspace(0,pi/2,n2);
z = linspace(r1,r2,n3);
for i = 1:n1
    for j = 1:n2
        for k = 1:n3
            X(i,j,k) = x(i);
            Y(i,j,k) = y(j);
            Z(i,j,k) = z(k);
        end
    end
end

ind = zeros(n1,n2,n3);
P = zeros(n1*n2*n3,3);
s = 1;
for j = 1:n2
    for i = 1:n1    
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
        if(abs(xs-L) < eps)
            bound = [bound; [-3, QP(i,index)-1]];
        end 
    end
end

mkdir('mesh')
dlmwrite('mesh/vertices.txt', P, 'delimiter', ' ', 'precision', 16);
dlmwrite('mesh/elements.txt', QP-1, 'delimiter', ' ');
dlmwrite('mesh/boundaryType.txt', bound, 'delimiter', ' ');
dlmwrite('mesh/elementType.txt', 63*ones(size(QP,1),1), 'delimiter', ' ');
dlmwrite('mesh/elementType6.txt', 6*ones(size(QP,1),1), 'delimiter', ' ');

if(1 == 1)
    figure;
    hold on;
    for k = 1:n3;
        for i = 1:n1
            plot3(X(i,:,k),Y(i,:,k),Z(i,:,k),'k');
        end
        for j = 1:n2
            plot3(X(:,j,k),Y(:,j,k),Z(:,j,k),'k');
        end
    end
    axis equal;
end


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











