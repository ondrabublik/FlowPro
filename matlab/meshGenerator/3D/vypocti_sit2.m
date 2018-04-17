function vypocti_sit2

load 'P'
load 'QP'
load 'QEP'

nq = length(QP);
Qs = zeros(nq,6,3);
for i = 1:nq
    for j = 1:6
        for s = 1:3
            Qs(i,j,s) = (P(QEP(i,j,1),s) + P(QEP(i,j,2),s) + P(QEP(i,j,3),s) + P(QEP(i,j,4),s))/4;
        end
    end
end

% hledani sousedu
QT = zeros(nq,3);
for i = 1:nq
    for j = 1:3
        QT(i,j) = (P(QP(i,1),j) + P(QP(i,2),j) +P(QP(i,3),j) +P(QP(i,4),j) +P(QP(i,5),j) +P(QP(i,6),j) +P(QP(i,7),j) +P(QP(i,8),j))/8;
    end
end

Qpom = zeros(nq*6,5);
s = 1;
for i = 1:nq
    for k = 1:6
        for j = 1:3
            Qpom(s,j) = Qs(i,k,j);
        end
        Qpom(s,4) = i;
        Qpom(s,5) = k;
        s = s + 1;
    end
end
Qpom = myqs(Qpom);

QQ = zeros(nq,6);
for i = 1:nq*6-1
    if(abs(Qpom(i,1)-Qpom(i+1,1)) < 1e-5)
        if(sqrt((Qpom(i,1)-Qpom(i+1,1))^2 + (Qpom(i,2)-Qpom(i+1,2))^2 + (Qpom(i,3)-Qpom(i+1,3))^2) < 1e-5)
            QQ(Qpom(i,4),Qpom(i,5)) = Qpom(i+1,4);
            QQ(Qpom(i+1,4),Qpom(i+1,5)) = Qpom(i,4);
        end
    end
end

QQ = QQ-1;
x_min = min(min(Qs(:,:,1)));
x_max = max(max(Qs(:,:,1)));
for i = 1:nq
    for k = 1:6
        if(QQ(i,k) < 0)
            QQ(i,k) = -1;
%             if(abs(Qs(i,k,1)-x_min) < 1e-3)
%                 QQ(i,k) = -2;
%             end
%             if(abs(Qs(i,k,1)-x_max) < 1e-3)
%                 QQ(i,k) = -3;
%             end
        end
    end
end

% reorganizace QEP
QEP_pom = zeros(nq,6*4);
for i = 1:nq
    for j = 1:6
        for k = 1:4
            QEP_pom(i,4*(j-1) + k) = QEP(i,j,k);
        end
    end
end
QEP = QEP_pom;

%___plneni geometrie_______________________________________________________
geometrie{1} = P(:,1);
geometrie{2} = P(:,2);
geometrie{3} = P(:,3);
geometrie{4} = QP-1;
geometrie{5} = QQ;
geometrie{6} = QEP-1;

save 'geometrie' geometrie;

dlmwrite('geometry/PXY.txt', P, 'delimiter', ' ', 'precision', 16);
save2File('geometry/TP.txt',QP-1)
save2File('geometry/TT.txt',QQ)
save2File('geometry/typ.txt',6*ones(length(QQ(:,1)),1))

figure
hold on
for i = 1:nq
    for k = 1:6
        if(QQ(i,k) < 0)
            switch QQ(i,k)
                case -1
                    plot3(Qs(i,k,1),Qs(i,k,2),Qs(i,k,3),'.k');
                case -2
                    plot3(Qs(i,k,1),Qs(i,k,2),Qs(i,k,3),'.g');
                case -3
                    plot3(Qs(i,k,1),Qs(i,k,2),Qs(i,k,3),'.r');
                case -4
                    plot3(Qs(i,k,1),Qs(i,k,2),Qs(i,k,3),'.m');
            end
        end
    end
end
axis equal;


function save2File(name,A)
fid = fopen(name,'w');
[m,n] = size(A);
for i = 1:m
    for j = 1:n
        fprintf(fid,'%i ',A(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);
