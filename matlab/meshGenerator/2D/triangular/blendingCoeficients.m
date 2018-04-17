function blendingCoeficients(cesta,PX,PY,TP,TT)
TP = TP(:,1:3) + 1;
TT = TT(:,1:3);

TE_ALE = zeros(size(TT));
for i = 1:length(TT(:,1))
    xs = sum(PX(TP(i,1:3)))/3;
    ys = sum(PY(TP(i,1:3)))/3;
    r = sqrt(xs^2+ys^2);
    for j = 1:3
        if(TT(i,j) == -1 && r < 2)
            TE_ALE(i,j) = 2;
        end
    end

%     for j = 1:3
%         if(TT(i,j) == -1 && ys > 0.1)
%             TE_ALE(i,j) = 2;
%         end
%         if(TT(i,j) == -1 && ys > -0.1 && ys < 0.1)
%             TE_ALE(i,j) = 3;
%         end
%         if(TT(i,j) == -1 && ys < -0.1)
%             TE_ALE(i,j) = 4;
%         end
%     end
end
n_profiles = max(max(TE_ALE))-1;

% reseni poissonovy rovnice
np = length(PX);
nt = length(TP(:,1));
Ptyp = zeros(np,1);
poc = zeros(np,1);
A = sparse(np,np);
plus = [2 3 1];
for i = 1:nt
    for j = 1:3
        jp = plus(j);
        S = sqrt((PX(TP(i,jp)) - PX(TP(i,j)))^2 + (PY(TP(i,jp)) - PY(TP(i,j)))^2);
        A(TP(i,jp),TP(i,j)) = 1/S;
        A(TP(i,jp),TP(i,jp)) = A(TP(i,jp),TP(i,jp)) - 1/S;
        poc(TP(i,jp)) =  poc(TP(i,jp)) + 1/S;
        
        if(TT(i,j) < 0)
            Ptyp(TP(i,jp)) = 1;
            Ptyp(TP(i,j)) = 1;
            if(TE_ALE(i,j) > 0)
                Ptyp(TP(i,jp)) = TE_ALE(i,j);
                Ptyp(TP(i,j)) = TE_ALE(i,j);
            end
        end
    end
end

for i = 1:np
    A(i,:) = A(i,:)/poc(i);
end

I = 1:np;
Iin = I(Ptyp == 0);
Iste = I(Ptyp == 1);
z = zeros(np,n_profiles);
for k = 1:n_profiles
    z(Iste,k) = 0;
    Ih = I(Ptyp == k+1);
    z(Ih,k) = 1;
    
    b = -A(Iin,I)*z(:,k);
    z(Iin,k) = A(Iin,Iin)\b;
    
    z(z(:,k) < 0,k) = 0;
    figure('color','w')
    tricontf(PX,PY,TP,z(:,k),10);
    axis equal
    box on
    set(gca,'fontsize',14)
    colorbar
end

dlmwrite(strcat(cesta, 'TEale.txt'), TE_ALE, 'delimiter', ' ');
dlmwrite(strcat(cesta, 'blendingFunctions.txt'), z, 'delimiter', ' ');






