function PXY = deform(alpha)

[geometry, simulation] = loadArgs;
geometryPath = strcat('../../simulations/', geometry, '/');
meshPath = strcat(geometryPath,'/mesh/');
PXY0 = load('initMesh/vertices.txt');
TP = load([meshPath,'elements.txt']);
TP = TP(:,1:3)+1;
TT = load([meshPath,'neighbors.txt']);
TT = TT(:,1:3);

np = length(PXY0(:,1));
nt = length(TT(:,1));
prof = zeros(np,1);
Ibound = zeros(np,1);
for i = 1:nt
    for j = 1:3
        if(TT(i,j) == -1)
            jp = mod(j,3)+1;
            prof(TP(i,j)) = 1;
            prof(TP(i,jp)) = 1;
        end
        if(TT(i,j) < 0)
            jp = mod(j,3)+1;
            Ibound(TP(i,j)) = 1;
            Ibound(TP(i,jp)) = 1;
        end
    end
end
Iin = 1:np;
Iin = Iin(Ibound == 0);
Ib = 1:np;
Ib = Ib(Ibound == 1);
Ibody = 1:np;
Ibody = Ibody(prof == 1);

w = zeros(np,1);
for i = 1:nt
    for j = 1:3
        jp = mod(j,3)+1;
        L = 1;%sqrt((PXY0(TP(i,j),1)-PXY0(TP(i,jp),1))^2 + (PXY0(TP(i,j),2)-PXY0(TP(i,jp),2))^2);
        if(TT(i,j) < 0)
            w(TP(i,j)) = w(TP(i,j)) + 1/L;
            w(TP(i,jp)) = w(TP(i,jp)) + 1/L;
        else
            w(TP(i,j)) = w(TP(i,j)) + 1/(2*L);
            w(TP(i,jp)) = w(TP(i,jp)) + 1/(2*L);
        end
    end
end

M = sparse(np,np);
for i = 1:nt
    for j = 1:3
        jp = mod(j,3)+1;
        L = 1;%sqrt((PXY0(TP(i,j),1)-PXY0(TP(i,jp),1))^2 + (PXY0(TP(i,j),2)-PXY0(TP(i,jp),2))^2);
        if(TT(i,j) > -1)
            M(TP(i,j),TP(i,jp)) = M(TP(i,j),TP(i,jp)) + (1/(2*L))/w(TP(i,j));
            M(TP(i,jp),TP(i,j)) = M(TP(i,jp),TP(i,j)) + (1/(2*L))/w(TP(i,jp));
        end
    end
end
for i = Ib
    M(i,:) = 0;
    M(i,i) = 1;
end
for i = Iin
    M(i,i) = M(i,i) - 1;
end
Min = M(Iin,Iin);
Mb = M(Iin,Ib);

PX = zeros(np,1);
PY = zeros(np,1);
PX(Ib) = PXY0(Ib,1);
PY(Ib) = PXY0(Ib,2);

% deformace site
PXbody = PX(Ibody);
PYbody = PY(Ibody);
PX(Ibody) = cos(alpha)*PXbody - sin(alpha)*PYbody;
PY(Ibody) = sin(alpha)*PXbody + cos(alpha)*PYbody;

PX(Iin) = Min\(-Mb*PX(Ib));
PY(Iin) = Min\(-Mb*PY(Ib));

PXY = [PX,PY];

% figure
% hold on
% triplot(TP,PXY0(:,1),PXY0(:,2),'color',[0.7,0.7,0.7]);
% % set(h,'facealpha',0.3);
% triplot(TP,PX,PY,'r')
% axis equal













