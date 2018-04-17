function PXYd = deformation(def)
global profileNURBS profile X A b

% initialization
if(nargin == 0)
    init;
    return
end

profileNURBSdef = profileNURBS;
dx = def(1:profileNURBSdef.nc);
dy = def(profileNURBSdef.nc+[1:profileNURBSdef.nc]);

profileNURBSdef = profileNURBS;
profileNURBSdef.xc = profileNURBS.xc + dx;
profileNURBSdef.yc = profileNURBS.yc + dy;
data = nurbsClosed(profileNURBSdef);
data = data(1:end-1,:);
I = 1:profileNURBSdef.np-1;
b(I,:) = data - profile;

cx = A\b(:,1);
cy = A\b(:,2);

PXY = load('initMesh/vertices.txt');
PXYd = PXY;

ntot = length(cx);
for i = 1:size(PXY,1)
    for j = 1:ntot
        rbf = 1-norm(PXY(i,:)-X(j,:));
        PXYd(i,1) = PXYd(i,1) + cx(j)*rbf;
        PXYd(i,2) = PXYd(i,2) + cy(j)*rbf;
    end
end


function init
global meshPath simulationPath optimisationPath paramFile profileNURBS profile X A b

[geometry, simulation] = loadArgs;
geometryPath = strcat('../../simulations/', geometry, '/');
meshPath = strcat(geometryPath, 'mesh/');
simulationPath = strcat(geometryPath, simulation, '/');
optimisationPath = strcat(geometryPath, simulation, '/optimisation/');
paramFile = [simulationPath,'parameters.txt'];

prof = load([meshPath,'profileNURBS.mat']);
profileNURBS = prof.profile;

d = load([meshPath,'profile.mat']);
profile = d.data;
d = load([meshPath,'box.mat']);
box = d.data;
np = size(profile,1);
nb = size(box,1);
ntot = np + nb;
X = [profile;box];
A = zeros(ntot,ntot);
b = zeros(ntot,2);
for i = 1:ntot
    for j = 1:ntot
        A(i,j) = 1-norm(X(i,:)-X(j,:));
    end
end
