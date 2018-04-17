function initDeform
global A
[geometry, simulation] = loadArgs;
geometryPath = strcat('../../simulations/', geometry, '/');
meshPath = strcat(geometryPath, '/mesh/');
profile = load([meshParh,'profile.mat']);
box = load([meshPath,'box.mat']);
np = size(profile,1);
nb = size(box,1);
nc = np + nb;
X = [profile;box];
A = zeros(nc,nc);
for i = 1:nc
    for j = 1:nc
        A(i,j) = 1-norm(X(i,:)-X(j,:));
    end
end