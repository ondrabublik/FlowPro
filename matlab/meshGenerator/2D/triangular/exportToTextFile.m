%path = '../galerkin/meshes/NACA0012_1_25deg/';
path = 'Desktop/martin0_1';

load geometrie
PX = geometrie{1};
PY = geometrie{2};
T = geometrie{3};
Neighbours = geometrie{4};

mkdir(path);
dlmwrite(strcat(path,'/points.txt'), [PX, PY], 'delimiter', ' ', 'precision', 16);
dlmwrite(strcat(path,'/triangles.txt'), T(:,1:3), 'delimiter', ' ', 'precision', 16);
dlmwrite(strcat(path,'/neighbours.txt'), Neighbours(:,1:3), 'delimiter', ' ', 'precision', 16);

T = T + 1;
T = T(:,1:3);
figure(1);
triplot(T, PX, PY);
axis equal
% axis([0 3 0 1]);
