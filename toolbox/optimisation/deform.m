function PXY = deform(alpha)
global optimisationPath

PXY0 = load([optimisationPath,'vertices.txt']);
PXY = zeros(size(PXY0));
PXY(:,1) = cos(alpha)*PXY0(:,1) - sin(alpha)*PXY0(:,2);
PXY(:,2) = sin(alpha)*PXY0(:,1) + cos(alpha)*PXY0(:,2);

r = 0.8;
R = 2;
rad = sqrt((PXY(:,1)-0.5).^2 + PXY(:,2).^2);
b = (rad-r)/(R-r);
b(rad > R) = 1;
b(rad < r) = 0;

PXY(:,1) = b.*PXY0(:,1) + (1-b).*PXY(:,1);
PXY(:,2) = b.*PXY0(:,2) + (1-b).*PXY(:,2);
