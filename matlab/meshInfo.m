function meshInfo
[meshPath, simulPath, outputPath] = getPath;
PXY = importdata(strcat(meshPath, 'vertices.txt'));
typ = firstDigit(importdata(strcat(meshPath, 'elementType.txt')));
display(['Number of elements: ',num2str(length(PXY(:,1)))]);
tri = importdata(strcat(meshPath, 'elements.txt'));

figure
triplot(tri+1,PXY(:,1),PXY(:,2))
axis equal