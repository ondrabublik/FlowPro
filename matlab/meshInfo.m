function meshInfo
[meshPath, simulPath, outputPath] = getPath;
PXY = importdata(strcat(meshPath, 'vertices.txt'));
typ = firstDigit(importdata(strcat(meshPath, 'elementType.txt')));
display(['Number of elements: ',num2str(length(PXY(:,1)))]);

bf = importdata(strcat(meshPath, 'blendingFunctions.txt'));
tri = importdata(strcat(meshPath, 'elements.txt'));

trimesh(tri+1,PXY(:,1),PXY(:,2),bf)