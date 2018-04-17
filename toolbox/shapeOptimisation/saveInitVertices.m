function saveInitVertices

[geometry, simulation] = loadArgs;
geometryPath = strcat('../../simulations/', geometry, '/');
meshPath = strcat(geometryPath, 'mesh/');
copyfile([meshPath,'vertices.txt'],'initMesh/vertices.txt');