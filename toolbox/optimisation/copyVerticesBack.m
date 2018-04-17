function copyVerticesBack
global meshPath optimisationPath

copyfile('initMesh/vertices.txt',[meshPath,'vertices.txt']);
copyfile('initMesh/vertices.txt',[optimisationPath,'vertices.txt']);