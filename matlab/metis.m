function metis(nDoms)

if ischar(nDoms)
    nDoms = str2double(nDoms);
end

meshFileName = 'meshfile.txt';

metisCommand = sprintf('mpmetis %s %d', meshFileName, nDoms);
metisOutputFile = sprintf('%s.epart.%d', meshFileName, nDoms);
metisOutputFile2 = sprintf('%s.npart.%d', meshFileName, nDoms);

[meshPath, simulationPath] = getPath;
TP = dlmread(strcat(meshPath, 'elements.txt'),' ') + 1;
typ = dlmread(strcat(meshPath, 'elementType.txt'));
filePath = strcat(simulationPath, 'part.txt');
if (nDoms == 1)
    fileID = fopen(filePath, 'w');
    for i = 1 : size(TP, 1)
        fprintf(fileID, '%d\n', 0);
    end
    fclose(fileID);
else 
   
    createMeshFile(TP, typ, meshFileName);
    try
        system(metisCommand);
        copyfile(metisOutputFile, filePath);
        delete(meshFileName, metisOutputFile, metisOutputFile2);
    catch
        fclose('all');
        delete(meshFileName);
        error('Metis not suported! ');
    end
end
%fprintf(1, 'the partition of the mesh has been saved into %s\n', filePath);

% PXY = importdata(strcat(meshPath, 'PXY.txt'));
% PX = PXY(:,1);
% PY = PXY(:,2);
% countDomSize(part, nDoms);
% plotMesh(PX, PY, TP, typ, part);
end

function createMeshFile(TP, typ, fileName)
    fclose('all');
    fid = fopen(fileName, 'w');
    fprintf(fid, '%d\n', size(TP, 1));
    for i = 1 : size(TP, 1)            
        for j = 1 : typ(i)
            fprintf(fid, '%d ', TP(i,j));                
        end
        fprintf(fid, '\n');
    end
    
    fclose(fid);
end
