function meshInfo
[meshPath, simulPath, outputPath] = getPath;
PXY = importdata(strcat(meshPath, 'vertices.txt'));
typ = firstDigit(importdata(strcat(meshPath, 'elementType.txt')));
display(['Number of elements: ',num2str(length(PXY(:,1)))]);
TP = nactiMat(strcat(meshPath, 'elements.txt'));
tri = makeTri(TP,typ);

if(tri ~= -1)
    figure
    triplot(tri(:,1:3),PXY(:,1),PXY(:,2))
    axis equal
end

function TP = nactiMat(str)
fid = fopen(str,'r');
n = 0;
max = 0;
while(feof(fid) == 0)
    line = fgetl(fid);
    v = str2num(line);
    if(length(v) > max)
        max = length(v);
    end
    n = n + 1;
end
fclose(fid);
TP = zeros(n,max);
fid = fopen(str,'r');
n = 1;
while(feof(fid) == 0)
    line = fgetl(fid);
    v = str2num(line);
    TP(n,1:length(v)) = v;
    n = n + 1;
end
fclose(fid);

function tri = makeTri(TP,typ)
nt = 0;
for i = 1:length(typ)
    if(typ(i) == 3)
        nt = nt + 1;
    end
    if(typ(i) == 4)
        nt = nt + 2;
    end
end
if(nt == 0)
    tri = -1;
    return;
end
tri = zeros(nt,3);
s = 1;
for i = 1:length(typ)
    if(typ(i) == 3)
        tri(s,:) = TP(i,1:3)+1;
        s = s + 1;
    end
    if(typ(i) == 4)
        tri(s,:) = TP(i,1:3)+1;
        tri(s+1,:) = TP(i,[1,3,4])+1;
        s = s + 2;
    end
end



