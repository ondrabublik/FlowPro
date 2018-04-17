function deform
global meshPath profileNURBS profile

deformation;
dx = zeros(profileNURBS.nc,1);
dy = zeros(profileNURBS.nc,1);
dy(1) = dy(1)-0.05;
PXYd = deformation([dx;dy]);

TP = loadTP([meshPath,'elements.txt'])+1;
PXY = load([meshPath,'vertices.txt']);
figure
hold on
triplot(TP(:,1:3),PXY(:,1),PXY(:,2),'color',[0.7,0.7,0.7])
triplot(TP(:,1:3),PXYd(:,1),PXYd(:,2))
plot(profile(:,1),profile(:,2),'.r');
axis equal

function TP = loadTP(str)
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
