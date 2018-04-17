function showDual

    [geometry, simulation] = loadArgs;
    geometryPath = strcat('../../simulations/', geometry, '/');
    meshPath = strcat(geometryPath, 'mesh/');
    optimisationPath = strcat(geometryPath, simulation, '/optimisation/');

    W = load([optimisationPath,'lambda.txt']);
    PXY = load(strcat(meshPath, 'vertices.txt'));
    typ = firstDigit(load(strcat(meshPath, 'elementType.txt')));
    TP = loadTP(strcat(meshPath, 'elements.txt'))+1; 
    
    PX = PXY(:,1);
    PY = PXY(:,2);    
    
    nPoints = size(PX, 1);    

    osy = [min(PX) max(PX) min(PY) max(PY)];

    Wb = prevedNaUzly(TP,W,typ,nPoints);
    
    % tvorba triangulace
    nq = length(TP(:,1));
    tri = [];
    for i = 1:nq
        for j = 2:typ(i)-1
            tri = [tri;TP(i,[1,j,j+1])];
        end
    end
    
    try
        for i = 1:10
            myContour(tri, PX, PY, Wb(:,i), ['dual ',num2str(i)], osy);
        end
    catch
    end
end

function myContour(tri, PX, PY, Quantity, name, osy)
    figure('name', name, 'color', 'w');
%     Quantity(Quantity > 20) = 20;
%     Quantity(Quantity < -15) = -15;
    [h1, h2] = tricontf(PX,PY,tri,Quantity,50);
%     tricontour(tri,PX,PY,Quantity,30)
    set(h2, 'linestyle', 'none');
    box on;
    axis equal;
    axis(osy);
    colorbar;
    colormap jet;
end

function Wb = prevedNaUzly(tri,W,typ,np)
    nt = length(W(:,1));
    nr = length(W(1,:));
    Wb = zeros(np,nr);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:typ(i)
            Wb(tri(i,j),:) = Wb(tri(i,j),:) + W(i,:);
            poc(tri(i,j)) = poc(tri(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i,:) = Wb(i,:)/poc(i);
    end
end


function color = colors(domIdx)
    switch mod(domIdx, 10)
        case 0
            color = 'blue';
        case 1
            color = 'magenta';
        case 2
            color = 'green';
        case 3
            color = 'red';
        case 4
            color = 'cyan';
        case 5
            color = 'black';
        case 6
            color = [188 143 143] / 255; % Rosy Brown        
        case 7
            color = [176 48 96] / 255; % maroon
        case 8
            color = [255 215 0] / 255; % yellow
        case 9
            color = [255 140 0] / 255;
        %             color = [210 105 30] / 255; % chocolate brown
    end
end

function t = firstDigit(typ)
    t = typ;
    for i = 1:length(typ)
        while(typ(i) > 0)
            t(i) = typ(i);
            typ(i) = fix(typ(i)/10);
        end
    end
end

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
end





