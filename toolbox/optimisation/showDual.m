function showDual

    [geometry, simulation] = loadArgs;
    geometryPath = strcat('../../simulations/', geometry, '/');
    meshPath = strcat(geometryPath, 'mesh/');
    optimisationPath = strcat(geometryPath, simulation, '/optimisation/');

    W = load([optimisationPath,'lambda.txt']);
    PXY = load(strcat(meshPath, 'vertices.txt'));
    TP = load(strcat(meshPath, 'elements.txt')); 
    typ = firstDigit(load(strcat(meshPath, 'elementType.txt')));
    
    PX = PXY(:,1);
    PY = PXY(:,2); 
    TP = TP + 1;    
    
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

function Wb = prevedNaUzly(TP,W,typ,np)
    
    nt = length(W(:,1));
    nr = length(W(1,:));
    Wb = zeros(np,nr);
    poc = zeros(np,1);
    for i = 1:nt
        for j = 1:typ(i)
            Wb(TP(i,j),:) = Wb(TP(i,j),:) + W(i,:);
            poc(TP(i,j)) = poc(TP(i,j)) + 1;
        end
    end
    
    for i = 1:np
        Wb(i,:) = Wb(i,:)/poc(i);
    end
end

function plotMesh(PX, PY, TP, typ, name, osy)
    figure('name', name, 'color', 'w');
%     triplot(TP, PX, PY, 'color', [0 0 0]);    
    hold on;
    for i = 1 : size(TP, 1)
        indexes = TP(i, [1:typ(i), 1]);
        plot(PX(indexes), PY(indexes), 'color', [0 0 0]);
    end
    box on;
    axis equal;
    axis(osy);
end

function countDomSize(part, nDoms)
    domSize = zeros(nDoms, 1);
    for i = 1 : length(part)
        domIdx = part(i) + 1;
        domSize(domIdx) = domSize(domIdx) + 1;
    end
    
    fprintf('size of domains: \n');
    for j = 1 : nDoms
        fprintf('%d: %d\n', j - 1, domSize(j));  % colors(j - 1)
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

function plotMeshMetis(PX, PY, TP, typ, part)
    figure;
    hold on;
    for i = 1 : size(TP, 1)
        indexes = TP(i, [1:typ(i), 1]);

        plot(PX(indexes), PY(indexes), 'color', colors(part(i)));
    end
    box on;
    axis equal;
    
    % triplot(TP, PX, PY, 'color', [0 0 0]);
end

function [W,PXY] = loadData(simulationPath, str)
    fid = fopen(str,'r');
    while 1
        line = fgetl(fid);
        if ~ischar(line)
            break;
        end
        
    end
    fclose(fid);
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







