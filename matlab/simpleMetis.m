function simpleMetis(nDoms,char)
    if ischar(nDoms)
        nDoms = str2double(nDoms);
    end

    [meshPath, simulation] = getPath();
    PXY = importdata(strcat(meshPath, 'vertices.txt'));
    TP = nactiMat(strcat(meshPath, 'elements.txt')) + 1;
    typ = firstDigit(importdata(strcat(meshPath, 'elementType.txt')));
    filePath = strcat(simulation, 'part.txt');

    if(char ~= 'k')
        if(char == 'x')
            PX = PXY(:,1);
        elseif(char == 'y')
            PX = PXY(:,2);
        else
            PX = PXY(:,3);
        end
        nt = length(TP(:,1));
        map = zeros(nt,1);
        if (nDoms ~= 1)
            PXs = zeros(nt,1);
            for i = 1:nt
                vn = verticesNumber(typ(i));
                for j = 1:vn
                    PXs(i) = PXs(i) + PX(TP(i,j));
                end
                PXs(i) = PXs(i)/vn;
            end
            [Y,I] = sort(PXs);
            mapInd = 0;
            s = 0;
            for i = 1:nt
                map(I(i)) = mapInd;
                s = s + 1;
                if(s > nt/nDoms)
                    mapInd = mapInd + 1;
                    s = 0;
                end
            end
        end
    else
        nt = length(TP(:,1));
        map = zeros(nt,1);
        PXs = zeros(nt,1);
        PYs = zeros(nt,1);
        nd = nt/nDoms;
        for i = 1:nt
            vn = verticesNumber(typ(i));
            for j = 1:vn
                PXs(i) = PXs(i) + PXY(TP(i,j),1);
                PYs(i) = PYs(i) + PXY(TP(i,j),2);
            end
            PXs(i) = PXs(i)/vn;
            PYs(i) = PYs(i)/vn;
        end
        xs = sum(PXs)/nt;
        ys = sum(PYs)/nt;
        R = sqrt((PXs-xs).^2 + (PYs-ys).^2);
        [Y,I] = sort(R);
        map(I([1:nd]')) = nDoms-1;

        alfa = acos((PXs-xs)./sqrt((PXs-xs).^2 + (PYs-ys).^2));
        alfa((PYs-ys) < 0) = 2*pi - alfa((PYs-ys) < 0);
        [Y,I] = sort(alfa);
        mapInd = 0;
        s = 0;
        for i = 1:nt
            if(map(I(i)) == 0) 
                map(I(i)) = mapInd;
                s = s + 1;
                if(s > nd)
                    mapInd = mapInd + 1;
                    s = 0;
                end
            end
        end
    
    end

%     plotMeshMetis(PXY(:,1), PXY(:,2), TP, typ, map)
    fileID = fopen(filePath, 'w');
    for i = 1 : size(TP, 1)
        fprintf(fileID, '%d\n', map(i));
    end
    fclose(fileID);
end


function n = verticesNumber(typ)
    switch(typ)
        case 3
            n = 3;
        case 4
            n = 4;
        case 5
            n = 4;
        case 6
            n = 8;
        case 7
            n = 5;
    end
end


function plotMeshMetis(PX, PY, TP, typ, part)
    tri = convert2Triangular(TP, typ);
    figure;
    hold on;
    colors = 'rbgkmycrbgkmyc';
    triplot(tri, PX, PY, 'color', [0.8,0.8,0.8]);
    nt = length(TP(:,1));
    PXs = zeros(nt,1);
    PYs = zeros(nt,1);
    for i = 1:nt
        PXs(i) = sum(PX(TP(i,1:typ(i))))/typ(i);
        PYs(i) = sum(PY(TP(i,1:typ(i))))/typ(i);
    end
    part = part+1;
    for i = 1 : max(part)
        indexes = 1:nt;
        indexes = indexes(part == i);
        plot(PXs(indexes), PYs(indexes), 'marker','.','color', colors(i),'markerfacecolor',colors(i),'linestyle','none');
    end
    box on;
    axis equal;
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
end

function tri = convert2Triangular(elements, elementType)

    nElems = size(elements, 1);
    
    if nargin > 1
        squares = elements(elementType == 4, :);        
    else
        squares = elements(elements(:,4) ~= 0, :);
    end
    
    tri = zeros(nElems + length(squares), 3);
    tri(1:nElems, :) = elements(:,1:3);
    for i = 1 : length(squares)
        tri(nElems+i, :) = [squares(i,1) squares(i,3) squares(i,4)];
    end    
end