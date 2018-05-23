function show(varargin)
    % show - export data for post-processing
    % input parameters: 
    % show <variableName> -i<iteration> -f<format>    
    %
    % formats: vtk (Paraview), txt (plain text)
    % 
    % variable names: depend on the model used
    %
    % example:
    %
    %   show mach pressure -i100 -ftxt
    
    str = 'java -d64 -Xmx8g -jar FlowPro.jar postprocessing ';
    for i = 1:nargin
        str = [str,' ', varargin{i}];
    end

    path = pwd;
    cd(getFlowProPath);
    system(str);
    cd(path)
%     return
    try
        [meshPath, simulPath, outputPath] = getPath;
        vertices = load([outputPath,'vertices.txt']);
        if(size(vertices,2) > 2)
            return;
        end
        elements = nactiMat([meshPath,'elements.txt']);
        elementType = load([meshPath,'elementType.txt']);

        % tvorba triangulace
        ne = length(elements(:,1));
        tri = [];
        for i = 1:ne
            for j = 2:elementType(i)-1
                tri = [tri;elements(i,[1,j,j+1])];
            end
        end    

        for i = 1:nargin
            if(varargin{i}(1) ~= '-')
                variableName = varargin{i};
                quantity = load([outputPath,[variableName,'.txt']]);
                myContour(tri, vertices, quantity, variableName);
            end
        end
    catch
        close all;
        display('Unable to show results in matlab!');
    end
end


function myContour(tri, vertices, Quantity, name)
    figure('name', name, 'color', 'w');
    [h1, h2] = tricontf(vertices(:,1),vertices(:,2),tri+1,Quantity,30);
    % tricontour(tri,PX,PY,Quantity,30)
    set(h2, 'linestyle', 'none');
    box on;
    axis equal;
    osy = [min(vertices(:,1)) max(vertices(:,1)) min(vertices(:,2)) max(vertices(:,2))];
    axis(osy);
    colorbar;
    colormap jet
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