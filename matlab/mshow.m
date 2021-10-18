function mshow(varargin)

[meshPath, simulPath, outputPath] = getPath;

args = '';
for i = 1:nargin
    quantityName = varargin{i};
    
    if strcmpi(quantityName,'mesh') || strcmpi(quantityName,'meshALE') || strcmpi(quantityName,'order') || ...
       strcmpi(quantityName,'bodies') || strcmpi(quantityName,'residuum') || ...
       strcmpi(quantityName,'av') || strcmpi(quantityName,'meshvelo') || strcmpi(quantityName,'y') || ...
       strcmpi(quantityName,'bf')
        continue
    end
    
%     filePath = [outputPath, quantityName, '.txt'];
%     if ~exist(filePath, 'file')    
    args = [args,' ', quantityName];
%     end
end

if ~isempty(args)
    show(args);
end

par = loadParam(simulPath);

polyOrder = 1;
for k = 1 : nargin
    q = varargin{k};

    if q(1) == '-'
        if q(2) == 'p'
            polyOrder = q(3:end);
        end
%         warning('cannot evaluate option ''%s'': switches are not allowed here', q)
        continue
    end

    switch lower(q)
        case 'mesh'
            showMesh;
            continue  
        case 'meshale'
            showMeshALE;
            continue  
        case 'order'
            showMesh;
            showOrder;
            continue
        case 'bodies'
            showBodies;
            continue        
        case 'residuum'
            showResiduum;
            continue
        case 'av'
            showArtificialViscosity;
            continue
        case 'meshvelo'
            showMeshVelo;
            continue  
        case 'y'
            showWallDistance;
            continue
        case 'bf'
            showBlendingFunctions;
            continue
    end

    if par.dimension ~= 2
        error('mshow supports only 2D view');
    end
    
    if polyOrder == 1
        elements = dlmread([meshPath, 'elements.txt'], ' ');
        elementType = load([meshPath, 'elementType.txt']);
    elseif polyOrder > 1
        elements = dlmread([outputPath, 'elements.txt'], ' ');
        elementType = 3 * ones(size(elements, 1), 1);
    end
    tri = convert2Triangular(elements, elementType);
    
    plotMe(q)
end

function plotMe(quantityName)
    vertices = load([outputPath,'vertices.txt']);
    filePath = [outputPath, quantityName, '.txt'];
    
    if ~exist(filePath, 'file')
        warning('cannot find file %s', filePath)
    end
    
    quantity = load(filePath);
    m = size(quantity, 2);
    if m == 1
        myContour(tri, vertices, quantity, quantityName);
    elseif m > 2
        magnitude = zeros(size(quantity, 1), 1);
        for j = 1 : m
            magnitude = magnitude + quantity(:,j).^2;
        end
        magnitude = sqrt(magnitude);
        myContour(tri, vertices, magnitude, quantityName);
%         figure('name', quantityName, 'color', 'w');
%         quiver(vertices(:,1), vertices(:,2), quantity(:,1), quantity(:,2));
%         axis equal
    end
end

end

function myContour(tri, vertices, Quantity, name)
    figure('name', name, 'color', 'w');
    [~, h2] = tricontf(vertices(:,1),vertices(:,2),tri+1,Quantity,30);
    % tricontour(tri,PX,PY,Quantity,30)
    set(h2, 'linestyle', 'none');
    box on;
    axis equal;
    osy = [min(vertices(:,1)) max(vertices(:,1)) min(vertices(:,2)) max(vertices(:,2))];
    axis(osy);
    colorbar;
    colormap jet
end

% splits any quadriliteral elements into two triangles for the
% visualisation purposes
function tri = convert2Triangular(elements, elementType)

    nElems = size(elements, 1);
    
    if nargin > 1
        squares = elements(firstDigit(elementType) == 4, :);        
    else
        squares = elements(elements(:,4) ~= 0, :);
    end
    
    tri = zeros(nElems + length(squares), 3);
    tri(1:nElems, :) = elements(:,1:3);
    for i = 1 : length(squares)
        tri(nElems+i, :) = [squares(i,1) squares(i,3) squares(i,4)];
    end    
end

function showMesh(xy, elems, neigh, type)
if nargin == 0
    [meshPath, ~, ~] = getPath;    

    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    neigh = dlmread([meshPath, 'neighbors.txt']);
    type = dlmread([meshPath, 'elementType.txt']);

    fprintf('Mesh has %d elements and %d vertices.\n', length(elems), length(xy));
end

x = xy(:,1);
y = xy(:,2);

tri = convert2Triangular(elems, type);

figure('color','w','name', 'mesh')
triplot(tri, xy(:,1), xy(:,2), 'k', 'linewidth', .1);

hold on
color = 'bgrm';
linewidth = 2;
for i = 1 : size(elems,1)
    for j = 1 : type(i)        
        if neigh(i,j) < 0
            jp = mod(j, type(i)) + 1;
            ind = [elems(i,j),elems(i,jp)];

            plot(x(ind),y(ind),'color',color(-neigh(i,j)),'linewidth',linewidth);
        end
        
%         xx = x(elems(i,j));
%         yy = y(elems(i,j));
%         xlims = [-30 -20];
%         ylims = [-5 5];
%         if xx >= xlims(1) && xx <= xlims(2) && yy >= ylims(1) && yy < ylims(2)
%             text(xx, yy, sprintf('%d', elems(i,j)))            
%         end        
    end
        
%     inds = elems(i,1:3);
%     if max(x(inds)) >= 0.1-eps
%         xx = mean(x(inds));    
%         yy = mean(y(inds));
%         text(xx, yy, sprintf('%d', i));
%     end

%     xlims = [-30 -5];
%     ylims = [-5 5];
%     xs = mean(x(elems(i,:)));
%     ys = mean(y(elems(i,:)));
%     if xs >= xlims(1) && xs <= xlims(2) && ys >= ylims(1) && ys < ylims(2)
%         text(xs, ys, sprintf('%d', i))
%     end
end

axis equal
box on

end

function showMeshALE(xy, elems, neigh, type)
    if nargin == 0
        [meshPath, ~, ~] = getPath;    

        xy = dlmread([meshPath, 'vertices.txt']);    
        elems = dlmread([meshPath, 'elements.txt'])+1;
        neigh = dlmread([meshPath, 'neighbors.txt']);
        neighALE = dlmread([meshPath, 'neighborsALE.txt']);
        type = dlmread([meshPath, 'elementType.txt']);

        fprintf('Mesh has %d elements and %d vertices.\n', length(elems), length(xy));
    end

    x = xy(:,1);
    y = xy(:,2);

    tri = convert2Triangular(elems, type);

    figure('color','w','name', 'mesh')
    triplot(tri, xy(:,1), xy(:,2), 'k', 'linewidth', .1);

    hold on
    color = 'bgrmycbgrmyc';
    linewidth = 2;
    
    for i = 1 : size(elems,1)
        for j = 1 : type(i)        
            if neigh(i,j) < 0
                jp = mod(j, type(i)) + 1;
                ind = [elems(i,j),elems(i,jp)];                                
                
                plot(x(ind),y(ind),'color',color(neighALE(i,j)),'linewidth',linewidth);
            end
        end
    end
    axis equal
    box on
end

function showOrder
    [meshPath, simulPath, ~] = getPath;

    xy = dlmread([meshPath, 'vertices.txt']);
    xy = xy(:, 1:2);
    elems = dlmread([meshPath, 'elements.txt'])+1;
    type = dlmread([meshPath, 'elementType.txt']);
    nElem = length(elems(:,1));
    xys = zeros(nElem,2);
    
    filePath = [simulPath, 'order.txt'];
    if exist(filePath, 'file')
        order = dlmread([simulPath, 'order.txt']);
    else
        par = loadParam();
        order = par.order * ones(size(elems, 1), 1);
        fprintf('File %s does not exist - global accuracy is %d.\n', filePath, par.order);
    end
    
    for i = 1:nElem
        for j = 1:type(i)
            xys(i,:) = xys(i,:) + xy(elems(i,j),:);
        end
        xys(i,:) = xys(i,:)/type(i);
    end
    
    col = 'brgmckybrgmcky';
    hold on
    for i = 1:10
        x = xys(order == i,1);
        y = xys(order == i,2);
        plot(x,y,'marker','.','color',col(i),'linestyle','none');
        
        n = length(x);
        if n > 0
            fprintf('%d elems of order %d\n', n, i);
        end
    end
end

function showBodies

[meshPath, ~, ~] = getPath;    

xy = dlmread([meshPath, 'vertices.txt']);    
ale = dlmread([meshPath, 'boundaryTypeALE.txt']);

color = 'brkmcgkybrkmcgkybrkmcgky';

figure('color','w','name', 'bodies')
hold on
linewidth = 2;
for i = 1 : length(ale)
    inds = ale(i,2:3)+1;
    plot(xy(inds,1), xy(inds,2), 'color', color(ale(i,1)), 'linewidth', linewidth)
end
axis equal
box on

end

function showResiduum

[~, simulPath, ~] = getPath;
residuum = load([simulPath,'residuum.txt']);

linewidth = 2;

figure('color','w','name','residuum')
semilogy(residuum(:,1),'color','k','linewidth',linewidth)
grid on
box on
ylabel('resid')
xlabel('iteration')

figure('color','w','name','residuum')
semilogy(residuum(:,2),residuum(:,1),'color','k','linewidth',linewidth)
grid on
box on
ylabel('resid')
xlabel('time')

figure('color','w','name','residuum')
semilogy(reorganise(residuum(:,3)/1000),residuum(:,1),'color','k','linewidth',linewidth)
grid on
box on
ylabel('resid')
xlabel('CPU [s]')

end

function showArtificialViscosity
    [meshPath, simulPath, ~] = getPath;   
    
    av = load([simulPath,'artificialViscosity.txt']);
    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    type = firstDigit(dlmread([meshPath, 'elementType.txt']));
    tri = convert2Triangular(elems, type);
    
    
    for m = 1:length(av(1,:))
        vav = zeros(size(xy,1),1);
        pom = zeros(size(xy,1),1);
        for i = 1:size(elems,1)
            for j = 1:type(i)
                vav(elems(i,j)) = vav(elems(i,j)) + av(i,m);
                pom(elems(i,j)) = pom(elems(i,j)) + 1;
            end
        end
        vav = vav./pom;

        figure('color', 'w');
        [~, h2] = tricontf(xy(:,1),xy(:,2),tri,vav,30);
        % tricontour(tri,PX,PY,Quantity,30)
        set(h2, 'linestyle', 'none');
        box on;
        axis equal;
        osy = [min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
        axis(osy);
        colorbar;
        colormap jet
    end
end

function showMeshVelo
    [meshPath, simulPath, ~] = getPath;   
    
    uxy = load([simulPath,'UXY.txt']);
    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    type = firstDigit(dlmread([meshPath, 'elementType.txt']));
    tri = convert2Triangular(elems, type);
    
    
    for m = 1:length(uxy(1,:))
        try
            figure('color', 'w');
            [~, h2] = tricontf(xy(:,1),xy(:,2),tri,uxy(:,m),30);
    %         tricontour(tri,xy(:,1),xy(:,2),uxy(:,m),30)
            set(h2, 'linestyle', 'none');
            box on;
            axis equal;
            osy = [min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
            axis(osy);
            colorbar;
            colormap jet
        catch
        end
    end
end

function showWallDistance
    [meshPath, ~, ~] = getPath;   
    
    dist = load([meshPath,'wallDistance.txt']);
    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    type = firstDigit(dlmread([meshPath, 'elementType.txt']));
    tri = convert2Triangular(elems, type);
    
    figure('color', 'w');
    [~, h2] = tricontf(xy(:,1),xy(:,2),tri,dist,30);
    % tricontour(tri,PX,PY,Quantity,30)
    set(h2, 'linestyle', 'none');
    box on;
    axis equal;
    osy = [min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
    axis(osy);
    colorbar;
    colormap jet
end

function x = reorganise(x)
    n = length(x);
    for i = 2:n
        if(x(i) < x(i-1))
            x(i:end) = x(i:end) + x(i-1)-x(i);
        end
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

function showBlendingFunctions
    [meshPath, ~, ~] = getPath;   
    
    bfs = load([meshPath,'blendingFunctions.txt']);
    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    type = firstDigit(dlmread([meshPath, 'elementType.txt']));
    tri = convert2Triangular(elems, type);
    
    for i = 1:length(bfs(1,:))
        figure('color', 'w');
        [~, h2] = tricontf(xy(:,1),xy(:,2),tri,bfs(:,i),30);
%         trimesh(tri,xy(:,1),xy(:,2),bfs(:,i));
        set(h2, 'linestyle', 'none');
        box on;
        axis equal;
        osy = [min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))];
        axis(osy);
        colorbar;
        colormap jet
    end
end
