function mshow(varargin)

[meshPath, simulPath, outputPath] = getPath;

args = '';
for i = 1:nargin
    quantityName = varargin{i};
    
    if strcmp(lower(quantityName),'mesh') || strcmp(lower(quantityName),'bodies') || strcmp(lower(quantityName),'residuum')
        continue
    end
    
    filePath = [outputPath, quantityName, '.txt'];
    if ~exist(filePath, 'file')    
        args = [args,' ', quantityName];
    end
end

if args ~= ''
    show(args);
end

par = loadParam(simulPath);

if par.dimension ~= 2
    error('mshow supports only 2D view');
end
elements = dlmread([meshPath, 'elements.txt'], ' ');
elementType = load([meshPath, 'elementType.txt']);

tri = convert2Triangular(elements, elementType);

for k = 1 : nargin
    q = varargin{k};

    if q(1) == '-'
        warning('cannot evaluate option ''%s'': switches are not allowed here', q)
        continue
    end
    
    switch lower(q)
        case 'mesh'
            showMesh;
            continue        
        case 'bodies'
            showBodies;
            continue        
        case 'residuum'
            showResiduum;
            continue
    end
    
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
        figure('name', quantityName, 'color', 'w');
        quiver(vertices(:,1), vertices(:,2), quantity(:,1), quantity(:,2));
        axis equal
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

function showMesh(xy, elems, neigh, type)
if nargin == 0
    [meshPath, ~, ~] = getPath;    

    xy = dlmread([meshPath, 'vertices.txt']);    
    elems = dlmread([meshPath, 'elements.txt'])+1;
    neigh = dlmread([meshPath, 'neighbors.txt']);
    type = dlmread([meshPath, 'elementType.txt']);

    fprintf('sit ma %d elementu\n', length(elems));
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
    end
end

axis equal
box on

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

function x = reorganise(x)
    n = length(x);
    for i = 2:n
        if(x(i) < x(i-1))
            x(i:end) = x(i:end) + x(i-1)-x(i);
        end
    end
end
end
