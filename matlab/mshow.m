function mshow(varargin)

args = '';
for i = 1:nargin
    args = [args,' ', varargin{i}];
end
show(args);

[meshPath, simulPath, outputPath] = getPath;

par = loadParam(simulPath);
vertices = load([outputPath,'vertices.txt']);
if par.dimension ~= 2 && size(vertices,2) ~= 2
    error('mshow supports only 2D view');
end
elements = dlmread([meshPath, 'elements.txt'], ' ');
elementType = load([meshPath, 'elementType.txt']);

tri = convert2Triangular(elements, elementType);

show(varargin{:});

for k = 1 : nargin
    q = varargin{k};

    if q(1) == '-'
        warning('cannot evaluate option ''%s'': switches are not allowed here', q)
        continue
    end
    
    plotMe(q)
end

function plotMe(quantityName)
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
