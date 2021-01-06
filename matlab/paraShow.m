function paraShow(varargin)

args = '';
for i = 1:nargin
    args = [args,' ',varargin{i}];
end
args = [args,' -fvtk'];
show(args);
[~, ~, outputPath] = getPath;
system(['start paraview ',outputPath,'results.vtk']);
