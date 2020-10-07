% getPath   Return mesh, simulation and output paths.
% [meshPath, simulPath, outputPath] = getPath returns path of the mesh,
% simulation and outpud folders for the current simulation.
% [meshPath, simulPath, outputPath] = getPath(geometry, simulation) returns
% path of the mesh, simulation nad outpud folders for the simulation
% specified by the name of the geometry and simulation.

function [meshPath, simulPath, outputPath] = getPath(geometry, simulation)
flowProPath = getFlowProPath;

if nargin == 0
    [geometry, simulation] = loadArgs;
end

meshPath = sprintf('%s/simulations/%s/%s/',flowProPath, geometry, 'mesh');
simulPath = sprintf('%s/simulations/%s/%s/',flowProPath, geometry, simulation);
outputPath = sprintf('%s/simulations/%s/%s/%s/',flowProPath, geometry, simulation, 'output');

if ~exist(meshPath, 'dir')
    warning('Geometry %s does not exist.', geometry);
end

if ~exist(simulPath, 'dir')
    warning('Simulation %s does not exist.', simulation);
%     simulPath = [];
end

% if ~exist(outputPath, 'dir')
%     warning('Output directory does not exist.');
%     outputPath = [];
% end
