%param  Print the state of computation when the data was saved.

function state

[geometry, simulation] = loadArgs;
flowProPath = getFlowProPath;

path = sprintf('%s/simulations/%s/', flowProPath, geometry);
if ~exist(path, 'dir')
    error('Geometry %s does not exist.', geometry);
end
path = sprintf('%s%s/', path, simulation);
if ~exist(path, 'dir')
    error('Simulation %s does not exist.', simulation);
end

path = strcat(path, 'state.txt');

type(path);