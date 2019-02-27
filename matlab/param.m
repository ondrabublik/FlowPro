%param  Open text file with parameters for the problem.

function param

[geometry, simulation] = loadArgs;
flowProPath = getFlowProPath;

path = sprintf('%s/simulations/%s/', flowProPath, geometry);
if ~exist(path, 'dir')
    error('Geometry %s does not exist.', geometry);
end
path = sprintf('%s%s/', path, simulation);
if ~exist(path, 'dir')
    answer = input(sprintf('Simulation %s does not exist, do you want to creat it [Y/n]? ', ...
             simulation), 's');
    
    switch lower(answer)
        case {'y', 'yes', ''}
            actualPath = pwd;
            cd(flowProPath)
            !java -jar FlowProManager.jar createparamfile
            cd(actualPath)
            edit([path,'parameters.txt']);
        otherwise
            return;
    end
else
    edit([path,'parameters.txt']);
end

[~, simulPath, ~] = getPath;
cd(simulPath);