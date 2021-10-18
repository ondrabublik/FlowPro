%   loadArgs   Load geometry and simulaton name from 'args.txt' file.
%   loadArgs returns two strings respectively with the geometry and
%   simulation name.

function [geometry, simulation] = loadArgs
    flowProPath = getFlowProPath;
        
    file = fopen([flowProPath,'/args.txt']);
    line = fgetl(file);
    fclose(file);
    
    strArray = strsplit(line, ' ');
    
    geometry = strArray{1};
    simulation = strArray{2};
end
