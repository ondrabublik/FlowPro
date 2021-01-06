%   loadArgs   Load geometry and simulaton name from 'args.txt' file.
%   loadArgs returns two strings respectively with the geometry and
%   simulation name.

function [geometry, simulation] = loadArgs
    flowProPath = getFlowProPath;
    
    str = fileread([flowProPath,'/args.txt']);
    
    % read first line from file and get rid of \n and \r characters
%     line = strsplit(str, '\n')
%     line = strsplit(line{1}, '\r');
    
    % split the line into an array of strings, which are separated by space
%     strArray = strsplit(line{1}, ' ');
    strArray = strsplit(str, ' ');
    
    geometry = strArray{1};
    simulation = strArray{2};
end
