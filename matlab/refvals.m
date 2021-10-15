function vals = refvals(simulationPath)
% refvals   Return reference values of the simulation in a structure.

if nargin == 0
    [~, simulationPath, ~] = getPath;
end
    
filePath = strcat(simulationPath, 'referenceValuse.txt');

vals = loadPropertiesFile(filePath);
