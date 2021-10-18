function  refVals = loadRefVals(simulationPath)
% refvals   Return reference values of the simulation in a structure.

if nargin == 0
    [~, simulationPath, ~] = getPath;
end
    
filePath = strcat(simulationPath, 'referenceValues.txt');

refVals = loadPropertiesFile(filePath);
