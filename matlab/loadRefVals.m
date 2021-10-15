function  refVals = loadRefVals(simulPath)

if nargin == 0
    [~, simulPath, ~] = getPath;
end
    
filePath = strcat(simulPath, 'referenceValues.txt');

refVals = loadPropertiesFile(filePath);
