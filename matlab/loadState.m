function  state = loadState(simulPath)

if nargin == 0
    [~, simulPath, ~] = getPath;
end
    
filePath = strcat(simulPath, 'state.txt');

state = loadPropertiesFile(filePath);
