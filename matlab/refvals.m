function map = refvals(geometry, simulation)
% refvals   Return reference values for the simulation.

if nargin == 0
    [geomPath, simulPath] = getPath;
else 
    [geomPath, simulPath] = getPath(geometry, simulation);
end

path = strcat(simulPath, 'referenceValues.txt');
fid = fopen(path);

keys = {};
vals = [];

tline = fgets(fid);
while ischar(tline)
    strarray = strsplit(tline, '#'); % get rid of comments
%     str = strcat(strarray{1}, ';');  % append simicolon to suppress the output
%     eval(str);
    
    strarray = strsplit(strarray{1}, '=');         
    
    if length(strarray) == 2
        keys = [keys, {strarray{1}}];
        vals = [vals, str2double(strarray{2})];
    end
    
    tline = fgets(fid);  
end

map = containers.Map(keys, vals);

fclose(fid);
