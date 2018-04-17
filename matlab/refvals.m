function [p, rho, v, t, l] = refvals(geometry, simulation)

if nargin == 0
    [geomPath, simulPath] = getPath;
else 
    [geomPath, simulPath] = getPath(geometry, simulation);
end

path = strcat(simulPath, 'referenceValues.txt');
fid = fopen(path);

tline = fgets(fid);
while ischar(tline)
    strarray = strsplit(tline, '#'); % get rid of comments
    str = strcat(strarray{1}, ';');  % append simicolon to suppress the output
    eval(str);
    
    tline = fgets(fid);    
end

fclose(fid);
