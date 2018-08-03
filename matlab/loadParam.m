function  par = loadParam(simulPath)

if nargin == 0
    [geomPath, simulPath, outputPath] = getPath;
end
    
path = strcat(simulPath, 'parameters.txt');
fid = fopen(path);

tline = fgets(fid);

keys = {};
vals = {};

while ischar(tline)
    strarray = strsplit(tline, '%'); % get rid of comments    
    str = strtrim(strarray{1});
    
    if ~isempty(str) && ~strcmp(str(end), '=')
        strarray = strsplit(str, '=');

        key = strtrim(strarray{1});
        valueStr = strtrim(strarray{2});
        value = str2num(valueStr);
        if isempty(value)
            value = valueStr;
        end

        keys = [keys, key];
        vals = [vals, value];
    end
        
    tline = fgets(fid);
end
fclose(fid);

par = cell2struct(vals, keys, 2);

if ~isfield(par, 'kappa')
    par.kappa = par.cp / par.cv;
elseif ~isfield(par, 'cp')
    par.cp = par.kappa * par.cv;
elseif ~isfield(par, 'cv')
    par.cv = par.kappa / par.cp;
end

