function props = loadPropertiesFile(filePath)

fid = fopen(filePath);

keys = {};
vals = {};

tline = fgets(fid);
while ischar(tline)
    strarray = strsplit(tline, '%'); % get rid of comments    
    tline = strtrim(strarray{1});
    strarray = strsplit(tline, '#');
    str = strtrim(strarray{1});
    
    if ~isempty(str) && ~strcmp(str(end), '=')
        strarray = strsplit(str, '=');

        key = strtrim(strarray{1});
        valueStr = strtrim(strarray{2});
        value = str2double(valueStr);
        if isnan(value) || isempty(value)
            value = valueStr;
        end

        keys = [keys, key];
        vals = [vals, value];
    end
        
    tline = fgets(fid);
end
fclose(fid);

props = cell2struct(vals, keys, 2);
