[geomPath, simulPath, outputPath] = getPath;
    
path = strcat(simulPath, 'parameters.txt');
fid = fopen(path);

tline = fgets(fid);

while ischar(tline)
    strarray = strsplit(tline, '%'); % get rid of comments
    str = strtrim(strarray{1});
    
    if ~isempty(str) && ~strcmp(str(end), '=')
        str = strcat(str, ';');  % append simicolon to suppress the output    
        try
            eval(str);
        end
    end
        
    tline = fgets(fid);
end

fclose(fid);
