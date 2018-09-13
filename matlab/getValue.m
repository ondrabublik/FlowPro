function val = getValue(propName,file)

[meshPath, simulPath, outputPath] = getPath;
if(nargin == 1)
    fileName = [simulPath,'parameters.txt'];
else
    fileName = [simulPath,file];
end

% Read txt into cell A
try
    fid = fopen(fileName,'r');
    tline = fgetl(fid);
    while ischar(tline)
        str = strsplit(tline,'=');
        if(strcmp(str{1},propName))
            val = looseComment(str{2});
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
catch
    disp(' File error! ')
    val = '';
end

function str = looseComment(str)
    for i = 1:length(str)
        if(str(i) == '%')
            str = str(1:i-1);
            break
        end
    end