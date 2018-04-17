function val = getValue(propName)

[meshPath, simulPath, outputPath] = getPath;
fileName = [simulPath,'parameters.txt'];

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