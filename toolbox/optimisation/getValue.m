function val = getValue(fileName,propName)

% Read txt into cell A
try
    fid = fopen(fileName,'r');
    tline = fgetl(fid);
    while ischar(tline)
        str = strsplit(tline,'=');
        if(strcmp(str{1},propName))
            val = str{2};
            break;
        end
        tline = fgetl(fid);
    end
    fclose(fid);
catch
    disp(' File error! ')
    val = '';
end