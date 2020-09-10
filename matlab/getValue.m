function val = getValue(propName, fileName)
% getValue   Return the value of a specified parameter in parameters.txt.
%   getValue(propName) returns a value of the property propName.

if nargin < 2
    fileName = 'parameters.txt';
end

[~, simulPath, ~] = getPath;
fileFullName = [simulPath, fileName];

% Read txt into cell A
val = '';
try
    fid = fopen(fileFullName,'r');
    tline = fgetl(fid);
    while ischar(tline)
        str = strsplit(tline,'=');
        if(strcmp(str{1},propName))
            val = strtrim(looseComment(str{2}));
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