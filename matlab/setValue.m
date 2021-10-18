function setValue(propName, propVal)
% setValue   Change value of a parameter in parameters.txt.

if isnumeric(propVal)
    if numel(propVal) == 1
        propVal = num2str(propVal);
    else
        propVal = matrix2str(propVal);
    end
end

[meshPath, simulPath, outputPath] = getPath;
fileName = [simulPath,'parameters.txt'];

% Read txt into cell A
fid = fopen(fileName,'r');
tline = fgetl(fid);
A{1} = [];
i = 1;
entryExists = false;
while ischar(tline)
    A{i} = tline;
    str = strsplit(tline,'=');
    if strcmp(str{1},propName)
        entryExists = true;
        A{i} = [propName,' = ', propVal];
    end
    
    tline = fgetl(fid);
    i = i + 1;
end
fclose(fid);

if ~entryExists
    A{i} = [propName,' = ', propVal];
end

fid = fopen(fileName, 'w');
for i = 1:length(A)
    fprintf(fid,'%s\n', A{i});
end
fclose(fid);


function str = matrix2str(A)
str = sprintf('%.7e', A(1));
for i = 2 : numel(A)
    str = sprintf('%s, %.7e', str, A(i));
end
