files = dir(fullfile('.', '*.txt'));

for i = 1 : length(files)
    name = files(i).name;
    arr = split(name, '.');
    shortName = arr{1};
    numstr = shortName(end-7:end);
    num = str2num(numstr);
    if num > 3000
        delete(name);
    end
end

%delete('a.txt')
