function makeFigures(quantity)

[~,simPath,~] = getPath;

list = dir([simPath,'/animation']);



for i = 1:length(list)
    s = getNumber(list(i).name);
    if s > -1
        eval(['show ', quantity, ' -i', num2str(s), ' -fvtk -n', num2str(s)]);
        close all;
        display(['step: ',num2str(s)]);
    end
end

function s = getNumber(name)
    s = -1;
    if(name(1) == 'W' && name(2) ~= 'e')
        s = str2num(name(2:end-4));
    end
