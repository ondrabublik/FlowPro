function makeFigures(quantity, step, endIter)

if(nargin < 3)
    endIter = 1e5;
else
    if isstring(endIter)
        endIter = str2double(endIter);
    end
end

if isstring(step)
    step = str2double(step);
end
s = 0;
while 1
    try
        eval(['show ', quantity,' -i', num2str(s), ' -fvtk -n',num2str(s)]);
        close all;
        display(['step: ',num2str(s)]);
        s = s + step;
    catch
        break;
    end
    if(s > endIter)
        break;
    end
end