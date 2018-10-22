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
        eval(['show ', quantity,' -i', num2str(s), ' -p1 -fvtk -n',num2str(s)]);
%         axis([-0.6 1.2 -1 1])
%         caxis([0 1.4]);
%         caxis([0 0.08]);
%         print([simulationPath,'/fig',num2str(1000000+s)],'-dpng','-r300');
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