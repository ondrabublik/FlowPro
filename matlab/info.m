function info
% Gives overview of matlab scripts and functions available in FlowPro.

flowProPath = getFlowProPath();
matlabPath = strcat(flowProPath, '/matlab/');
scripts = dir(strcat(matlabPath, '/*.m'));

for i = 1 : length(scripts)
    name = strrep(scripts(i).name, '.m', '');
    
    lines = splitlines(help(name));
    [~, linehelp] = strtok(lines{1});
    
    fprintf('%s%s\n', pad(name, 22, ' '), strtrim(linehelp));
end
