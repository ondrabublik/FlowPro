function list
flowProPath = getFlowProPath;

geomList = dir([flowProPath,'/simulations']);
for i = 1 : length(geomList)
    if geomList(i).isdir && ~strcmp(geomList(i).name, '.') && ~strcmp(geomList(i).name , '..')
        fprintf(1, '%s: ', geomList(i).name);
        
        simulList = dir(strcat(flowProPath,'/simulations/', geomList(i).name));
        for j = 1 : length(simulList)
            if simulList(j).isdir && ~strcmp(simulList(j).name, '.') ...
                    && ~strcmp(simulList(j).name , '..') && ~strcmp(simulList(j).name , 'mesh')
                fprintf(1, '%s ', simulList(j).name);
            end
        end
        
        fprintf(1, '\n');
    end
end
