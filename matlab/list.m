function list
flowProPath = getFlowProPath;

searchDir([flowProPath,'/simulations/'], '');

end

function searchDir(absPath, relPath)
    dirContent = dir(absPath);

    for i = 1 : length(dirContent)
        if dirContent(i).isdir && ~strcmp(dirContent(i).name, '.') && ~strcmp(dirContent(i).name , '..')
            name = dirContent(i).name;
            if exist([absPath, name, '/', 'mesh'], 'dir')
                fprintf(1, '%s%s: ', relPath, name);
                searchGeom([absPath, name, '/'])
            else
                searchDir([absPath, name, '/'], [relPath, name, '/']);
            end            
        end
    end
end

function searchGeom(absPath)

    simulList = dir(absPath);
    for j = 1 : length(simulList)
        if simulList(j).isdir && ~strcmp(simulList(j).name, '.') ...
                && ~strcmp(simulList(j).name , '..') && ~strcmp(simulList(j).name , 'mesh')
            fprintf(1, '%s ', simulList(j).name);
        end
    end

    fprintf(1, '\n');
end
