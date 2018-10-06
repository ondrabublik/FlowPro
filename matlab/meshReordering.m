function meshReordering

    [meshPath, simulPath, outputPath] = getPath;

    fileNeigh = [meshPath,'neighbors.txt'];
    if ~(exist(fileNeigh, 'file') == 2)
        steps = getValue('steps');
        path = pwd;
        cd(getFlowProPath);
        system('java -d64 -Xmx12g -jar FlowPro.jar local');
        cd(path);
        setValue('steps',steps);
    end
    
    TP = loadStructure(fileNeigh);
        n = size(TP,1);
        A = sparse(n,n);
        for i = 1:n
            v = str2num(TP{i}) + 1;
            for j = 1:length(v);
                if(v(j) > 0)
                    A(i,v(j)) = 1;
                end
            end
        end
        figure
        subplot(1,2,1)
        title('before reordering')
        spy(A)
        r = symrcm(A);
        subplot(1,2,2)
        title('after reordering')
        spy(A(r,r))

        fileElements = [meshPath,'elements.txt'];
        TE = loadStructure(fileElements);
        fid = fopen(fileElements,'w');
        for i = 1:size(TE,1)
            fprintf(fid,'%s\n',TE{r(i)});
        end
        fclose(fid);
        
        fileType = [meshPath,'elementType.txt'];
        typ = load(fileType);
        fid = fopen(fileType,'w');
        for i = 1:size(typ,1)
            fprintf(fid,'%d\n',typ(r(i)));
        end
        fclose(fid);
        
        delete([meshPath,'neighbors.txt']);
        delete([meshPath,'wallDistance.txt']);
end


function s = loadStructure(str)
    fid = fopen(str,'r');
    n = 0;
    while(feof(fid) == 0)
        fgetl(fid);
        n = n + 1;
    end
    fclose(fid);
    
    s = cell(n);
    fid = fopen(str,'r');
    i = 1;
    while(feof(fid) == 0)
        line = fgetl(fid);
        s{i} = line;
        i = i + 1;
    end
    fclose(fid);   
end


