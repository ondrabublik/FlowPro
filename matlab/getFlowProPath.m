function path = getFlowProPath
    fid = fopen('flowProPath.txt','r');
    path = fgetl(fid);
    fclose(fid);