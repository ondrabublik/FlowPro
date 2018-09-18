function startFlowProMatlab
path = pwd;
addpath([path,'/matlab']);
fid = fopen('matlab/flowProPath.txt','w');
fprintf(fid,'%s\n',path);
fclose(fid);
