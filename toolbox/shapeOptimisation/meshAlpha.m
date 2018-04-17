function meshAlpha(alpha,name)
verticesDeform = deformation(alpha);

fid = fopen(name,'w');
fprintf(fid,'%6.16f %6.16f\n',verticesDeform');
fclose(fid);


