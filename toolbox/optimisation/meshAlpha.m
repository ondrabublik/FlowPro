function meshAlpha(alpha,str)
verticesDeform = deform(alpha);

fid = fopen(str,'w');
fprintf(fid,'%6.16f %6.16f\n',verticesDeform');
fclose(fid);