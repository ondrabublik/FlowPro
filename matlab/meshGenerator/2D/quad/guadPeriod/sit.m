function sit

n1 = 30;
n2 = 20;
lx = 3;

X = zeros(n1+1,n2+1);
Y = zeros(n1+1,n2+1);

mkdir('mesh')
fid = fopen('mesh/vertices.txt','w');
for i = 1:n1+1
    for j = 1:n2+1
        X(i,j) = lx*(i-1)/n1;
        Y(i,j) = (j-1)/n2;
        fprintf(fid,'%12.8f %12.8f\n',X(i,j),Y(i,j));
    end
end
fclose(fid);

ind = zeros(n1+1,n2+1);
s = 0;
for i = 1:n1+1
    for j = 1:n2+1
        ind(i,j) = s;
        s = s + 1;
    end
end

fid = fopen('mesh/elements.txt','w');
for i = 1:n1
    for j = 1:n2
        fprintf(fid,'%d %d %d %d\n',ind(i,j),ind(i+1,j),ind(i+1,j+1),ind(i,j+1));
    end
end
fclose(fid);

ind = zeros(n1,n2);
s = 0;
for i = 1:n1
    for j = 1:n2
        ind(i,j) = s;
        s = s + 1;
    end
end

fid = fopen('mesh/neighbors.txt','w');
for i = 1:n1
    for j = 1:n2
        if(j > 1)
            Q1 = ind(i,j-1);
        else
            Q1 = -1;
        end
        
        if(j < n2)
            Q3 = ind(i,j+1);
        else
            Q3 = -1;
        end
            
        if(i > 1)
            Q4 = ind(i-1,j);
        else
            Q4 = ind(n1,j);
        end
        
        if(i < n1)
            Q2 = ind(i+1,j);
        else
            Q2 = ind(1,j);
        end
        
        fprintf(fid,'%d %d %d %d\n',Q1,Q2,Q3,Q4);
    end
end
fclose(fid);

fid = fopen('mesh/elementType.txt','w');
for i = 1:(n1*n2)
    fprintf(fid,'4\n');
end
fclose(fid);

fid = fopen('mesh/TEshift.txt','w');
for i = 1:n1
    for j = 1:n2
        Q1 = 0;
        Q2 = 0;
        Q3 = 0;
        Q4 = 0;
               
        if(i == 1)
            Q4 = 1;
        end
        
        if(i == n1)
            Q2 = 2;
        end
        
        fprintf(fid,'%d %d %d %d\n',Q1,Q2,Q3,Q4);
    end
end
fclose(fid);

fid = fopen('mesh/shift.txt','w');
fprintf(fid,'%d %d\n',lx,0);
fprintf(fid,'%d %d\n',-lx,0);
fclose(fid);









