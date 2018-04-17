function matlLinSol

while 1
    z = load('zamek.txt');
    if(z == -1)
        break;
    end

    if(z == 0)
        A = load('A.txt');
        b = load('b.txt');
        n = length(b);

        A = sparse(A(:,1)+1,A(:,2)+1,A(:,3),n,n);
        tic
        U = A\b;
        toc

        fid = fopen('x.txt','w');
        fprintf(fid,'%12.8f\n',U);
        fclose(fid);

        fid = fopen('zamek.txt','w');
        fprintf(fid,'%i\n',1);
        fclose(fid);
    end

    pause(0.1)
end