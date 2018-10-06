function runExternLinearSolverCRS2
close all;
clc

% start flowpro
run
disp('pausing ...')
pause(2)
disp('runing ...')

javaaddpath(pwd);
import java.net.Socket
import java.io.*

s = Socket('localhost', 5767);
ois = ObjectInputStream(BufferedInputStream(s.getInputStream));
ois.readObject;
oos = ObjectOutputStream(BufferedOutputStream(s.getOutputStream));
oos.writeObject('OK');
oos.flush();

while 1
    try
        tag = ois.readObject;
    catch
        disp('Connection closed!');
        close all;
        break;
    end
    switch tag
        case 'init'

        case 'solve'
            % read data from FlowPro
            I = ois.readObject();
            J = ois.readObject();
            H = ois.readObject();
            diag = ois.readObject();
            ilu0 = ois.readObject();
            rhs = ois.readObject();
            n = length(rhs);
            
            I = expand(I,n);
            I = cast(I,'double');
            J = cast(J,'double');
            A = sparse(I+1,J+1,H,n,n);
            ILU0 = sparse(I+1,J+1,ilu0,n,n);

            U = triu(ILU0);
            L = tril(ILU0);
            L = L - spdiags(spdiags(L,0),0,n,n) + spdiags(ones(n,1),0,n,n);
            
            % UMFPACK ==================================
%             tic
%             x = A\rhs;
%             toc
            % UMFPACK ==================================
            
            % gmres =================================================
%             tic
%             setup.type = 'nofill';
%             [L,U] = ilu(A,setup);

            [x,flag] = gmres(A,rhs,30,1e-2,5,L,U);
            disp(num2str(flag))
            
%             toc
            % gmres =================================================
            
            % write data to FlowPro
            oos.writeUnshared(x);
            oos.flush;
            %-------------------------
    end
    oos.writeObject('OK');
    oos.flush;
    oos.reset();
end

function J = expand(I,n)
J = zeros(n,1);
s = 1;
for i = 1:length(I)-1
    for j = I(i):I(i+1)-1
        J(s) = i-1;
        s = s + 1;
    end
end



