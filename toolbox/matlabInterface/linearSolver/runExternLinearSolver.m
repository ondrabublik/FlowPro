function runExternLinearSolver
close all;

% start flowpro
run
pause(5)

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
            rhs = ois.readObject();
            n = length(rhs);
            A = sparse(I+1,J+1,H,n,n);
            
            tic
            x = A\rhs;
            toc

% %             setup.type = 'crout';
% %             setup.milu = 'row';
% %             setup.droptol = 0.1;
% %             [L,U] = ilu(A,setup);
%             tic
%             x = bicgstab(A,rhs,1e-4,200);
%             toc
            
            % write data to FlowPro
            oos.writeUnshared(x);
            oos.flush;
            %-------------------------
    end
    oos.writeObject('OK');
    oos.flush;
    oos.reset();
end
