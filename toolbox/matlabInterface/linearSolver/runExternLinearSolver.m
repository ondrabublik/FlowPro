function runExternLinearSolver
close all;

% start flowpro
run
disp('pausing ...')
pause(5)
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
            rhs = ois.readObject();
            n = length(rhs);
            A = sparse(I+1,J+1,H,n,n);
            
            % UMFPACK ==================================
            tic
            x = A\rhs;
            toc
            % UMFPACK ==================================
            
            
            %amg ========================================================
%             opt = amgset;                               % default option file
%             opt = amgset(opt,'coarsest',10);            % set the number of levels
%             opt = amgset(opt,'PreCond','pcg');          % set the Krylov method
%             opt = amgset(opt,'PrintOnScreen','off');    % turn off the option to print the log on the screen 
%             opt = amgset(opt,'SaveCsn','off');           % save the set of coarse-grid points
%             opt = amgset(opt,'CsnType','amg');          % choose the coarsening method
% 
%             % Initial vector
%             x = rand(length(A),1);
% 
%             % Solve a linear system
%             x = amg(A,x,rhs,opt);
%             norm(A*x-rhs)
            %amg ========================================================
            
            % gmres =================================================
%             setup.type = 'nofill';
%             setup.droptol = 0;
%             setup.udiag = 0;
%             tic
%             [L,U] = ilu(A,setup);
%             toc
%             tic
%             x = gmres(A,rhs,10,1e-3,5,L,U);
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
