function runExternLinearSolverCRS
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
            ilu = ois.readObject();
            rhs = ois.readObject();
            
%             n = length(rhs);
%             A = sparse(I+1,J+1,H,n,n);
            
            % UMFPACK ==================================
%             tic
%             x = A\rhs;
%             toc
            % UMFPACK ==================================
            
            
            
            % gmres =================================================
            setup.type = 'nofill';
            setup.droptol = 0;
            setup.udiag = 0;
            tic
%             x = gmres(@(x)prod(x,I,J,H),rhs,10,1e-3,10);
            x = gmres(@(x)prod(x,I,J,H),rhs,20,1e-2,10,@(x)mfun(x,I,J,ilu,diag));
            toc
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

function y = prod(x,I,J,H)
y = zeros(size(x));
n = length(x);
for i = 1:n
    y(i) = 0; 
    for j = I(i)+1:I(i+1)
        y(i) = y(i) + H(j) * x(J(j)+1);
    end
end

function x = mfun(b,I,J,ilu,diag)
x = zeros(size(b));
y = zeros(size(b));
n = length(b);
for i = 1:n
    sum = 0;
    for j = I(i)+1:diag(i)
        sum = sum + ilu(j) * y(J(j)+1);
    end
    y(i) = b(i) - sum;
end
    
for i = n:-1:1
    sum = 0; 
    for j = diag(i) + 1:I(i+1)
        sum = sum + ilu(j) * x(J(j)+1); 
    end
    x(i) = (b(i) - sum) / ilu(diag(i)+1); 
end





