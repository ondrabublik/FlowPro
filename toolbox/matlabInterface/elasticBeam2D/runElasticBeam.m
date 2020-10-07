function runElasticBeam

% start flowpro
args testElasticBeam3
run

javaaddpath(pwd);
import java.net.Socket
import java.io.*

s = Socket('localhost', 5767);
ois = ObjectInputStream(BufferedInputStream(s.getInputStream));
testMsg = ois.readObject;
fprintf('Test message received: \"%s\"\n', testMsg);
oos = ObjectOutputStream(BufferedOutputStream(s.getOutputStream));
oos.writeObject('OK');
oos.flush();

figure
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
            L = 0.5;
            n = 20;
            h = L/n;
            E = 0.3;
            Is = 1;
            nu = 2*n;
            M = zeros(nu,nu);
            K = zeros(nu,nu);

            rho = 10;
            Ap = 0.1;
            Fref = 100;
            y = linspace(0,L,n)';
         
            u = zeros(nu-2,1);
            uOld = zeros(nu-2,1);
            uOld2 = zeros(nu-2,1);

            a = 4818;
            g = 729 * h;
            c = 1482;
            d = -321 * h;
            e = 172 * h * h;
            f = -73 * h * h;
            Me = rho*Ap*h/12600*[a, g, c, d; g, e, -d, f; c, -d, a, -g; d, f, -g, e];
            Ke = E*Is/h^3*[12 6*h -12 6*h; 6*h 4*h^2 -6*h 2*h^2; -12 -6*h 12 -6*h; 6*h 2*h^2 -6*h 4*h^2];
            Fe = h/2*[1,h/6,1,-h/6]';
            k = 1:4;
            l = 1:4;
            for i = 1:n-1
                M((i-1)*2 + k, (i-1)*2 + l) = M((i-1)*2 + k, (i-1)*2 + l) + Me(k,l);
                K((i-1)*2 + k, (i-1)*2 + l) = K((i-1)*2 + k, (i-1)*2 + l) + Ke(k,l);
            end
            I = 3:nu;
            M = M(I,I);
            K = K(I,I);
            alfa = 0;
            beta = 0;
            B = alfa*M + beta*K;

            t = 0;
            
            XL = [-0.01*ones(size(y)),y];
            XR = [0.01*ones(size(y)),y];
            
            boundaryCoord = [zeros(size(y)),y];
            
            figure
            J = 1:2:nu-2;
            subplot(1,2,1)
            tisk = plot([0;u(J)],y);
            subplot(1,2,2)
            hold on;
            tisk2 = plot(0);
            tisk3 = plot(0);
%             axis([-0.5 0.5 0 0.5]);
            
        case 'def'
            % read data from FlowPro
            bodyIndex = ois.readInt();
            forceCoord = ois.readObject();
            forceRbfCoef = ois.readObject();
            t = ois.readDouble();
            dt = ois.readDouble();
            %-------------------------
            
            fxL = Fref*(rbf(XL,forceCoord,forceRbfCoef(1,:)));
            fxR = Fref*(rbf(XR,forceCoord,forceRbfCoef(1,:)));
            fx = fxL + fxR;
%             fy = Fref*(rbf(XL,forceCoord,forceRbfCoef(2,:))-rbf(XR,forceCoord,forceRbfCoef(2,:)));
            
            F = zeros(nu,1);
            for i = 1:n-1
                F((i-1)*2 + k) = F((i-1)*2 + k) + (fx(i+1)+fx(i))/2*Fe(k);
            end
            F = F(I);
            A = M/dt^2 + B/dt + K;
            u = A\(M*(2*uOld - uOld2)/dt^2 + B*uOld/dt + F);
            deformation = zeros(2,n);
            J = 1:2:nu-2;
            deformation(1,:) = [0,u(J)'];
            
            set(tisk,'Xdata',[0;u(J)]);
            set(tisk2,'Ydata',fxL);
            set(tisk3,'Ydata',-fxR);
            drawnow;
            
            % write data to FlowPro
            oos.writeObject(boundaryCoord);
            oos.writeObject(deformation);
            oos.flush;
            %-------------------------
        case 'next'
            uOld2 = uOld;
            uOld = u;
            t = t + dt;
    end
    oos.writeObject('OK');
    oos.flush;
end

function z = rbf(X,Y,c)
    n = size(X,1);
    m = size(Y,1);
    z = zeros(n,1);
    for i = 1:n
        for j = 1:m
            z(i) = z(i) + (1-norm(X(i,:)-Y(j,:)))*c(j);
        end
    end



    
    
    
    

