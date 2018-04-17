function runDoor

% start flowpro
run

javaaddpath(pwd);
import java.net.Socket
import java.io.*

s = Socket('localhost', 5767);
ois = ObjectInputStream(BufferedInputStream(s.getInputStream));
ois.readObject;
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
            y0 = L/2;
            n = 20;
         
            alfa = [0,0];
            omega = [0,0];

            y = linspace(y0,L,n);
            x = zeros(size(y));
            h = y(2)-y(1);
            
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
            
            xa = x*cos(alfa) - y*sin(alfa);
            ya = x*sin(alfa) + y*cos(alfa);
            
            
            Mx = Fref*(rbf([xa,ya],forceCoord,forceRbfCoef(1,:))).*ya*h;
            My = Fref*(rbf([xa,ya],forceCoord,forceRbfCoef(2,:))).*xa*h;
            
            M = sum(Mx+My);
             = A\(M*(2*uOld - uOld2)/dt^2 + B*uOld/dt + F);
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



    
    
    
    

