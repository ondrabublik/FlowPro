function runElastic

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
            n = 50;
            fi = linspace(0,2*pi,n+1)';
            fi = fi(1:n);
            boundaryCoord = [0.1*cos(fi)+0.5,0.1*sin(fi)+0.5];
            deformation = zeros(2,n);
            
        case 'def'
            % read data from FlowPro
            k = ois.readInt();
            forceCoord = ois.readObject();
            forceRbfCoef = ois.readObject();
            t = ois.readDouble();
            dt = ois.readDouble();
            %-------------------------
            
            fx = rbf(boundaryCoord,forceCoord,forceRbfCoef(1,:));
            fy = rbf(boundaryCoord,forceCoord,forceRbfCoef(2,:));
            
            subplot(2,1,1)
            hold off
            plot(fx)
            hold on
            plot(fy,'r')
            axis([0,length(fx),-0.01,0.01])
            
            deformation(2,:) = sin(2*pi*t)*(boundaryCoord(:,2)-0.5)/2;
            
            tisk = boundaryCoord + deformation';
            subplot(2,1,2)
            plot(tisk(:,1),tisk(:,2))
            axis equal
            axis([0.3 0.7 0.3 0.7]);
            drawnow
            
            % write data to FlowPro
            oos.writeObject(boundaryCoord);
            oos.writeObject(deformation);
            oos.flush;
            %-------------------------
        case 'next'
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



    
    
    
    

