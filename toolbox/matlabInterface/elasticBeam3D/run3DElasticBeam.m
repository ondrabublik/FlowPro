function run3DElasticBeam

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
            L = 1;
            n = 15;
            y = linspace(0,L,n)';
            z = linspace(0,1,n);
            
            boundaryCoord = [y,zeros(size(y))];
            figure
            col = 'rbg';
        case 'def'
            % read data from FlowPro
            bodyIndex = ois.readInt();
            forceCoord = ois.readObject();
            forceRbfCoef = ois.readObject();
            t = ois.readDouble();
            dt = ois.readDouble();
            %-------------------------
            
            deformation = zeros(3,length(z));
            deformation(2,:) = 0.1*z.^2*sin(t + pi*bodyIndex/2);
            
            if(bodyIndex == 0)
                hold off
            else
                hold on
            end
            plot(deformation(2,:),'color',col(bodyIndex+1));
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



    
    
    
    

