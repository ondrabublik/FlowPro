function runTurekhronBenchmark

% start flowpro
args turekhron
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
            iter = 1;
            fun = @(z,t) z/20 .* sin(2*pi*z) * sin(2*pi*t);                        
            
        case 'def'        
            % read data from FlowPro
            bodyIndex = ois.readInt();
            forceCoord = ois.readObject();
            forceRbfCoef = ois.readObject();
            t = ois.readDouble();
            dt = ois.readDouble();
            %-------------------------                                 
            
            xy = forceCoord;
            maxx = max(max(forceCoord(:,1)));
            [~, inds] = sort(xy(:,1));
            xy = xy(inds, :);                        
            sideIndexes = xy(:,1) < maxx - 1e-5;
            endIndexes = xy(:,1) > maxx - 1e-5;
            dx = (xy(3,1) - xy(1,1)) / 2;
            xy(sideIndexes, 1) = xy(sideIndexes, 1) + dx;
            xy(endIndexes, 2) = sum(xy(endIndexes, 2)) / 2;
            
            if iter == 1
                plot(forceCoord(:,1), forceCoord(:,2), '.')
                hold on
                axis equal
                plot(xy(:,1), xy(:,2), 'ro');
                
                xy0 = xy;
                x = xy0(:,1)';
                n = length(x);
                deformation = zeros(2,n);
                dyOld = zeros(1,n);
                minx = min(x);
                maxx = max(x);
                x0 = (x - minx) / (maxx - minx);
                figure
                hold on
                tisk = plot(x, 0*x, 'o');
                tisk2 = plot(x, 0*x, '.');
            end
            iter = iter + 1;
            
            boundaryCoord = xy;                                                
            deformation(2,:) = fun(x0, t);            

%             set(tisk, 'Xdata', xy0(:,1)' + deformation(1,:));
            set(tisk, 'Ydata', xy0(:,2)' + dyOld);
%             set(tisk2, 'Xdata', forceCoord(:,1));
            set(tisk2, 'Ydata', boundaryCoord(:,2));
            dyOld = deformation(2,:);
            drawnow;
            
            % write data to FlowPro
            oos.writeObject(boundaryCoord);
            oos.writeObject(deformation);
            oos.flush;
            %-------------------------
        case 'next'
%             uOld2 = uOld;
%             uOld = u;
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



    
    
    
    

