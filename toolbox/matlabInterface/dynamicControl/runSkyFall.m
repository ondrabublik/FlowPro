function runSkyFall

% start flowpro
run
pause(3)

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
%         close all;
        break;
    end
    switch tag
        case 'init'
            nBody = ois.readInt();
            disp([num2str(nBody),' bodies'])
            x = zeros(nBody,1);
            y = zeros(nBody,1);
            alfa = zeros(nBody,1);
            
            u = zeros(nBody,1);
            v = zeros(nBody,1);
            om = zeros(nBody,1);
            uo = zeros(nBody,1);
            vo = zeros(nBody,1);
            omo = zeros(nBody,1);
            Fxo = zeros(nBody,1);
            Fyo = zeros(nBody,1);
            Mo = zeros(nBody,1);
            
            m = 10;
            I = 100;
            
            xt = 0;
            yt = 0;
            alfat = 0;
            fxt = 0;
            fyt = 0;
            Mt = 0;
            tt = 0;
            
            figure
            subplot(1,2,1)
            hold on
            tiskFx = plot(tt,fxt);
            tiskFy = plot(tt,fyt,'r');
            tiskM = plot(tt,Mt,'g');
            subplot(1,2,2)
            hold on
            tiskx = plot(tt,xt);
            tisky = plot(tt,yt,'r');
            tiskalfa = plot(tt,alfat,'g');
        case 'def'
            % read data from FlowPro
            t = ois.readDouble();
            dt = ois.readDouble();
            Fx = ois.readObject();
            Fy = ois.readObject();
            M = ois.readObject();
            %-------------------------
            Fthrust = -1;
            u(2) = uo(2) + dt*(3/2*Fx(2)-1/2*Fxo(2) + Fthrust*sin(alfa(2)))/m;
            v(2) = vo(2) + dt*(3/2*Fy(2)-1/2*Fyo(2) + Fthrust*cos(alfa(2)))/m;
            om(2) = omo(2) + dt*(3/2*M(2)-1/2*Mo(2))/I;
            x(2) = x(2) + dt*(3/2*u(2)-1/2*uo(2));
            y(2) = y(2) + dt*(3/2*v(2)-1/2*vo(2));
            alfa(2) = alfa(2) + dt*(3/2*om(2)-1/2*omo(2));
            
            % domain
            x(1) = x(2);
            y(1) = y(2);
            alfa(1) = alfa(2);
            
            % kormidlo
            h = 2;
            x(3) = x(2) + h*sin(alfa(2));
            y(3) = y(2) + h*cos(alfa(2));
            alfa(3) = alfa(2);
            
            uo = u;
            vo = v;
            omo = om;
            Fxo = Fx;
            Fyo = Fy;
            Mo = M;
            
            % write data to FlowPro
            oos.writeObject(x);
            oos.writeObject(y);
            oos.writeObject(alfa);
            oos.flush;
            %-------------------------
            
            tt = [tt,t];
            xt = [xt,x(2)];
            yt = [yt,y(2)];
            alfat = [alfat,alfa(2)];
            fxt = [fxt,Fxo(2)];
            fyt = [fyt,Fyo(2)];
            Mt = [Mt,M(2)];
            
            set(tiskx,'Xdata',tt,'Ydata',xt);
            set(tisky,'Xdata',tt,'Ydata',yt);
            set(tiskalfa,'Xdata',tt,'Ydata',alfat);
            set(tiskFx,'Xdata',tt,'Ydata',fxt);
            set(tiskFy,'Xdata',tt,'Ydata',fyt);
            set(tiskM,'Xdata',tt,'Ydata',Mt);
            drawnow;
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



    
    
    
    

