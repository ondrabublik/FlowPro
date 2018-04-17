function runRope
global n h k kt m b G

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
            L = 0.4;
            n = 20;
            h = L/(n-1);
            x = linspace(0,L,n)';
            y = zeros(n,1);
            X = [x,y];
            X0 = X;
            
            u = zeros(4*n,1);
            for i = 1:n
                u(4*(i-1) + 1) = x(i);
                u(4*(i-1) + 2) = y(i);
            end
            un = u;

            m = 0.1/n;
            k = 1000*n/10;
            kt = 50*n/10;
            b = 0;
            G = [0,-0.005];

            Ix = 1:4:4*n;
            Iy = 2:4:4*n;

            t = 0;

            figure
            subplot(2,1,1)
            tisk = plot(u(Ix),u(Iy),'marker','o');
            axis equal
            axis([-0.1 0.7 -0.2 0.2]);
            
            dR = sparse(4*n,4*n);
            E = speye(4*n,4*n);
            hd = 1e-6;
            Pref = 500;
            
        case 'def'
            %-------------------------
            % read data from FlowPro
            bodyIndex = ois.readInt();
            forceCoord = ois.readObject();
            forceRbfCoef = ois.readObject();
            t = ois.readDouble();
            dt = ois.readDouble();
            %-------------------------
            
            dX1 = zeros(size(X));
            dX2 = zeros(size(X));
            J = 1:2;
            for i = 1:n
                if(i == n)
                    r = u(4*(i-1)+J) - u(4*(i-2)+J);
                else
                    r = u(4*i+J) - u(4*(i-1)+J);
                end
                r = r/norm(r);
                dX1(i,:) = 0.005*[-r(2),r(1)];
                dX2(i,:) = -0.005*[-r(2),r(1)];
            end
            
            pL = Pref*rbf(X+dX1,forceCoord,forceRbfCoef(1,:));
            pR = Pref*rbf(X+dX2,forceCoord,forceRbfCoef(1,:));
            dp = pR-pL;
            
            fx = zeros(n,1);
            fy = zeros(n,1);
            for i = 1:n
                n1 = [0,0];
                n2 = [0,0];
                if(i < n)
                    nor = u(4*i+J) - u(4*(i-1)+J);
                    n1 = [-nor(2),nor(1)]/norm(nor);
                end
                if(i > 1)
                    nor = u(4*(i-1)+J) - u(4*(i-2)+J);
                    n2 = [-nor(2),nor(1)]/norm(nor);
                end
                
                fx(i) = fx(i) + h/2*dp(i)*n1(1);
                fy(i) = fy(i) + h/2*dp(i)*n1(2);
                
                fx(i) = fx(i) + h/2*dp(i)*n2(1);
                fy(i) = fy(i) + h/2*dp(i)*n2(2);
            end
            
            r = fun(u,fx,fy);
            for i = 1:4*n
                u(i) = u(i) + hd;
                dR(:,i) = (fun(u,fx,fy) - r)/hd;
                u(i) = u(i) - hd;
            end

            un = u + (E/dt - dR)\r;

            set(tisk,'Xdata',u(Ix));
            set(tisk,'Ydata',u(Iy));
            subplot(2,1,2)
            hold off
            plot(pL)
            hold on
            plot(pR,'r')
            drawnow;
            
            Xold = X;
            X = [un(Ix),un(Iy)];
            
            disp(['FSI error = ',num2str(sum(sum(abs(X-Xold))))]);
            
            %-------------------------
            % write data to FlowPro
            oos.writeObject(X0);
            oos.writeObject((X-X0)');
            oos.flush;
            %-------------------------
        case 'next'
            u = un;
            t = t + dt;
            title(t);
            disp(['t = ',num2str(t)]);
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

    
function R = fun(u,fx,fy)
global n h k kt m b G

   R = zeros(4*n,1);
   I = 1:2;
   for i = 2:n
       rm = u(4*(i-2)+I) - u(4*(i-1)+I);
       rn = norm(rm);
       rm = rm/rn*(rn-h);
       F = k*rm;
       if(i < n)
           r = (u(4*(i-2)+I) + u(4*i+I))/2 - u(4*(i-1)+I);
           rp = u(4*i+I) - u(4*(i-1)+I);
           rn = norm(rp);
           rp = rp/rn*(rn-h);
           F = F + kt*r + k*rp;
       end
       R(4*(i-1) + 1) = u(4*(i-1) + 3);
       R(4*(i-1) + 2) = u(4*(i-1) + 4);
       R(4*(i-1) + 3) = (F(1) + G(1) - b*u(4*(i-1) + 3) + fx(i))/m;
       R(4*(i-1) + 4) = (F(2) + G(2) - b*u(4*(i-1) + 4) + fy(i))/m;
   end
   R(1) = 0;
   R(2) = 0;
    
    
    

