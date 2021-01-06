function [t,bodies] = showBody3(varargin)
global bodies nBody t

if(nargin == 0)
    action = 'init';
else
    action = varargin{1};
end

switch action
    case 'init'
        [~, simulPath, ~] = getPath;
        body = load(strcat(simulPath, 'bodiesDynamic.txt'));
        [~, n] = size(body);
        nBody = (n-1)/6;

        % structure
        t = body(:,1);
        bodies{nBody} = [];
        for i = 1:nBody
            bodies{i}.x = body(:,6*(i-1)+2);
            bodies{i}.y = body(:,6*(i-1)+3);
            bodies{i}.z = body(:,6*(i-1)+4);
            bodies{i}.Fx = body(:,6*(i-1)+5);
            bodies{i}.Fy = body(:,6*(i-1)+6);
            bodies{i}.Fz = body(:,6*(i-1)+7);
        end
        
        figure('color', 'w')
        axes('position', [0.1,0.1,0.55,0.85]);
        box on;
        axis equal;
        a = 0.025;
        b = 0.1;
        c = 0.05;
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'x',...
            'Position', [0.9-a-2*b,0.9-c,0.1,0.1],'Callback', 'showBody3(1,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'y',...
            'Position', [0.9-a-b,0.9-c,0.1,0.1],'Callback', 'showBody3(1,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'z',...
            'Position', [0.9-a,0.9-c,0.1,0.1],'Callback', 'showBody3(1,3)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'ux',...
            'Position', [0.9-a-2*b,0.8-c,0.1,0.1],'Callback', 'showBody3(2,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'uy',...
            'Position', [0.9-a-b,0.8-c,0.1,0.1],'Callback', 'showBody3(2,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'uz',...
            'Position', [0.9-a,0.8-c,0.1,0.1],'Callback', 'showBody3(2,3)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'Fx',...
            'Position', [0.9-a-2*b,0.7-c,0.1,0.1],'Callback', 'showBody3(3,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'Fy',...
            'Position', [0.9-a-b,0.7-c,0.1,0.1],'Callback', 'showBody3(3,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'Fz',...
            'Position', [0.9-a,0.7-c,0.1,0.1],'Callback', 'showBody3(3,3)');
        
        uicontrol('Style', 'edit', 'units','normalized', 'string', '0',...
            'Position', [0.9-a-2*b,0.6-c,0.1,0.1], 'tag', 'xMin');
        uicontrol('Style', 'edit', 'units','normalized', 'string', '1',...
            'Position', [0.9-a-b,0.6-c,0.1,0.1], 'tag', 'xMax');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'xAxis',...
            'Position', [0.9-a,0.6-c,0.1,0.1], ...
            'Callback', 'showBody3(''xAxis'')');
        
        uicontrol('Style', 'edit', 'units','normalized', 'string', '0',...
            'Position', [0.9-a-2*b,0.5-c,0.1,0.1], 'tag', 'yMin');
        uicontrol('Style', 'edit', 'units','normalized', 'string', '1',...
            'Position', [0.9-a-b,0.5-c,0.1,0.1], 'tag', 'yMax');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'yAxis',...
            'Position', [0.9-a,0.5-c,0.1,0.1], ...
            'Callback', 'showBody3(''yAxis'')');
        
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'close',...
            'Position', [0.9-a,0,0.1,0.1],'Callback', 'close(gcf)');
        
        uicontrol('Style', 'radiobutton', 'units','normalized', 'String', 'reference values',...
            'Position', [0.9-a-2*b,0.4-c,0.3,0.1],'tag','refVals','value', 1);
        
        uicontrol('Style', 'radiobutton', 'units','normalized', 'String', 'FFT',...
            'Position', [0.9-a-2*b,0.3-c,0.3,0.1],'tag','fft','value', 0);
        
    case 'xAxis'    
        xMin = str2double(get(findobj('tag','xMin'),'string'));
        xMax = str2double(get(findobj('tag','xMax'),'string'));
        xlim([xMin xMax]);
        
    case 'yAxis'    
        yMin = str2double(get(findobj('tag','yMin'),'string'));
        yMax = str2double(get(findobj('tag','yMax'),'string'));
        ylim([yMin yMax]);
        
    otherwise
        pRef = 1;
        tRef = 1;
        lRef = 1;
        if(get(findobj('tag','refVals'),'value'))
            [geometry, simulation] = loadArgs;
            [pRef, ~, ~, tRef, lRef] = refvals(geometry, simulation);
        end
        
        fftVal = get(findobj('tag','fft'),'value');
        
        i = varargin{1};
        j = varargin{2};
        col = 'brgkmycbrgkmycbrgkmyc';
        switch i
            case 1
                % plot displacements
                switch j
                    case 1
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,lRef*bodies{i}.x,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('x [m]','fontsize',14);
                    case 2
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,lRef*bodies{i}.y,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('y [m]','fontsize',14);
                    case 3
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,lRef*bodies{i}.z,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('z [m]','fontsize',14);
                end

            case 2
                % plot velocity
                I = 1:length(t)-1;
                legenda = '';
                switch j
                    case 1
                        hold off
                        for i = 1:nBody
                            u = (bodies{i}.x(I+1)-bodies{i}.x(I))./(t(I+1)-t(I))*lRef/tRef;
                            plot(tRef*t(I),u,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('u [m/s]','fontsize',14);

                    case 2
                        hold off
                        for i = 1:nBody
                            v = (bodies{i}.y(I+1)-bodies{i}.y(I))./(t(I+1)-t(I))*lRef/tRef;
                            plot(tRef*t(I),v,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('v [m/s]','fontsize',14);

                    case 3
                        hold off
                        for i = 1:nBody
                            w = (bodies{i}.z(I+1)-bodies{i}.z(I))./(t(I+1)-t(I))*lRef/tRef;
                            plot(tRef*t(I),w,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('w [m/s]','fontsize',14);
                end

            case 3
                fRef = pRef*lRef*lRef;
                switch j
                    case 1
                        % plot forces
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,bodies{i}.Fx*fRef,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('Fx [N]','fontsize',14);

                    case 2
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,bodies{i}.Fy*fRef,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('Fy [N]','fontsize',14);

                    case 3
                        hold off
                        for i = 1:nBody
                            plot(tRef*t,bodies{i}.Fz*fRef,'color',col(i));
                            hold on
                        end
                        xlabel('t [s]','fontsize',14);
                        ylabel('Fz [N]','fontsize',14);
                end
        end
        legenda = '';
        for i = 1:nBody
            legenda = [legenda;'Body',num2str(i)];
        end
        legend(legenda);
        box on;
        grid on;
end

function [f,P1] = fourier(t,X)
    I = 1:length(t)-1;
    dt = t(I+1)-t(I);
    dtAvg = sum(dt)/length(dt);
    Fs = 1/dtAvg;
    Y = fft(X);
    L = length(Y);
    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L;
    
    
    
    
    
    
    
    
    
    
    