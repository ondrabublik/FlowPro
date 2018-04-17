function [t,bodies] = showBody(varargin)
global bodies nBody t
action = '';
if(nargin == 0)
    action = 'init';
else
    action = varargin{1};
end
switch action
    case 'init'
        [meshPath, simulPath, outputPath] = getPath;
        [geometry, simulation] = loadArgs;

        body = load(strcat(simulPath, 'bodiesDynamic.txt'));
        [iter, n] = size(body);
        nBody = (n-1)/6;

        [pRef, rhoRef, uRef, tRef, lRef] = refvals(geometry, simulation);
        % structure
        t = tRef*body(:,1);
        bodies{nBody} = [];
        for i = 1:nBody
            bodies{i}.x = lRef*body(:,6*(i-1)+2);
            bodies{i}.y = lRef*body(:,6*(i-1)+3);
            bodies{i}.alfa = body(:,6*(i-1)+4);
            bodies{i}.Fx = pRef*lRef*body(:,6*(i-1)+5);
            bodies{i}.Fy = pRef*lRef*body(:,6*(i-1)+6);
            bodies{i}.momentum = pRef*lRef^2*body(:,6*(i-1)+7);
        end
        
        figure('color', 'w')
        axes('position', [0.1,0.1,0.55,0.85]);
        box on;
        axis equal;
        a = 0.025;
        b = 0.1;
        c = 0.05;
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'x',...
            'Position', [0.9-a-2*b,0.9-c,0.1,0.1],'Callback', 'showBody(1,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'y',...
            'Position', [0.9-a-b,0.9-c,0.1,0.1],'Callback', 'showBody(1,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'alfa',...
            'Position', [0.9-a,0.9-c,0.1,0.1],'Callback', 'showBody(1,3)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'ux',...
            'Position', [0.9-a-2*b,0.8-c,0.1,0.1],'Callback', 'showBody(2,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'uy',...
            'Position', [0.9-a-b,0.8-c,0.1,0.1],'Callback', 'showBody(2,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'omega',...
            'Position', [0.9-a,0.8-c,0.1,0.1],'Callback', 'showBody(2,3)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'Fx',...
            'Position', [0.9-a-2*b,0.7-c,0.1,0.1],'Callback', 'showBody(3,1)');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'Fy',...
            'Position', [0.9-a-b,0.7-c,0.1,0.1],'Callback', 'showBody(3,2)'); 
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'M',...
            'Position', [0.9-a,0.7-c,0.1,0.1],'Callback', 'showBody(3,3)');
        
        uicontrol('Style', 'edit', 'units','normalized', 'string', '0',...
            'Position', [0.9-a-2*b,0.6-c,0.1,0.1], 'tag', 'xMin');
        uicontrol('Style', 'edit', 'units','normalized', 'string', '1',...
            'Position', [0.9-a-b,0.6-c,0.1,0.1], 'tag', 'xMax');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'xAxis',...
            'Position', [0.9-a,0.6-c,0.1,0.1], ...
            'Callback', 'showBody(''xAxis'')');
        
        uicontrol('Style', 'edit', 'units','normalized', 'string', '0',...
            'Position', [0.9-a-2*b,0.5-c,0.1,0.1], 'tag', 'yMin');
        uicontrol('Style', 'edit', 'units','normalized', 'string', '1',...
            'Position', [0.9-a-b,0.5-c,0.1,0.1], 'tag', 'yMax');
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'yAxis',...
            'Position', [0.9-a,0.5-c,0.1,0.1], ...
            'Callback', 'showBody(''yAxis'')');
        
        uicontrol('Style', 'pushbutton', 'units','normalized', 'String', 'close',...
            'Position', [0.9-a,0,0.1,0.1],'Callback', 'close(gcf)');
        
    case 'xAxis'    
        xMin = str2double(get(findobj('tag','xMin'),'string'));
        xMax = str2double(get(findobj('tag','xMax'),'string'));
        xlim([xMin xMax]);
        
    case 'yAxis'    
        yMin = str2double(get(findobj('tag','yMin'),'string'));
        yMax = str2double(get(findobj('tag','yMax'),'string'));
        ylim([yMin yMax]);
        
    otherwise
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
                            plot(t,bodies{i}.x,'color',col(i));
                            hold on
                        end
                        ylabel('x [m]','fontsize',14);
                    case 2
                        hold off
                        for i = 1:nBody
                            plot(t,bodies{i}.y,'color',col(i));
                            hold on
                        end
                        ylabel('y [m]','fontsize',14);
                    case 3
                        hold off
                        for i = 1:nBody
                            plot(t,bodies{i}.alfa,'color',col(i));
                            hold on
                        end
                        ylabel('angle [rad]','fontsize',14);
                end

            case 2
                % plot velocity
                I = 1:length(t)-1;
                legenda = '';
                switch j
                    case 1
                        hold off
                        for i = 1:nBody
                            u = (bodies{i}.x(I+1)-bodies{i}.x(I))./(t(I+1)-t(I));
                            plot(t(I),u,'color',col(i));
                            hold on
                        end
                        ylabel('u [m/s]','fontsize',14);

                    case 2
                        hold off
                        for i = 1:nBody
                            v = (bodies{i}.y(I+1)-bodies{i}.y(I))./(t(I+1)-t(I));
                            plot(t(I),v,'color',col(i));
                            hold on
                        end;
                        ylabel('v [m/s]','fontsize',14);

                    case 3
                        hold off
                        for i = 1:nBody
                            omega = (bodies{i}.alfa(I+1)-bodies{i}.alfa(I))./(t(I+1)-t(I));
                            plot(t(I),omega,'color',col(i));
                            hold on
                        end
                        ylabel('omega [rad/s]','fontsize',14);
                end

            case 3
                switch j
                    case 1
                        % plot forces
                        hold off
                        for i = 1:nBody
                            plot(t,bodies{i}.Fx,'color',col(i));
                            hold on
                        end
                        ylabel('Fx [N]','fontsize',14);

                    case 2
                        hold off
                        for i = 1:nBody
                            plot(t,bodies{i}.Fy,'color',col(i));
                            hold on
                        end
                        ylabel('Fy [N]','fontsize',14);

                    case 3
                        hold off
                        for i = 1:nBody
                            plot(t,bodies{i}.momentum,'color',col(i));
                            hold on
                        end
                        ylabel('momentum [Nm]','fontsize',14);
                end
        end
        xlabel('t [s]','fontsize',14);
        legenda = '';
        for i = 1:nBody
            legenda = [legenda;'Body',num2str(i)];
        end
        legend(legenda);
        box on;
        grid on;
end
    