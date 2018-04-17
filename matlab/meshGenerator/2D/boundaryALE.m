function boundaryALE(action)
global P TP type obd Xobd Yobd bound;
warning off;

% graficke uzivatelske prostredi
if (nargin < 1)
    action ='initialize';
end;

if strcmp(action,'initialize')
    clc;
    
    figure(...
        'Name','Definition of ALE boundaries', ...
        'Position',[170 150 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.85 0.9], ...
        'Tag', 'fig');
    
    %====================================
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.82 0.1 0.05], ...
        'visible', 'on', ...
        'String','choose points', ...
        'Callback','boundaryALE(''choosePoints'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.9 0.75 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'neuman','stationary wall','1','2','3','4','5','6','7','8','9','10'}, ...
        'tag','typeALE');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.7 0.1 0.05], ...
        'visible', 'on', ...
        'String','set type', ...
        'Callback','boundaryALE(''setType'')');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.18 0.1 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','save boundary', ...
        'Callback','boundaryALE(''saveBoundary'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Close', 'Callback','boundaryALE(''konec'')');
    
    % inicializacni funkce
    cd ../..
    [meshPath, simulPath, outputPath] = getPath;
    P = load([meshPath,'vertices.txt']);
    TP = load([meshPath,'elements.txt'])+1;
    type = load([meshPath,'elementType.txt']);
    cd meshGenerator/2D
    
    % vykresluje sit
    maxX = max(P(:,1));
    minX = min(P(:,1));
    maxY = max(P(:,2));
    minY = min(P(:,2));
    axis([(minX-0.1) (maxX+0.1) (minY-0.1) (maxY+0.1)]);
    axis('equal');
    kresliSit;
    bound = najdiHranicniBody; % vraci indexy bodu na hranici
    % vykresluje body hranice
    hold on;
    plot(bound(:,1),bound(:,2),'.','Color','g');
    hold off;
    
    % inicializacni funkce    
    obd = [];
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'choosePoints') % vybira body hranice
    if(isempty(obd) == 0)
        delete(obd);
    end
    [Xobd,Yobd] = ginput(2);
    hold on;
    obd = plot([Xobd(1),Xobd(2),Xobd(2),Xobd(1),Xobd(1)],[Yobd(1),Yobd(1),Yobd(2),Yobd(2),Yobd(1)],'r');
    hold off;

%--------------------------------------------------------------------------     
elseif strcmp(action,'setType') % prirazuje zvoleny typ hranici
    if(Xobd(2) < Xobd(1))
        pom = Xobd(1);
        Xobd(1) = Xobd(2);
        Xobd(2) = pom;
    end
    if(Yobd(2) < Yobd(1))
        pom = Yobd(1);
        Yobd(1) = Yobd(2);
        Yobd(2) = pom;
    end
    I = 1:size(bound,1);
    IB = I(bound(:,1) > Xobd(1) & bound(:,1) < Xobd(2) & bound(:,2) > Yobd(1) & bound(:,2) < Yobd(2));
    hran = get(findobj('tag','typeALE'), 'Value');
    hold on;
    bound(IB,5) = hran;
    col = 'gbrkmcyrkmcy';
    plot(bound(IB,1),bound(IB,2),'.','Color',col(hran));
    hold off;
    
%--------------------------------------------------------------------------  
elseif strcmp(action,'saveBoundary') % vypocitava vsechny hodnoty site potrebne pro vypocet
    cla;
    hold on;
    kresliSit;
    col = 'gbrkmcyrkmcy';
    hold on;
    for i = 1:size(bound,1)
        plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'color',col(bound(i,5)),'linewidth',2);
    end
    hold off;
    
    cd ../..
    [meshPath, simulPath, outputPath] = getPath;
    bound(:,[3,4,5]) = bound(:,[3,4,5])-1;
    bound = bound(:,[5,3,4]);
    bound = bound(bound(:,1) > 0, :);
    dlmwrite(strcat(meshPath, 'boundaryTypeALE.txt'), bound, 'delimiter', ' ');
    cd meshGenerator/2D
    
    
%--------------------------------------------------------------------------     
elseif strcmp(action,'konec') % ukoncuje program
    close(clf);
end
 
function bound = najdiHranicniBody
global P TP type

    m = length(P(:,1));
    n = length(TP(:,1));
    S = zeros(m,1);
    for j = 1:n;
        if(type(j) == 3)
            S(TP(j,1)) = S(TP(j,1)) + abs(acos(((P(TP(j,2),1)-P(TP(j,1),1)).*(P(TP(j,3),1)-P(TP(j,1),1)) + (P(TP(j,2),2)-P(TP(j,1),2)).*(P(TP(j,3),2)-P(TP(j,1),2)))./(((P(TP(j,2),1)-P(TP(j,1),1)).^2 + (P(TP(j,2),2)-P(TP(j,1),2)).^2).^(1/2).*((P(TP(j,3),1)-P(TP(j,1),1)).^2 + (P(TP(j,3),2)-P(TP(j,1),2)).^2).^(1/2))));
            S(TP(j,2)) = S(TP(j,2)) + abs(acos(((P(TP(j,3),1)-P(TP(j,2),1)).*(P(TP(j,1),1)-P(TP(j,2),1)) + (P(TP(j,3),2)-P(TP(j,2),2)).*(P(TP(j,1),2)-P(TP(j,2),2)))./(((P(TP(j,3),1)-P(TP(j,2),1)).^2 + (P(TP(j,3),2)-P(TP(j,2),2)).^2).^(1/2).*((P(TP(j,1),1)-P(TP(j,2),1)).^2 + (P(TP(j,1),2)-P(TP(j,2),2)).^2).^(1/2))));
            S(TP(j,3)) = S(TP(j,3)) + abs(acos(((P(TP(j,1),1)-P(TP(j,3),1)).*(P(TP(j,2),1)-P(TP(j,3),1)) + (P(TP(j,1),2)-P(TP(j,3),2)).*(P(TP(j,2),2)-P(TP(j,3),2)))./(((P(TP(j,1),1)-P(TP(j,3),1)).^2 + (P(TP(j,1),2)-P(TP(j,3),2)).^2).^(1/2).*((P(TP(j,2),1)-P(TP(j,3),1)).^2 + (P(TP(j,2),2)-P(TP(j,3),2)).^2).^(1/2))));
        else
            v1 = [P(TP(j,2),1) - P(TP(j,1),1),P(TP(j,2),2) - P(TP(j,1),2)];
            v2 = [P(TP(j,3),1) - P(TP(j,2),1),P(TP(j,3),2) - P(TP(j,2),2)];
            v3 = [P(TP(j,4),1) - P(TP(j,3),1),P(TP(j,4),2) - P(TP(j,3),2)];
            v4 = [P(TP(j,1),1) - P(TP(j,4),1),P(TP(j,1),2) - P(TP(j,4),2)];
            vv1 = sqrt(v1(1)^2 + v1(2)^2);
            vv2 = sqrt(v2(1)^2 + v2(2)^2);
            vv3 = sqrt(v3(1)^2 + v3(2)^2);
            vv4 = sqrt(v4(1)^2 + v4(2)^2);
            S(TP(j,1)) = S(TP(j,1)) + abs(acos(-(v1(1)*v4(1) + v1(2)*v4(2))/(vv1*vv4)));
            S(TP(j,2)) = S(TP(j,2)) + abs(acos(-(v2(1)*v1(1) + v2(2)*v1(2))/(vv2*vv1)));
            S(TP(j,3)) = S(TP(j,3)) + abs(acos(-(v3(1)*v2(1) + v3(2)*v2(2))/(vv3*vv2)));
            S(TP(j,4)) = S(TP(j,4)) + abs(acos(-(v4(1)*v3(1) + v4(2)*v3(2))/(vv4*vv3)));
        end
    end

    jplusT = [2,3,1];
    jplusQ = [2,3,4,1];
    eps = 0.1;
    bound = [];
    for i = 1:n
        for j = 1:type(i)
            if(type(i) == 3)
                jp = jplusT(j);
            else
                jp = jplusQ(j);
            end
            if(S(TP(i,j)) < 2*pi-eps && S(TP(i,jp)) < 2*pi-eps)
                bound = [bound;(P(TP(i,j),[1,2]) + P(TP(i,jp),[1,2]))/2, TP(i,j), TP(i,jp), 1];
            end
        end
    end
    
    boundPom = myqs(bound);
    bound = [];
    s = 1;
    nBound = length(boundPom(:,1));
    while s <= nBound
        if(s < nBound)
            if(norm(boundPom(s,[1,2]) - boundPom(s+1,[1,2])) < 1e-5)
                s = s + 2;
            else
                bound = [bound; boundPom(s,:)];
                s = s + 1;
            end
        else
            bound = [bound; boundPom(s,:)];
            s = s + 1;
        end
    end

% _________________________________________________________________________
function kresliSit
global P TP type;
    hold on;
    for i = 1:length(TP(:,1))
        if(type(i) == 3)
            plot([P(TP(i,1),1), P(TP(i,2),1), P(TP(i,3),1), P(TP(i,1),1)], [P(TP(i,1),2), P(TP(i,2),2), P(TP(i,3),2), P(TP(i,1),2)],'g');
        else
            plot([P(TP(i,1),1), P(TP(i,2),1), P(TP(i,3),1), P(TP(i,4),1), P(TP(i,1),1)], [P(TP(i,1),2), P(TP(i,2),2), P(TP(i,3),2), P(TP(i,4),2), P(TP(i,1),2)],'g');
        end
    end
    hold off;

    
function out = myqs(data)
if(isempty(data) ~= 1)
    n = length(data(:,1));
    if(n < 2)
        out = data;
    else
        ind = floor(n/2);
        p = data(ind,:);
        L = zeros(size(data));
        iL = 1;
        R = zeros(size(data));
        iR = 1;
        if(n > 1)
            for i = 1:n
                if(i ~= ind)
                    if(je_mensi(data(i,:),p))
                        L(iL,:) = data(i,:);
                        iL = iL + 1;
                    else
                        R(iR,:) = data(i,:);
                        iR = iR + 1;
                    end
                end
            end
            L = myqs(L(1:iL-1,:));
            R = myqs(R(1:iR-1,:));
            out = [L; p; R];
        end
    end
else
    out = [];
end

function c = je_mensi(a,b)
c = 0;
if(b(1) - a(1) > 1e-5)
    c = 1;
elseif(abs(a(1) -b(1)) < 1e-5)
    if(b(2) - a(2) > 1e-5)
        c = 1;
    end
end    
    
    