% graficke uzivatelske rozhrani

function generator(action)
global geom node cnect prvek P T typ bound obd Xobd Yobd;

if (nargin < 1)
    action ='initialize';
end

if strcmp(action,'initialize')
    if exist('geom.mat','file')
        load 'geom';
    else
        geom = [];
    end
    figure(...
        'Name','Unstructured mesh generator', ...
        'Position',[120 50 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.6 0.9],'box','on', ...
        'Tag', 'fig');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.95 0.4 0.04],'visible', 'on', 'background','blue', 'foregroundcolor', 'white','String', 'Geometry definition');
    
    % tlacitka pro tvorbu geometrie
    % obdelnik
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','xld');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.67 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','xld');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.71 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','yld');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.73 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','yld');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.77 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','d');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.79 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','d');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.83 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','v');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.85 0.90 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','v');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.89 0.90 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','alpha');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.92 0.90 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','alfa_sq');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.90 0.05 0.04], ...
        'visible', 'on', ...
        'String','rectangle', ...
        'Callback','generator(''obdelnik'')');
    
    % elipsa
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sx');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.67 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sx');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.70 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sy');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.72 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sy');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','a');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.77 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','a');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.80 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','b');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.82 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','b');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.85 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','alpha');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.87 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','alfa_el');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.90 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','n');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.92 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','10','tag','n_el');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.85 0.05 0.04], ...
        'visible', 'on', ...
        'String','ellipse', ...
        'Callback','generator(''elipsa'')');
    
    % nacteni geometrie
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.80 0.12 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','File name:');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.77 0.80 0.06 0.04],'BackGroundColor','w','visible', 'on', 'String','...','tag','soub');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.83 0.80 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','px');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.85 0.80 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','px');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.89 0.80 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','py');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.91 0.80 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','py');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.80 0.05 0.04], ...
        'visible', 'on', ...
        'String','load', ...
        'Callback','generator(''nacti'')');
    
    % tlacitko pro prijmuti prvku
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.73 0.1 0.05], ...
        'visible', 'on', ...
        'String','accept', ...
        'Callback','generator(''prijmiPrvek'')');
    
    % tlacitko pro smazani geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.78 0.73 0.12 0.05], ...
        'visible', 'on', ...
        'String','clear', ...
        'Callback','generator(''smazgeometrii'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.68 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Refining mesh');
    
    % tlacitka pro zahusteni site
    uicontrol(...
        'Style','radiobutton', ...
        'Units','normalized', ...
        'Position',[0.8 0.62 0.1 0.05], ...
        'BackGroundColor',[0.8,0.8,0.8], ...
        'visible', 'on', ...
        'HitTest', 'off', ...
        'String','Refine?', ...
        'Tag','zahustit');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.62 0.1 0.05], ...
        'visible', 'on', ...
        'String','refine function', ...
        'Callback','generator(''zadejfunkci'')');
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.56 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Mesh generation');
    
    uicontrol(...
        'Style','radiobutton', ...
        'Units','normalized', ...
        'Position',[0.67 0.50 0.07 0.04], ...
        'BackGroundColor',[0.8,0.8,0.8], ...
        'visible', 'on', ...
        'HitTest', 'off', ...
        'String','Load?', ...
        'Tag','loadMesh');
    
    % tlacitko pro generovani site
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.50 0.04 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','hmax');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.79 0.50 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0.1','tag','hmax');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.84 0.50 0.08 0.05], ...
        'visible', 'on', ...
        'String','generate', ...
        'Callback','generator(''generator'')');
    
     % tlacitko pro vyhlazeni site
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.92 0.50 0.08 0.05], ...
        'visible', 'on', ...
        'String','smooth', ...
        'Callback','generator(''vyhladsit'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.44 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Boundary conditions');
    
    % tlacitka pro zadavani okrajovych podminek
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.65 0.38 0.1 0.05], ...
        'visible', 'on', ...
        'String','choose points', ...
        'Callback','generator(''vyberbody'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.75 0.38 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'inlet','outlet','wall','inviscid wall','userdef'}, ...
        'tag','druhZobraz');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.85 0.38 0.15 0.05], ...
        'visible', 'on', ...
        'String','set type', ...
        'Callback','generator(''priradtyp'')');
    
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.32 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Mesh calculation');
    
    % textova pole pro informace o siti
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.26 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'elements:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.26 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocTroj');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.21 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'points:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.21 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocBod');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.16 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'average quality:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.16 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'prumKval');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.11 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'minimum quality:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.11 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'minKval');
    
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.87 0.21 0.13 0.1],'BackGroundColor',[0.8,0.8,0.8], 'String', 'geometry name:');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.87 0.23 0.13 0.05],'BackGroundColor','w','visible', 'on', 'String','...','tag','geometrie');
    
    % tlacitko pro celkovy tisk a vypocteni geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.87 0.17 0.13 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','calculate mesh', ...
        'Callback','generator(''vypoctisit'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Close', 'Callback','generator(''konec'')');
    
    % inicializacni funkce
    zobrazHranici;
    obd = [];
    
elseif strcmp(action,'obdelnik') % pridava obdelnik do dane geometrie
    zobrazHranici;
    x = str2double(get(findobj('tag','xld'),'string'));
    y = str2double(get(findobj('tag','yld'),'string'));
    d = str2double(get(findobj('tag','d'),'string'));
    v = str2double(get(findobj('tag','v'),'string'));
    alfa = str2double(get(findobj('tag','alfa_sq'),'string'));
    alfa = alfa/180*pi; 
    
    Xp = [x; x+d; x+d; x];
    Yp = [y; y; y+v; y+v];
    X = cos(alfa)*Xp + sin(alfa)*Yp;
    Y = -sin(alfa)*Xp + cos(alfa)*Yp;
    prvek = [X,Y];
    tiskni(prvek);
 
    
elseif strcmp(action,'elipsa') % pridava elipsu do dane geometrie
    zobrazHranici;
    sx = str2double(get(findobj('tag','sx'),'string'));
    sy = str2double(get(findobj('tag','sy'),'string'));
    a = str2double(get(findobj('tag','a'),'string'));
    b = str2double(get(findobj('tag','b'),'string'));
    alfa = str2double(get(findobj('tag','alfa_el'),'string'));
    n = str2double(get(findobj('tag','n_el'),'string'));
    
    I = 1:n;
    Xp = a*cos(2*pi*I/n)';
    Yp = b*sin(2*pi*I/n)';
    X = cos(alfa)*Xp + sin(alfa)*Yp + sx';
    Y = -sin(alfa)*Xp + cos(alfa)*Yp + sy';
    prvek = [X,Y];
    tiskni(prvek);
    
    
elseif strcmp(action,'nacti') % nacita libovolnou polygonialni geometrii
    zobrazHranici;
    jmenoFce = get(findobj('tag','soub'),'string');
    px = str2double(get(findobj('tag','px'),'string'));
    py = str2double(get(findobj('tag','py'),'string'));
%     fce = str2func(jmenoFce);
%     data = fce();
    eval(['data=', jmenoFce, ';']);
    
    prvek = [data(:,1)+px, data(:,2)+py];
    tiskni(prvek);
    
    
elseif strcmp(action,'prijmiPrvek') % prijima prvek geometrie
    if(isempty(prvek) ~= 1) 
        if(isempty(geom) == 1)    
            geom{1} = prvek;
        else
            geom{length(geom)+1} = prvek;
        end
    end
   prvek = [];
   zobrazHranici;
    
   
elseif strcmp(action,'smazgeometrii') % maze geometrii
    geom = [];
    obd = [];
    cla;
   
    
elseif strcmp(action,'generator') % generujesit
    obd = [];
    n = length(geom);
    
    % generovani vektoru uzlu
    node = [];
    for i = 1:n
        node = [node;geom{i}];
    end
    
    % generovani vektoru propojeni
    cnect = [];
    for i = 1:n
        n1 = length(geom{i}(:,1));
        n2 = length(cnect);
        cnect = [cnect; (1:n1-1)' + n2,(2:n1)' + n2;  n1 + n2, n2+1];
    end
    
    hdata.hmax = str2double(get(findobj('tag','hmax'),'string'));
    zahustit = get(findobj('tag','zahustit'),'Value');
    if(zahustit == 1)
        hdata.fun = @fun;
    end
    
    if(get(findobj('tag','loadMesh'),'Value') == 0)
        [P,T,junk,stat] = meshfaces(node,cnect,[],hdata,[]); % vypujceny kod
    else
        load mesh;
        P = mesh{1};
        T = mesh{2};
        stat.Triangles = size(T,1);
        stat.Nodes = size(P,1);
        stat.Mean_quality = -1;
        stat.Min_quality = -1;
    end
    
    set(findobj('tag','pocTroj'),'String',num2str(stat.Triangles));
    set(findobj('tag','pocBod'),'String',num2str(stat.Nodes));
    set(findobj('tag','prumKval'),'String',num2str(stat.Mean_quality));
    set(findobj('tag','minKval'),'String',num2str(stat.Min_quality));
    
%     load mesh;
%     P = mesh{1};
%     T = mesh{2};
    
    typ = zeros(length(P(:,1)),1);
    %P je matice, ktera ma v prvnim radku souradnice bodu X, ve druhem sloupci je souradnice Y, ve
    %tretim sloupci je okajova podminka bodu

    % urceni vnitrnich bodu
    nt = length(T(:,1));
    np = length(P(:,1));
    uhel = zeros(nt,3);
    I = 1:nt;
    uhel(I,1) = abs(acos(((P(T(I,2),1)-P(T(I,1),1)).*(P(T(I,3),1)-P(T(I,1),1)) + (P(T(I,2),2)-P(T(I,1),2)).*(P(T(I,3),2)-P(T(I,1),2)))./(((P(T(I,2),1)-P(T(I,1),1)).^2 + (P(T(I,2),2)-P(T(I,1),2)).^2).^(1/2).*((P(T(I,3),1)-P(T(I,1),1)).^2 + (P(T(I,3),2)-P(T(I,1),2)).^2).^(1/2))));
    uhel(I,2) = abs(acos(((P(T(I,3),1)-P(T(I,2),1)).*(P(T(I,1),1)-P(T(I,2),1)) + (P(T(I,3),2)-P(T(I,2),2)).*(P(T(I,1),2)-P(T(I,2),2)))./(((P(T(I,3),1)-P(T(I,2),1)).^2 + (P(T(I,3),2)-P(T(I,2),2)).^2).^(1/2).*((P(T(I,1),1)-P(T(I,2),1)).^2 + (P(T(I,1),2)-P(T(I,2),2)).^2).^(1/2))));
    uhel(I,3) = abs(acos(((P(T(I,1),1)-P(T(I,3),1)).*(P(T(I,2),1)-P(T(I,3),1)) + (P(T(I,1),2)-P(T(I,3),2)).*(P(T(I,2),2)-P(T(I,3),2)))./(((P(T(I,1),1)-P(T(I,3),1)).^2 + (P(T(I,1),2)-P(T(I,3),2)).^2).^(1/2).*((P(T(I,2),1)-P(T(I,3),1)).^2 + (P(T(I,2),2)-P(T(I,3),2)).^2).^(1/2))));

    % uhel prirazeny jednotlivym bodum
    up = zeros(np,1);
    for i = 1:nt
        up(T(i,1)) = up(T(i,1)) + uhel(i,1);
        up(T(i,2)) = up(T(i,2)) + uhel(i,2);
        up(T(i,3)) = up(T(i,3)) + uhel(i,3);
    end

    jplus = [2,3,1];
    eps = 0.1;
    bound = [];
    for i = 1:nt
        for j = 1:3
            jp = jplus(j);
            if(up(T(i,j)) < 2*pi-eps && up(T(i,jp)) < 2*pi-eps)
                bound = [bound;(P(T(i,j),[1,2]) + P(T(i,jp),[1,2]))/2, T(i,j), T(i,jp), 0];
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
    
    % tisk
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    plot(bound(:,1),bound(:,2),'.','Color','b');
    xmax = max(P(:,1));
    xmin = min(P(:,1));
    ymax = max(P(:,2));
    ymin = min(P(:,2));
    axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
    axis('equal');
    hold off;

 
elseif strcmp(action,'vyhladsit') % vyhlazuje sit
    [P,T] = smoothmesh(P,T,10);
    
    % urceni vnitrnich bodu
    nt = length(T(:,1));
    np = length(P(:,1));
    uhel = zeros(nt,3);
    I = 1:nt;
    uhel(I,1) = abs(acos(((P(T(I,2),1)-P(T(I,1),1)).*(P(T(I,3),1)-P(T(I,1),1)) + (P(T(I,2),2)-P(T(I,1),2)).*(P(T(I,3),2)-P(T(I,1),2)))./(((P(T(I,2),1)-P(T(I,1),1)).^2 + (P(T(I,2),2)-P(T(I,1),2)).^2).^(1/2).*((P(T(I,3),1)-P(T(I,1),1)).^2 + (P(T(I,3),2)-P(T(I,1),2)).^2).^(1/2))));
    uhel(I,2) = abs(acos(((P(T(I,3),1)-P(T(I,2),1)).*(P(T(I,1),1)-P(T(I,2),1)) + (P(T(I,3),2)-P(T(I,2),2)).*(P(T(I,1),2)-P(T(I,2),2)))./(((P(T(I,3),1)-P(T(I,2),1)).^2 + (P(T(I,3),2)-P(T(I,2),2)).^2).^(1/2).*((P(T(I,1),1)-P(T(I,2),1)).^2 + (P(T(I,1),2)-P(T(I,2),2)).^2).^(1/2))));
    uhel(I,3) = abs(acos(((P(T(I,1),1)-P(T(I,3),1)).*(P(T(I,2),1)-P(T(I,3),1)) + (P(T(I,1),2)-P(T(I,3),2)).*(P(T(I,2),2)-P(T(I,3),2)))./(((P(T(I,1),1)-P(T(I,3),1)).^2 + (P(T(I,1),2)-P(T(I,3),2)).^2).^(1/2).*((P(T(I,2),1)-P(T(I,3),1)).^2 + (P(T(I,2),2)-P(T(I,3),2)).^2).^(1/2))));

    % uhel prirazeny jednotlivym bodum
    up = zeros(np,1);
    for i = 1:nt
        up(T(i,1)) = up(T(i,1)) + uhel(i,1);
        up(T(i,2)) = up(T(i,2)) + uhel(i,2);
        up(T(i,3)) = up(T(i,3)) + uhel(i,3);
    end

    jplus = [2,3,1];
    eps = 0.1;
    bound = [];
    for i = 1:nt
        for j = 1:3
            jp = jplus(j);
            if(up(TP(i,j)) < 2*pi-eps && up(TP(i,jp)) < 2*pi-eps)
                bound = [bound;(P(TP(i,j),[1,2]) + P(TP(i,jp),[1,2]))/2, TP(i,j), TP(i,jp), 0];
            end
        end
    end
    
    % tisk
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    plot(bound(:,1),bound(:,2),'.','Color','b');
    xmax = max(P(:,1));
    xmin = min(P(:,1));
    ymax = max(P(:,2));
    ymin = min(P(:,2));
    axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
    axis('equal');
    hold off;
    
elseif strcmp(action,'vyberbody') % vybira body hranice
    if(isempty(obd) == 0)
        delete(obd);
    end
    [Xobd,Yobd] = ginput(2);
    hold on;
    obd = plot([Xobd(1),Xobd(2),Xobd(2),Xobd(1),Xobd(1)],[Yobd(1),Yobd(1),Yobd(2),Yobd(2),Yobd(1)],'r');
    hold off;

    
elseif strcmp(action,'priradtyp') % prirazuje zvoleny typ hranici
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
    hran = get(findobj('tag','druhZobraz'), 'Value');
    hold on;
    switch hran
        case 1
            bound(IB,5) = -2;
            plot(bound(IB,1),bound(IB,2),'.','Color','g');
       case 2
            bound(IB,5) = -3;
            plot(bound(IB,1),bound(IB,2),'.','Color','r');
       case 3
            bound(IB,5) = -1;
            plot(bound(IB,1),bound(IB,2),'.','Color','k');
       case 4
            bound(IB,5) = -4;
            plot(bound(IB,1),bound(IB,2),'.','Color','m');
       case 5
            bound(IB,5) = -5;
            plot(bound(IB,1),bound(IB,2),'.','Color','y');
    end 
    hold off;

elseif strcmp(action,'zadejfunkci') % vypocitava vsechny hodnoty site potrebne pro vypocet
    open('fun.m');
 
elseif strcmp(action,'vypoctisit') % vypocitava vsechny hodnoty site potrebne pro vypocet
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'b');
    for i = 1:size(bound,1)
        switch bound(i,5)
            case -2
                plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'g','linewidth',2);
            case -3
                plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'r','linewidth',2);
            case -1
                plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'k','linewidth',2);
            case -4
                plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'m','linewidth',2);
            case -5
                plot([P(bound(i,3),1),P(bound(i,4),1)],[P(bound(i,3),2),P(bound(i,4),2)],'y','linewidth',2);
        end
    end
    hold off;
    
    % vypocet zakladnich geometrickych vztahu potrebnych pro vypocet
    ulozGeometrii(P,T,bound);

    
elseif strcmp(action,'konec') % ukoncuje program
    save 'geom' geom;
    close(clf);
end


function zobrazHranici
    global geom;
    cla;
    if(isempty(geom) == 0)
        n = length(geom);
        xmax = -1000;
        xmin = 1000;
        ymax = -1000;
        ymin = 1000;
        hold on;
        for i = 1:n
            plot([geom{i}(:,1);geom{i}(1,1)],[geom{i}(:,2);geom{i}(1,2)],'k','linewidth',2);
            pxmax = max(geom{i}(:,1));
            pxmin = min(geom{i}(:,1));
            pymax = max(geom{i}(:,2));
            pymin = min(geom{i}(:,2));
            if(pxmax > xmax)
                xmax = pxmax;
            end
            if(pxmin < xmin)
                xmin = pxmin;
            end
            if(pymax > ymax)
                ymax = pymax;
            end
            if(pymin < ymin)
                ymin = pymin;
            end
        end
        axis([xmin-0.3 xmax+0.3 ymin-0.3 ymax+0.3]);
        axis('equal');
        hold off;
    else
        cla;
    end
 
    
function tiskni(prvek)
    prvek = [prvek; prvek(1,1) prvek(1,2)];
    hold on;
    plot(prvek(:,1),prvek(:,2),'r');
    hold off;
    
    
function ulozGeometrii(PXY,TP,bound)

jmenoGeometrie = get(findobj('tag','geometrie'), 'string');
if strcmp(jmenoGeometrie, '...')
    msgbox('nejdrive zadejte jmeno geometrie');
    return;
elseif jmenoGeometrie(1) == '.'
    msgbox('neplatne jmeno geometrie');
    return;
end

cesta = strcat('../../../../simulations/', jmenoGeometrie, '/mesh/');

if exist(cesta, 'dir')
    rozhodnuti = questdlg('Geometrie jiz existuje, chcete ji zmenit?', ...
	'Geometrie jiz existuje', 'ano', 'ne', 'ne');
    % rozhodnuti je prazdny kdyz uzivatel stiskne krizek
    if strcmp(rozhodnuti, 'ne') || isempty(rozhodnuti)
        return;
    end
else
    mkdir(cesta)
end

h = waitbar(0, 'probiha vypocet site...');
bound(:,[3,4]) = bound(:,[3,4])-1;
dlmwrite(strcat(cesta, 'vertices.txt'), PXY, 'delimiter', ' ', 'precision', 16);
dlmwrite(strcat(cesta, 'elements.txt'), TP-1, 'delimiter', ' ', 'precision', '%i');
dlmwrite(strcat(cesta, 'boundaryType.txt'), bound(:,[5,3,4]), 'delimiter', ' ', 'precision', '%i');
dlmwrite(strcat(cesta, 'elementType.txt'), 3*ones(size(TP,1),1), 'delimiter', ' ', 'precision', '%i');

% zapis nazev geometrie do souboru args.txt
fout = fopen('../../../../args.txt', 'w');
fprintf(fout, '%s %s\n', jmenoGeometrie, 'default');
fclose(fout);

msg = sprintf('geometry %s was saved into the directory simulations/%s', ...
               jmenoGeometrie, jmenoGeometrie);
msgbox(msg);

waitbar(1,h);
close(h);



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


