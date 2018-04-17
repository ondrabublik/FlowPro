function GUIsitTroj_mod(action)
global geom node cnect prvek P T typ IH obd Xobd Yobd;
% graficke uzivatelske prostredi
if (nargin < 1)
    action ='initialize';
end;

if strcmp(action,'initialize')
    load 'geom';
    figure(...
        'Name','Program pro generovani nestrukturovanych siti', ...
        'Position',[170 50 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.6 0.9],'box','on', ...
        'Tag', 'fig');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.95 0.4 0.04],'visible', 'on', 'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani geometrie');
    
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
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.90 0.05 0.04], ...
        'visible', 'on', ...
        'String','obdelnik', ...
        'Callback','GUIsitTroj_mod(''obdelnik'')');
    
    % elipsa
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sx');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.67 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sx');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.70 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','sy');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.72 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','sy');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','a');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.77 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','a');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.80 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','b');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.82 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','1','tag','b');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.85 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','alfa');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.87 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','0','tag','alfa_el');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.90 0.85 0.02 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','n');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.92 0.85 0.03 0.04],'BackGroundColor','w','visible', 'on', 'String','10','tag','n_el');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.95 0.85 0.05 0.04], ...
        'visible', 'on', ...
        'String','elipsa', ...
        'Callback','GUIsitTroj_mod(''elipsa'')');
    
    % nacteni geometrie
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.80 0.12 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','Zadej jmeno souboru');
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
        'String','Nacti', ...
        'Callback','GUIsitTroj_mod(''nacti'')');
    
    % tlacitko pro prijmuti prvku
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.92 0.73 0.08 0.05], ...
        'visible', 'on', ...
        'String','Prijmi prvek', ...
        'Callback','GUIsitTroj_mod(''prijmiPrvek'')');
    
    % tlacitko pro smazani geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.84 0.73 0.08 0.05], ...
        'visible', 'on', ...
        'String','Smaz geometrii', ...
        'Callback','GUIsitTroj_mod(''smazgeometrii'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.68 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Zahusteni site');
    
    % tlacitka pro zahusteni site
    uicontrol(...
        'Style','radiobutton', ...
        'Units','normalized', ...
        'Position',[0.8 0.62 0.1 0.05], ...
        'BackGroundColor',[0.8,0.8,0.8], ...
        'visible', 'on', ...
        'HitTest', 'off', ...
        'String','Zahustit?', ...
        'Tag','zahustit');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.62 0.1 0.05], ...
        'visible', 'on', ...
        'String','Zadej funkci', ...
        'Callback','GUIsitTroj_mod(''zadejfunkci'')');
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.56 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro vytvoreni site');
    
    % tlacitko pro generovani site
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.50 0.04 0.04],'BackGroundColor',[0.8,0.8,0.8],'visible', 'on', 'String','hmax');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.79 0.50 0.04 0.04],'BackGroundColor','w','visible', 'on', 'String','0.1','tag','hmax');
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.84 0.50 0.08 0.05], ...
        'visible', 'on', ...
        'String','Generuj sit', ...
        'Callback','GUIsitTroj_mod(''generujsit'')');
    
     % tlacitko pro vyhlazeni site
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.92 0.50 0.08 0.05], ...
        'visible', 'on', ...
        'String','Vyhlad sit', ...
        'Callback','GUIsitTroj_mod(''vyhladsit'')');
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.44 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani okrajovych podminek');
    
    % tlacitka pro zadavani okrajovych podminek
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.70 0.38 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyber body', ...
        'Callback','GUIsitTroj_mod(''vyberbody'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.8 0.38 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'vstup','vystup','stena','nevazka stena'}, ...
        'tag','druhZobraz');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.38 0.1 0.05], ...
        'visible', 'on', ...
        'String','Prirad typ hranice', ...
        'Callback','GUIsitTroj_mod(''priradtyp'')');
    
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.32 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Vypocteni geometrie');
    
    % textova pole pro informace o siti
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.26 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Pocet trojuhelniku:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.26 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocTroj');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.21 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Pocet bodu:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.21 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'pocBod');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.16 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Prumerna kvalita:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.16 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'prumKval');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.11 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Minimalni kvalita:');
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.11 0.1 0.05],'BackGroundColor',[0.8,0.8,0.8], 'String', '0', 'Tag', 'minKval');
    
    % tlacitko pro celkovy tisk a vypocteni geometrie
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.26 0.1 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','Vypocti sit', ...
        'Callback','GUIsitTroj_mod(''vypoctisit'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Konec', 'Callback','GUIsitTroj_mod(''konec'')');
    
    % inicializacni funkce
    zobrazHranici;
    obd = [];
    
elseif strcmp(action,'obdelnik') % pridava obdelnik do dane geometrie
    zobrazHranici;
    x = str2double(get(findobj('tag','xld'),'string'));
    y = str2double(get(findobj('tag','yld'),'string'));
    d = str2double(get(findobj('tag','d'),'string'));
    v = str2double(get(findobj('tag','v'),'string'));
    
    prvek = [x,y; x+d,y; x+d, y+v; x, y+v];
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
    soub = get(findobj('tag','soub'),'string');
    px = str2double(get(findobj('tag','px'),'string'));
    py = str2double(get(findobj('tag','py'),'string'));
    eval(['load ',soub]);
    
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
   
    
elseif strcmp(action,'generujsit') % generujesit
%     obd = [];
%     n = length(geom);
%     
%     % generovani vektoru uzlu
%     node = [];
%     for i = 1:n
%         node = [node;geom{i}];
%     end
%     
%     % generovani vektoru propojeni
%     cnect = [];
%     for i = 1:n
%         n1 = length(geom{i}(:,1));
%         n2 = length(cnect);
%         cnect = [cnect; (1:n1-1)' + n2,(2:n1)' + n2;  n1 + n2, n2+1];
%     end
%     
%     hdata.hmax = str2double(get(findobj('tag','hmax'),'string'));
%     zahustit = get(findobj('tag','zahustit'),'Value');
%     if(zahustit == 1)
%         hdata.fun = @fun;
%     end
    
%     [P,T,junk,stat] = meshfaces(node,cnect,[],hdata,[]); % vypujceny kod
    
%     set(findobj('tag','pocTroj'),'String',num2str(stat.Triangles));
%     set(findobj('tag','pocBod'),'String',num2str(stat.Nodes));
%     set(findobj('tag','prumKval'),'String',num2str(stat.Mean_quality));
%     set(findobj('tag','minKval'),'String',num2str(stat.Min_quality));
    
    load P;
    load T;
    
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

    % IH jsou indexy bodu na hranici
    I = 1:length(P(:,1));
    eps = 0.1;
    IH = I(up(I) < 2*pi-eps);

    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    plot(P(IH,1),P(IH,2),'.','Color','b');
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

    % IH jsou indexy bodu na hranici
    I = 1:length(P(:,1));
    eps = 0.1;
    IH = I(up(I) < 2*pi-eps);
    
    % tisk
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'g');
    plot(P(IH,1),P(IH,2),'.','Color','b');
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
    IB = IH(P(IH,1) > Xobd(1) & P(IH,1) < Xobd(2) & P(IH,2) > Yobd(1) & P(IH,2) < Yobd(2));
    hran = get(findobj('tag','druhZobraz'), 'Value');
    hold on;
    switch hran
        case 1
            typ(IB) = -1;
            plot(P(IB,1),P(IB,2),'.','Color','g');
        case 2
            typ(IB) = -2;
            plot(P(IB,1),P(IB,2),'.','Color','r');
       case 3
            typ(IB) = 1;
            plot(P(IB,1),P(IB,2),'.','Color','k');
       case 4
            typ(IB) = 2;
            plot(P(IB,1),P(IB,2),'.','Color','m');
    end 
    hold off;

elseif strcmp(action,'zadejfunkci') % vypocitava vsechny hodnoty site potrebne pro vypocet
    open('fun.m');
 
elseif strcmp(action,'vypoctisit') % vypocitava vsechny hodnoty site potrebne pro vypocet
    % generovani stran
    E = [T(:,[1,2]); T(:,[2,3]); T(:,[3,1])];
    E = unique(sort(E,2),'rows');
    ne = length(E(:,1));
    typE = zeros(ne,1);
    np = length(P(:,1));
    pomP = zeros(np,1); % stupen spojeni s ostatnimi body hranice - musi byt 2
    
    I = 1:ne;
    Istn = I(typ(E(I,1)) == 1 & typ(E(I,2)) == 1);
    typE(Istn) = 1;
    for ip = Istn
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    IstnV = I(typ(E(I,1)) == 2 & typ(E(I,2)) == 2);
    typE(IstnV) = 2;
    for ip = IstnV
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    Ivst = I(typ(E(I,1)) == -1 & typ(E(I,2)) == -1);
    typE(Ivst) = -1;
    for ip = Ivst
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    Ivys = I((typ(E(I,1)) == -2 & typ(E(I,2)) == -2) | (typ(E(I,1)) == -1 & typ(E(I,2)) == -2) | (typ(E(I,1)) == -2 & typ(E(I,2)) == -1));
    typE(Ivys) = -2;
    for ip = Ivys
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    Iste = I(((typ(E(I,1)) == 1 & typ(E(I,2)) == -1) | (typ(E(I,1)) == -1 & typ(E(I,2)) == 1) | (typ(E(I,1)) == 1 & typ(E(I,2)) == -2) | (typ(E(I,1)) == -2 & typ(E(I,2)) == 1)));
    typE(Iste) = 1;
    for ip = Iste
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    IsteV = I(((typ(E(I,1)) == 2 & typ(E(I,2)) == -1) | (typ(E(I,1)) == -1 & typ(E(I,2)) == 2) | (typ(E(I,1)) == 2 & typ(E(I,2)) == -2) | (typ(E(I,1)) == -2 & typ(E(I,2)) == 2)));
    typE(IsteV) = 2;
    for ip = IsteV
        pomP(E(ip,1)) = pomP(E(ip,1)) + 1;
        pomP(E(ip,2)) = pomP(E(ip,2)) + 1;
    end
    
    I = 1:ne;
    Ikoliz = I(pomP(E(I,1)) == 3 & pomP(E(I,2)) == 3);
    typE(Ikoliz) = 0;
    
    cla;
    hold on;
    triplot(T,P(:,1),P(:,2),'b');
    for i = 1:ne
        switch typE(i)
            case 1
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'k','linewidth',2);
            case 2
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'m','linewidth',2);
            case (-1)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'g','linewidth',2);
            case (-2)
                plot([P(E(i,1),1),P(E(i,2),1)],[P(E(i,1),2),P(E(i,2),2)],'r','linewidth',2);
        end
    end
    hold off;
    
    % vypocet zakladnich geometrickych vztahu potrebnych pro vypocet
    nq = length(T(:,1));
    Q = zeros(nq,4);
    I = 1:nq;
    J = 1:3;
    Q(I,J) = T(I,J);
    
    E = [E,typE];
    P = [P,typ];
    
    vypoctiGeometriiM(P,Q,E);

    
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
    
    
function vypoctiGeometriiM(P,Q,E)
h = waitbar(0,'Please wait...');
PX = P(:,1);
PY = P(:,2);

np = length(P(:,1));
nq = length(Q(:,1));
ne = length(E(:,1));

% matice obsahujici vsechny potrebne hodnoty pro vypocet
TP = zeros(nq,3);
TT = zeros(nq,3);
TE = zeros(nq,3);

I = 1:ne;
Etyp = E(I,3);

% vytvoreni matice trojuhelniku/ctvercu
Qr = sparse(np,np);
for i = 1:nq
    Qr(Q(i,1),Q(i,2)) = i;
    Qr(Q(i,2),Q(i,3)) = i;
    Qr(Q(i,3),Q(i,1)) = i;
end

waitbar(0.25,h);

% vytvoreni matice stran
Er = sparse(np,np);
for i = 1:ne
    Er(E(i,1),E(i,2)) = i;
    Er(E(i,2),E(i,1)) = i;
end

waitbar(0.5,h);

% plneni matice Tg
for i = 1:nq
    TP(i,1) = Q(i,1)-1; % indexy bodu
    TP(i,2) = Q(i,2)-1;
    TP(i,3) = Q(i,3)-1;
    TT(i,1) = Qr(Q(i,2),Q(i,1))-1; % indexy sousednich ctvercu
    TT(i,2) = Qr(Q(i,3),Q(i,2))-1;
    TT(i,3) = Qr(Q(i,1),Q(i,3))-1;
    TE(i,1) = Er(Q(i,1),Q(i,2)); % index prvni strany
    TE(i,2) = Er(Q(i,2),Q(i,3));
    TE(i,3) = Er(Q(i,3),Q(i,1));
end

waitbar(0.75,h);

for i = 1:nq
    for j = 1:3
        switch(Etyp(TE(i,j)))
            case 1 % vazka stena
                TT(i,j) = -1;
            case -1 % vstup
                TT(i,j) = -2;
            case -2 % vystup
                TT(i,j) = -3;
            case 2 % nevazka stena
                TT(i,j) = -4;
        end
    end
end

geometrie{1} = PX;
geometrie{2} = PY;
geometrie{3} = TP;
geometrie{4} = TT;

save 'geometrie' geometrie;

waitbar(1,h);
close(h);









