function generator_quad(action)
global P Q typ bound obd Xobd Yobd;
% graficke uzivatelske prostredi
if (nargin < 1)
    action ='initialize';
end;

if strcmp(action,'initialize')
    load P;
    load Q;
    
    figure(...
        'Name','Program pro generovani nestrukturovanych siti', ...
        'Position',[170 50 1000 600], ...
        'Color',[0.8,0.8,0.8]);
    
    axes(...
        'Units','normalized', ...
        'Position',[0.04 0.05 0.6 0.9], ...
        'Tag', 'fig');
  
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.84 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Tlacitka pro zadavani okrajovych podminek');
    
    % tlacitka pro zadavani okrajovych podminek
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.70 0.78 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyber body', ...
        'Callback','generator_quad(''vyberbody'')');

    uicontrol(...
        'Style','popupmenu', ...
        'Units','normalized', ...
        'Position',[0.8 0.78 0.1 0.05], ...
        'BackGroundColor','w', ...
        'String',{'stena','vstup','vystup','nevazka stena'}, ...
        'tag','druhZobraz');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.78 0.1 0.05], ...
        'visible', 'on', ...
        'String','Prirad typ hranice', ...
        'Callback','generator_quad(''priradtyp'')');
    
    
    
    %====================================
    % zelena hranice
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.65 0.72 0.4 0.04],'background','blue', 'foregroundcolor', 'white','String', 'Vypocteni geometrie');
    
    % tlacitko pro celkovy tisk, vyhlazeni a vypocteni geometrie
    uicontrol('Style','text', 'Units','normalized', 'Position',[0.75 0.66 0.1 0.04],'BackGroundColor',[0.8,0.8,0.8], 'String', 'Jmeno geometrie:');
    uicontrol('Style','edit', 'Units','normalized', 'Position',[0.85 0.66 0.15 0.05],'BackGroundColor','w','visible', 'on', 'String','...','tag','geometrie');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.65 0.66 0.1 0.05], ...
        'visible', 'on', ...
        'String','Vyhlad sit', ...
        'Callback','generator_quad(''vyhladsit'')');
    
    uicontrol(...
        'Style','push', ...
        'Units','normalized', ...
        'Position',[0.9 0.6 0.1 0.05], ...
        'visible', 'on', ...
        'background','red', ...
        'String','Vypocti sit', ...
        'Callback','generator_quad(''vypoctisit'')');
    
    % Tlacitko pro konec
    uicontrol('Style','push', 'Units','normalized', 'Position',[0.9 0 0.1 0.05], 'String', 'Konec', 'Callback','generator_quad(''konec'')');
    
    % vykresluje sit
    maxX = max(P(:,1));
    minX = min(P(:,1));
    maxY = max(P(:,2));
    minY = min(P(:,2));
    axis([(minX-0.1) (maxX+0.1) (minY-0.1) (maxY-0.1)]);
    axis('equal');
    kresliSit;
    bound = najdiHranicniBody(P,Q); % vraci indexy bodu na hranici
    % vykresluje body hranice
    hold on;
    plot(bound(:,1),bound(:,2),'.','Color','b');
    hold off;
    
    % inicializacni funkce
    typ = zeros(length(P(:,1)),1);
    obd = [];
    
    
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
            bound(IB,5) = -1;
            plot(bound(IB,1),bound(IB,2),'.','Color','k');
       case 2
            bound(IB,5) = -2;
            plot(bound(IB,1),bound(IB,2),'.','Color','g');
       case 3
            bound(IB,5) = -3;
            plot(bound(IB,1),bound(IB,2),'.','Color','r');
       case 4
            bound(IB,5) = -4;
            plot(bound(IB,1),bound(IB,2),'.','Color','m');
    end
    hold off;

elseif strcmp(action,'vypoctisit') % vypocitava vsechny hodnoty site potrebne pro vypocet
    cla;
    hold on;
    kresliSit;
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
        end
    end
    hold off;
    
    % vypocet zakladnich geometrickych vztahu potrebnych pro vypocet
    ulozGeometrii(P,Q,bound);
    

elseif strcmp(action,'vyhlad sit') % vyhlazuje sit
    np = length(P(:,1));
    nq = length(Q(:,1));
    E = [Q(:,[1,2]); Q(:,[2,3]); Q(:,[3,4]); Q(:,[4,1])];
    E = unique(sort(E,2),'rows');
    ne = length(E(:,1));
    
    PPpoc = zeros(np,1);
    PPpom{np} = 0;
    for i = 1:ne
        PPpoc(E(i,1)) = PPpoc(E(i,1)) + 1;
        PPpoc(E(i,2)) = PPpoc(E(i,2)) + 1;  
        PPpom{E(i,1)} = [PPpom, E(i,2)];
        PPpom{E(i,2)} = [PPpom, E(i,1)];
    end
    PPpocmax = max(Ppoc);
    PP = zeros(np,PPpocmax);
    for i = 1:np
        for j = 1:PPpoc(i)
            PP(i,j) = PPpom{i}(j);
        end
    end
    
    % vlastni vyhlazeni
    I = 1:np;
    I = I(I ~= IH);
    for k = 1:3
        for i = I
            sx = 0;
            sy = 0;
            for j = 1:PPpoc(i)
                sx = sx + P(PP(i,j),1);
                sy = sy + P(PP(i,j),2);
            end
            sx = sx/PPpoc(i);
            sy = sy/PPpoc(i);
            P(i,1) = P(i,1)+sx;
            P(i,2) = P(i,2)+sy;
        end
    end
    kresliSit;
    
    
elseif strcmp(action,'konec') % ukoncuje program
    close(clf);
end
 
    
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
dlmwrite(strcat(cesta, 'meshIndexes.txt'), TP-1, 'delimiter', ' ');
dlmwrite(strcat(cesta, 'boundaryType.txt'), bound(:,[5,3,4]), 'delimiter', ' ');
dlmwrite(strcat(cesta, 'elementType.txt'), 4*ones(size(TP,1),1), 'delimiter', ' ');

% zapis nazev geometrie do souboru args.txt
fout = fopen('../../../args.txt', 'w');
fprintf(fout, '%s %s\n', jmenoGeometrie, 'default');
fclose(fout);

msg = sprintf('geometrie s nazvem %s byla ulozena do adresare simulations/%s', ...
               jmenoGeometrie, jmenoGeometrie);
msgbox(msg);

waitbar(1,h);
close(h);


function kresliSit
global P Q;
    hold on;
    for i = 1:length(Q(:,1))
        plot([P(Q(i,1),1), P(Q(i,2),1), P(Q(i,3),1), P(Q(i,4),1), P(Q(i,1),1)], [P(Q(i,1),2), P(Q(i,2),2), P(Q(i,3),2), P(Q(i,4),2), P(Q(i,1),2)],'g');
    end
    hold off;

function bound = najdiHranicniBody(P,Q)
m = length(P(:,1));
n = length(Q(:,1));
S = zeros(m,1);
for j = 1:n;
    v1 = [P(Q(j,2),1) - P(Q(j,1),1),P(Q(j,2),2) - P(Q(j,1),2)];
    v2 = [P(Q(j,3),1) - P(Q(j,2),1),P(Q(j,3),2) - P(Q(j,2),2)];
    v3 = [P(Q(j,4),1) - P(Q(j,3),1),P(Q(j,4),2) - P(Q(j,3),2)];
    v4 = [P(Q(j,1),1) - P(Q(j,4),1),P(Q(j,1),2) - P(Q(j,4),2)];
    vv1 = sqrt(v1(1)^2 + v1(2)^2);
    vv2 = sqrt(v2(1)^2 + v2(2)^2);
    vv3 = sqrt(v3(1)^2 + v3(2)^2);
    vv4 = sqrt(v4(1)^2 + v4(2)^2);
    S(Q(j,1)) = S(Q(j,1)) + abs(acos(-(v1(1)*v4(1) + v1(2)*v4(2))/(vv1*vv4)));
    S(Q(j,2)) = S(Q(j,2)) + abs(acos(-(v2(1)*v1(1) + v2(2)*v1(2))/(vv2*vv1)));
    S(Q(j,3)) = S(Q(j,3)) + abs(acos(-(v3(1)*v2(1) + v3(2)*v2(2))/(vv3*vv2)));
    S(Q(j,4)) = S(Q(j,4)) + abs(acos(-(v4(1)*v3(1) + v4(2)*v3(2))/(vv4*vv3)));
end

jplus = [2,3,4,1];
eps = 0.1;
bound = [];
for i = 1:n
    for j = 1:4
        jp = jplus(j);
        if(S(Q(i,j)) < 2*pi-eps && S(Q(i,jp)) < 2*pi-eps)
            bound = [bound;(P(Q(i,j),[1,2]) + P(Q(i,jp),[1,2]))/2, Q(i,j), Q(i,jp), 0];
        end
    end
end

function s = obsah(A,B,C)
v1 = B-A;
v2 = C-A;
s = abs(v1(1)*v2(2) - v1(2)*v2(1))/2;
