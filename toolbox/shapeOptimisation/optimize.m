function optimize(continueComput)
global profileNURBS meshPath optimisationPath paramFile hDer
    
    close all
    clc
    if(nargin == 0)
        continueComput = 0;
    end

    % initialise
    deformation;
    currentPath = pwd;
    setValue(paramFile,'numberOfAlpha',num2str(2*profileNURBS.nc));
    setValue(paramFile,'steps',num2str(200));
    
    % optimalizacni parametry
    hDer = 1e-4;
    k = 10;
    alpha = zeros(2*profileNURBS.nc,1);
    alphaEvo = alpha;
    IEvo = [];
    optStep = 9;
    optTol = 1e-5;
    nDomains = 1;
    order = 1;
    dAlphaOld = zeros(size(alpha));
    
    % omezeni
    area0 = getArea(alpha); 
    centre0 = getCentre(alpha);
    minThick = 0.005;
    
    figure('color','w')
    subplot(2,2,1)
    xlabel('iteration','fontsize',14);
    ylabel('alpha','fontsize',14);
    set(gca,'fontsize',14);
    box on;
    grid on;
    title('\alpha','fontsize',14)
    
    subplot(2,2,3)
    plot(IEvo,'color','k','linewidth',2);
    xlabel('iteration','fontsize',14);
    ylabel('functional','fontsize',14);
    set(gca,'fontsize',14);
    box on;
    grid on;
    title('f','fontsize',14)
    
    col = 'rgbkmyrgbkmyrgbkmyrgbkmyrgbkmy';
    for opt = 1:optStep
        if(opt == 1 && continueComput == 0)
            setValue(paramFile,'continueComputation','false');
            copyVerticesBack;
            %meshAlpha(alpha,[optimisationPath,'vertices.txt']);
        else
            setValue(paramFile,'continueComputation','true');
        end
        % run computation
        setValue(paramFile,'order',num2str(order));
        cd(getFlowProPath);
        if(nDomains > 1)
            simpleMetis(nDomains,'x');
            system(['java -d64 -Xmx8g -jar FlowPro.jar master ', num2str(nDomains)])
        else
            system('java -d64 -Xmx8g -jar FlowPro.jar local');
        end
        cd(currentPath);
        
        % compute PXY(alpha+h)
        for i = 1:length(alpha)
            betha = alpha;
            betha(i) = betha(i) + hDer;
            meshAlpha(betha,[optimisationPath,'vertices',num2str(i),'.txt']);
            disp(['Deformation parameter: ',num2str(i)]);
        end
        
        % export residuals
        setValue(paramFile,'continueComputation','true');
        cd(getFlowProPath);
        optimalisationExport;
        cd(currentPath);
        
        % compute gradient using adjoint method
        dI = getFunctionalGradient(length(alpha));
        lpf = 1;
        dAlpha = lpf*limit(k*dI) + (1-lpf)*dAlphaOld;
        dAlphaOld = dAlpha;
        alpha = alpha + dAlpha;
        
        % area projection
        dArea = zeros(profileNURBS.nc,1);
        dCentre = zeros(profileNURBS.nc,2);
        dThick = zeros(profileNURBS.nc,1);
        for j = 1:80
            area = getArea(alpha);
            centre = getCentre(alpha);
            thick = getMinimalThickness(alpha);
            for i = 1:length(alpha)
                betha = alpha;
                betha(i) = betha(i) + hDer;
                dArea(i) = (getArea(betha)-area)/hDer;
                dCentre(i,:) = (getCentre(betha)-centre)/hDer;
                dThick(i) = (getMinimalThickness(betha)-thick)/hDer;
            end
            lambdaArea = area-area0;
            lambdaCentre = centre-centre0;
            lambdaThick = 10*min(thick - minThick, 0);
            alpha = alpha - lambdaArea*dArea - lambdaCentre(1)*dCentre(:,1) - lambdaCentre(2)*dCentre(:,2) - lambdaThick*dThick;
            disp(['iter: ',num2str(j),' area = ', num2str(lambdaArea),', ',num2str(area),', centre = ', num2str(norm(lambdaCentre)), ', thick = ', num2str(lambdaThick),', ',num2str(thick)]);
        end
        
        % make new mesh and save to computational geometry folder
        meshAlpha(alpha,[meshPath,'vertices.txt']);
        
        % convergence process
        alphaEvo = [alphaEvo, alpha];
        IEvo = [IEvo; load([optimisationPath,'I0.txt'])];
        
        save alpha alphaEvo;
        save fun IEvo;
        
%         if opt > 1 && abs(IEvo(end)-IEvo(end-1)) < optTol
%             disp('Optimalisation converged! :-)');
%             break;
%         end
        
        % tisk
        subplot(2,2,1)
        hold off;
        plot(alphaEvo(1,:),'color',col(1));
        hold on;
        for j = 2:length(alpha)
            plot(alphaEvo(j,:),'color',col(j));
        end
        subplot(2,2,3)
        plot(IEvo,'color','k','linewidth',2);
        drawnow;
        
        subplot(1,2,2)
        hold off
        data = nurbsClosed(profileNURBS);
        plot(data(:,1),data(:,2),'color','k');
        hold on
        profilePom = profileNURBS;
        profilePom.xc = profilePom.xc + alpha(1:profilePom.nc);
        profilePom.yc = profilePom.yc + alpha(profilePom.nc + [1:profilePom.nc]);
        data = nurbsClosed(profilePom);
        plot(data(:,1),data(:,2),'color','r');
        drawnow;
        
        % copy to eryx
        path = 'D:\Soubory\Ostatni\MatlabSSH\ssh2_v2_m1_r6\ssh2_v2_m1_r6';
        print([path,'/figs/fig1'],'-dpng','-r300');
        pathLoc = pwd;
        cd(path);
        copyToEryx;
        cd(pathLoc);
        
        show mach pressure density velocity -fvtk
    end
end

function out = limit(in)
    maxim = 0.1;
    out = in;
    for i = 1:length(in)
        if(abs(in(i)) > maxim)
            out(i) = maxim*sign(in(i));
        end
    end
end


