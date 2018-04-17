function optimize(continueComput)
global meshPath optimisationPath
    
    if(nargin == 0)
        continueComput = 0;
    end

    [geometry, simulation] = loadArgs;
    geometryPath = strcat('../../simulations/', geometry, '/');
    meshPath = strcat(geometryPath, '/mesh/');
    simulationPath = strcat(geometryPath, simulation, '/');
    optimisationPath = strcat(geometryPath, simulation, '/optimisation/');
    paramFile = [simulationPath,'parameters.txt'];

    % vlastni optimalizace
    h = 1e-3;
    k = 0.5;
    alpha = 0;
    alphaEvo = [alpha];
    IEvo = [];
    derI = [];
    optStep = 20;
    optTol = 1e-4;
    
    % initMesh
    figure('color','w')
    subplot(2,1,1)
    plot(alphaEvo,'color','k','linewidth',2);
    xlabel('iteration','fontsize',14);
    ylabel('alpha','fontsize',14);
    set(gca,'fontsize',14);
    box on;
    grid on;
    title('\alpha','fontsize',14)
    
    subplot(2,1,2)
    plot(IEvo,'color','k','linewidth',2);
    xlabel('iteration','fontsize',14);
    ylabel('functional','fontsize',14);
    set(gca,'fontsize',14);
    box on;
    grid on;
    title('f','fontsize',14)
    
    currentPath = pwd;
    for opt = 1:optStep
        if(opt == 1 && continueComput == 0)
            setValue(paramFile,'continueComputation','false');
            copyVerticesBack;
            meshAlpha(0,[optimisationPath,'vertices.txt']);
        else
            setValue(paramFile,'continueComputation','true');
        end
        % run computation
        cd(getFlowProPath);
        nDomains = 2;
        if(nDomains > 1)
            simpleMetis(nDomains,'x');
            system(['java -d64 -Xmx8g -jar FlowPro.jar master ', num2str(nDomains)])
        else
            system('java -d64 -Xmx8g -jar FlowPro.jar local');
        end
        cd(currentPath)
        
        % compute PXY(alpha+h)
        meshAlpha(alpha+h,[optimisationPath,'/vertices',num2str(1),'.txt']);
        
        % export residuals
        cd(getFlowProPath);
        optimalisationExport;
        cd(currentPath);
        
        % compute gradient using adjoint method
        dI = getFunctionalGradient(length(alpha),h);
        dAlpha = limit(k*dI);
        alpha = alpha + dAlpha;
        
        % make new mesh and save to computational geometry folder
        meshAlpha(alpha,[meshPath,'vertices.txt']);
        
        % convergence process
        derI = [derI;dI]
        alphaEvo = [alphaEvo;alpha*180/pi]
        IEvo = [IEvo;load([optimisationPath,'I0.txt'])]
        
%         if opt > 1 && abs(IEvo(end)-IEvo(end-1)) < optTol
%             disp('Optimalisation converged! :-)');
%             break;
%         end
        
        % tisk
        subplot(2,1,1)
        plot(alphaEvo,'color','k','linewidth',2);
        subplot(2,1,2)
        plot(IEvo,'color','k','linewidth',2);
        drawnow;
        
        show mach pressure density velocity -fvtk
    end
    
    % show results
    cd(getFlowProPath);
    show mach
    cd(currentPath);
end

function out = limit(in)
    max = 0.2;
    if(abs(in) > max)
        out = max*sign(in);
    else
        out = in;
    end
end


