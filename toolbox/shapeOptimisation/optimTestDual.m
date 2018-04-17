function optimTestDual
global profileNURBS optimisationPath paramFile hDer

    order = 1;
    nDomains = 1;

    % initialise
    deformation;
    currentPath = pwd;
    setValue(paramFile,'numberOfAlpha',num2str(2*profileNURBS.nc));
    setValue(paramFile,'steps',num2str(500));
    
    % optimalizacni parametry
    hDer = 1e-3;
    alpha = zeros(2*profileNURBS.nc,1);
    
    setValue(paramFile,'continueComputation','false');
    copyVerticesBack;
    
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
    save dual dI;
    
    
    
    
        