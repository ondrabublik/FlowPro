function optimTestFD
global profileNURBS optimisationPath meshPath paramFile hDer

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
    
    % export residuals
    setValue(paramFile,'continueComputation','true');
    cd(getFlowProPath);
    optimalisationExport;
    cd(currentPath);
    I0 = load([optimisationPath,'I0.txt']);

    dI = zeros(size(alpha));
    for i = 1:length(alpha)
        setValue(paramFile,'continueComputation','false');
        copyVerticesBack;
        betha = alpha;
        betha(i) = betha(i) + hDer;
        meshAlpha(betha,[meshPath,'vertices.txt']);
        disp(['Deformation parameter: ',num2str(i)]);
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
        
        % export residuals
        setValue(paramFile,'continueComputation','true');
        cd(getFlowProPath);
        optimalisationExport;
        cd(currentPath);
        
        Ia = load([optimisationPath,'I0.txt']);
        dI(i) = (Ia-I0)/hDer;
    end
    save fd dI;
    
    
    
    
        