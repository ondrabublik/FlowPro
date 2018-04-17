function optimizeExport
global profileNURBS meshPath optimisationPath paramFile hDer
    
    close all
    clc

    saveInitVertices;
    
    % initialise
    deformation;
    currentPath = pwd;
    setValue(paramFile,'numberOfAlpha',num2str(2*profileNURBS.nc));
    
    % optimalizacni parametry
    hDer = 1e-3;
    alpha = zeros(2*profileNURBS.nc,1);
     
    % run computation
    cd(getFlowProPath);   
    system('java -d64 -Xmx8g -jar FlowPro.jar local');
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
    
    save dI dI;
    
%     % make new mesh and save to computational geometry folder
%     meshAlpha(alpha,[meshPath,'vertices.txt']);
% 
%     show mach pressure density velocity -fvtk


