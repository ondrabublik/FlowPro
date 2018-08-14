function show(varargin)
    % show - export data for post-processing
    % input parameters: 
    % show <variableName> -i<iteration> -f<format>    
    %
    % formats: vtk (Paraview), txt (plain text)
    % 
    % variable names: depend on the model used
    %
    % example:
    %
    %   show mach pressure -i100 -ftxt
    
    str = 'java -d64 -Xmx8g -jar FlowPro.jar postprocessing ';
    for i = 1:nargin
        str = [str,' ', varargin{i}];
    end
    
    path = pwd;
    cd(getFlowProPath);
    system(str);
    cd(path)
end
