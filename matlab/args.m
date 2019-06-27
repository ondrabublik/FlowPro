%args   Change or print geometry and simulation name.
%   args(geometry, simulation) saves geometry and simulation name
%   (parametrs for the JAVA application) into a text file.
%
%   args(geometry) sets simulation as default.
%
%   args prints current geometry and simulation name.

function args(geometry, simulation)

flowProPath = getFlowProPath;

    switch nargin
        case 0
            [geometry, simulation] = loadArgs;
            fprintf(1, 'Geometry name:   %s\n', geometry);
            fprintf(1, 'Simulation name: %s\n', simulation);
            return;
            
        case 1
            disp('Setting simulation to default.');
            simulation = 'default';                                   
    end
    
    path = sprintf('%s/simulations/%s/', flowProPath, geometry);
    if ~exist(path, 'dir')
        error('Geometry %s does not exist.', geometry);
    end
    path = strcat(path, simulation);
    if ~exist(path, 'dir')
        fprintf(1, 'Simulation %s/%s does not exist, ', geometry, simulation);
        fprintf(1, 'the simulation will be created after running the command "param".\n');
    end

    fout = fopen([flowProPath, '/args.txt'], 'w');
    fprintf(fout, '%s %s\n', geometry, simulation);
    fclose(fout);
end
