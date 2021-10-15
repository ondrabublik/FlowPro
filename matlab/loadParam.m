function  par = loadParam(simulationPath)
% loadParam   Return parameters of the simulation in a structure.

if nargin == 0
    [~, simulationPath, ~] = getPath;
end
    
filePath = strcat(simulationPath, 'parameters.txt');

par = loadPropertiesFile(filePath);

% if ~isfield(par, 'kappa')
%     par.kappa = par.cp / par.cv;
% elseif ~isfield(par, 'cp')
%     par.cp = par.kappa * par.cv;
% elseif ~isfield(par, 'cv')
%     par.cv = par.kappa / par.cp;
% end
