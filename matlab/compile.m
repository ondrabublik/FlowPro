path = pwd;
cd(getFlowProPath);
system('gradle build');
cd(path)

updatePackage
