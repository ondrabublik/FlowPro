function updatePackage

path = pwd;
cd(getFlowProPath);
system('java -jar FlowProManager.jar update');
cd(path)