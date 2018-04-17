function getParam

path = pwd;
cd(getFlowProPath);
system('java -jar FlowPro.jar getparameters');
cd(path)
