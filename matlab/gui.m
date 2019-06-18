function gui

path = pwd;
cd(getFlowProPath);
system('java -jar FlowProManager.jar gui');
cd(path)