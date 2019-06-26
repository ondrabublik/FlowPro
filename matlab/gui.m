function gui
%gui Run FlowPro GUI.

path = pwd;
cd(getFlowProPath);
system('java -jar FlowProManager.jar gui');
cd(path)