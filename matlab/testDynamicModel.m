function testDynamicModel

path = pwd;
cd(getFlowProPath);
system('java -d64 -Xmx8g -jar FlowPro.jar testDynamicModel 0.01');
cd(path)