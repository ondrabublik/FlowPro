function testDynamicModel(dt)

if nargin == 0
    dt = '0.05';
end

path = pwd;
cd(getFlowProPath);
system(['java -d64 -Xmx8g -jar FlowPro.jar testDynamicModel ', dt]);
cd(path)