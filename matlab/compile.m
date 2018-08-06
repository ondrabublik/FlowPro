path = pwd;
cd(getFlowProPath);
% copyfile('../FlowProUser/dist/FlowProUser.jar','lib/FlowProUser.jar')
%copyfile('../FlowProAPI/dist/FlowProAPI.jar','lib/FlowProAPI.jar')
system('ant');
cd(path)

%updatePackage
