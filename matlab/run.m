function run(nDomains)
% run   Run FlowPro.

if nargin == 0
    %metis(1);
    path = pwd;
    cd(getFlowProPath);
    if isunix
        system('java -d64 -Xmx12g -jar FlowPro.jar master 0');
    else
        %system('java -d64 -Xmx8g -Xss100m -jar FlowPro.jar local');
        system('start run.bat');
    end
    cd(path)
else    
    if ischar(nDomains)
        nDomains = str2double(nDomains);
    end
    
%     try
%         metis(nDomains);
%     catch
        disp('Metis does not work. Running simple script for domain division.');
        simpleMetis(nDomains,'x');
%     end
    
    path = pwd;
    cd(getFlowProPath);

%     system(['java -d64 -Xmx8g -jar FlowPro.jar master ',num2str(nDomains)]);
    system(['start java -d64 -Xmx8g -jar FlowPro.jar master ', num2str(nDomains)]);

    cd(path)
end
