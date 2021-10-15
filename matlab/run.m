function run(runInTerminal, nDomains)
% run   Run FlowPro.

if nargin == 0
    runInTerminal = true;
end

if nargin <= 1
    %metis(1);
    path = pwd;
    cd(getFlowProPath);
    if isunix
        system('java -d64 -Xmx12g -jar FlowPro.jar master 0');
    else
%         system('java -d64 -Xmx12g -Xss100m -jar FlowPro.jar master 0');
%         system('start java -d64 -Xmx12g -Xss100m -jar FlowPro.jar master 0');
        if runInTerminal
            system('start run.bat');
        else
            system('java -d64 -Xmx12g -Xss10m -jar FlowPro.jar master 0');
        end
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

    command = ['java -d64 -Xmx8g -jar FlowPro.jar master ', num2str(nDomains)];
    if runInTerminal
        command = ['start ', command];
    end

    system(command);
        
    cd(path)
end
