function optimalisationExport(str)
    path = pwd;
    cd(getFlowProPath);
    system(['java -d64 -Xmx8g -jar FlowPro.jar optimalisationExport ', str]);
    cd(path)