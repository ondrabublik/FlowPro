function remote(arg)

    switch arg
        case 'zip'
            ziptask;
        case 'run'
            cd ..
            !java -jar Fetcher.jar multifetch 147.228.126.27 -pclist matlab/pclist.txt
%             !start java -jar Fetcher.jar multifetch 147.228.126.27 -pclist matlab/pclist.txt
            cd matlab
    end

    function ziptask
        [geometry, simulation] = loadArgs;

        geomDir = sprintf('simulations/%s/', geometry);
        simulDir = sprintf('%s%s/', geomDir, simulation);
        meshDir = sprintf('%smesh/', geomDir);

        jarName = 'FlowPro.jar';

        zipContent = {jarName, 'lib', 'gauss', 'gaussTri', 'gaussTetra', 'args.txt', ...
             'testkeystore.ks', strcat(geomDir,'*.txt'), ...
             strcat(simulDir, '*.txt'),strcat(meshDir, '*.txt')};
        
        zipPath = sprintf('../tasks/%s$%s.zip', strrep(geometry, '/', '$'), strrep(simulation, '/', '$'));
        zip(zipPath, zipContent, '..');

        fprintf(1, 'zip file was created in directory %s\n', zipPath);
    end
end