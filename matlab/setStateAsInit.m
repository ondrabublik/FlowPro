function setStateAsInit

[~, simulPath, ~] = getPath;
copyfile([simulPath,'W.txt'],[simulPath,'initW.txt']);