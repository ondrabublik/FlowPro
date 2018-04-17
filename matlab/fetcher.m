function fetcher(command)

path = pwd;
cd(getFlowProPath);
system(sprintf('java -jar Fetcher.jar %s', command));
cd(path)
