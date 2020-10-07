function argsHere

splitPath = strsplit(pwd,'\');
for i = 1:length(splitPath)
    if(strcmp(splitPath{i},'simulations'))
        name = '';
        for j = (i+1):length(splitPath)-1
            name = [name,splitPath{j},'\'];
        end
        sim = splitPath{length(splitPath)};
    end
end

args(name,sim)
