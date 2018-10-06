function generator(a, b, n)

nodes = linspace(a, b, n)';
elements = [(0:n-2)' (1:n-1)'];
neighbours = [[-2; (0:n-3)'] [(1:n-2)'; -3]];
type = 2 * ones(n-1, 1);

jmenoGeometrie = input('Zadejte nazev geometrie: ', 's');

cesta = strcat('../../../simulations/', jmenoGeometrie, '/mesh/');

if exist(cesta, 'dir')
    rozhodnuti = questdlg('Geometrie jiz existuje, chcete ji zmenit?', ...
	'Geometrie jiz existuje', 'ano', 'ne', 'ne');
    % rozhodnuti je prazdny kdyz uzivatel stiskne krizek
    if strcmp(rozhodnuti, 'ne') || isempty(rozhodnuti)
        return;
    end
else
    mkdir(cesta)
end

dlmwrite(strcat(cesta, 'vertices.txt'), nodes, 'delimiter', ' ', 'precision', 16);
dlmwrite(strcat(cesta, 'elements.txt'), elements, 'delimiter', ' ');
dlmwrite(strcat(cesta, 'neighbors.txt'), neighbours, 'delimiter', ' ');
dlmwrite(strcat(cesta, 'elementType.txt'), type, 'delimiter', ' ');

