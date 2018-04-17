function deklaraceMatic(P,Q,E)
% P,Q,E charakterizuji geometrii site
global nq QP QQ QE QSx QSy QO ne Enx Eny EQ EP Etyp np PX PY PQpoc PQ Ptyp;

np = P{1};
PX = P{2};
PY = P{3};
PQpoc = P{4}; % pocet sousednich ctvercu bodu site
PQ = P{5};    % indexy sousednich ctvercu bodu
Ptyp = P{6};  % typ bodu site

nq = Q{1};
QP = Q{2}; % - indexy bodu tvoricich ctverec
QQ = Q{3}; % - indexy sousednich ctvercu
QE = Q{4}; % - indexy stran tvoricich ctverec
QSx = Q{5}; % - Sx
QSy = Q{6}; % - Sy
QO = Q{9}; %  - obsah ctverce

ne = E{1};
EP = E{2}; % - indexy bodu tvoricich hranu
EQ = E{3}; % - indexy sousednich ctvercu
Enx = E{4}; % - normaly strany
Eny = E{5};
Etyp = E{6}; % - typ hrany