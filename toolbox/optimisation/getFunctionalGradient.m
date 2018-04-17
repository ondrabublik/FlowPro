function dI = getFunctionalGradient(nAlpha,h)
global optimisationPath

% jacobian matrix
inp = load([optimisationPath,'J.txt']);
Jac = sparse(inp(:,2)+1,inp(:,1)+1,inp(:,3));

% RHS
dIdw = load([optimisationPath,'dIdw.txt']);

% solve dual
lambda = Jac\dIdw;

% dR/da
R = load([optimisationPath,'R0.txt']);
I = load([optimisationPath,'I0.txt']);

dI = zeros(nAlpha,1);
for i = 1:nAlpha
    Rh = load([optimisationPath,'R',num2str(i),'.txt']);
    dRda = (Rh-R)/h;
    Ih = load([optimisationPath,'I',num2str(i),'.txt']);
    dIda = (Ih-I)/h;
    dI(i) = lambda'*dRda + dIda;
end

printLambda(lambda);


function printLambda(lambda)
global optimisationPath

% jacobian matrix
inp = load([optimisationPath,'printMatrix.txt']);
I = inp(:,1)+1;
J = inp(:,2)+1;
H = inp(:,3);
P = sparse(I,J,H);

lambda = P*lambda;
nEqs = 4;
lambda = reshape(lambda,nEqs,length(lambda)/nEqs);

fid = fopen([optimisationPath,'lambda.txt'],'w');
fprintf(fid,'%6.9f %6.9f %6.9f %6.9f\n',lambda);
fclose(fid);





