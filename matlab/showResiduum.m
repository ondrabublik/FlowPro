function showResiduum

[meshPath, simulPath, outputPath] = getPath;
residuum = load([simulPath,'residuum.txt']);

figure('color','w')
semilogy(residuum(:,1),'color','k','linewidth',3)
grid on
box on
ylabel('resid')
xlabel('iteration')

figure('color','w')
semilogy(residuum(:,2),residuum(:,1),'color','k','linewidth',3)
grid on
box on
ylabel('resid')
xlabel('time')

figure('color','w')
semilogy(reorganise(residuum(:,3)/1000),residuum(:,1),'color','k','linewidth',3)
grid on
box on
ylabel('resid')
xlabel('CPU [s]')

function x = reorganise(x)
n = length(x);
for i = 2:n
    if(x(i) < x(i-1))
        x(i:end) = x(i:end) + x(i-1)-x(i);
    end
end