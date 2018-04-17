function showDeformation
[geometry, simulation] = loadArgs;
geometryPath = strcat('../../simulations/', geometry, '/');
meshPath = strcat(geometryPath, 'mesh/');
prof = load([meshPath,'profileNURBS.mat']);
profileNURBS = prof.profile;

load alpha;
nAnim = 300;
nt = size(alphaEvo,2);
t = linspace(1,nt,nAnim);
da = zeros(2*profileNURBS.nc,nAnim);
for i = 1:2*profileNURBS.nc
    da(i,:) = spline(1:nt,alphaEvo(i,:),t);
end

data = nurbsClosed(profileNURBS);
figure('color','w')
hold on
tisk = plot(data(:,1),data(:,2),'color','k','linewidth',2);
I = [1:profileNURBS.nc,1];
tisk2 = plot(profileNURBS.xc(I),profileNURBS.yc(I),'marker','s','markerfacecolor','r');
axis equal
axis([-1 0 -0.5 0.5])
box on;
grid on;
for op = 1:nAnim
    profileNURBSdef = profileNURBS;
    profileNURBSdef.xc = profileNURBS.xc + da(1:profileNURBSdef.nc,op);
    profileNURBSdef.yc = profileNURBS.yc + da(profileNURBSdef.nc+[1:profileNURBSdef.nc],op);
    data = nurbsClosed(profileNURBSdef);
    set(tisk,'Xdata',data(:,1),'Ydata',data(:,2));
    set(tisk2,'Xdata',profileNURBSdef.xc(I),'Ydata',profileNURBSdef.yc(I));
    pause(0.01);
end