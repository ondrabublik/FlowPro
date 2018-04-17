function xs = getCentre(alpha)
global profileNURBS

profileNURBSdef = profileNURBS;
dx = alpha(1:profileNURBSdef.nc);
dy = alpha(profileNURBSdef.nc+[1:profileNURBSdef.nc]);

profileNURBSdef = profileNURBS;
profileNURBSdef.xc = profileNURBS.xc + dx;
profileNURBSdef.yc = profileNURBS.yc + dy;
data = nurbsClosed(profileNURBSdef);

xs = [sum(data(:,1)),sum(data(:,2))]/size(data,1);
    