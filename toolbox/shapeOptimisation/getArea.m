function A = getArea(alpha)
global profileNURBS

profileNURBSdef = profileNURBS;
dx = alpha(1:profileNURBSdef.nc);
dy = alpha(profileNURBSdef.nc+[1:profileNURBSdef.nc]);

profileNURBSdef = profileNURBS;
profileNURBSdef.xc = profileNURBS.xc + dx;
profileNURBSdef.yc = profileNURBS.yc + dy;
data = nurbsClosed(profileNURBSdef);

A = 0;
for i = 1:size(data,1)-1
    A = A + (data(i+1,2)-data(i,2))*(data(i,1) + data(i+1,1))/2;
end
A = abs(A);
    