function minL = getMinimalThickness(alpha)
global profileNURBS

profileNURBSdef = profileNURBS;
dx = alpha(1:profileNURBSdef.nc);
dy = alpha(profileNURBSdef.nc+[1:profileNURBSdef.nc]);

profileNURBSdef = profileNURBS;
profileNURBSdef.xc = profileNURBS.xc + dx;
profileNURBSdef.yc = profileNURBS.yc + dy;
data = nurbsClosed(profileNURBSdef);
x = data(1:end-1,1);
y = data(1:end-1,2);

tri = delaunay(x,y);

nt = size(tri,1);
xs = zeros(nt,1);
ys = zeros(nt,1);

for i = 1:nt
    xs(i) = sum(x(tri(i,:)))/3;
    ys(i) = sum(y(tri(i,:)))/3;
end

in = inpolygon(xs,ys,x,y);

minL = 1e5;
for i = 1:nt
    if(in(i))
        a = sqrt((x(tri(i,2))-x(tri(i,1)))^2 + (y(tri(i,2))-y(tri(i,1)))^2);
        b = sqrt((x(tri(i,3))-x(tri(i,2)))^2 + (y(tri(i,3))-y(tri(i,2)))^2);
        c = sqrt((x(tri(i,1))-x(tri(i,3)))^2 + (y(tri(i,1))-y(tri(i,3)))^2);
        L = 2*sqrt(a*b*c/(a+b+c));
        if(L < minL)
            minL = L;
        end
    end
end
