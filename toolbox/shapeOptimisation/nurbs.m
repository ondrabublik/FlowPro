function data = nurbs(profile)

np = profile.np;
p = profile.p;
n = profile.nc;
px = profile.xc;
py = profile.yc;
w = profile.w;

t = linspace(0,1,np);
x = zeros(size(t));
y = zeros(size(t));
sum = zeros(size(t));
u = [zeros(1,p), linspace(0,1,n-p+1), ones(1,p)]; % knots vector
for i = 1:n
    bsp = zeros(size(t));
    for j = 1:length(t)
        bsp(j) = bspline_basis(i, p, u, t(j));
    end
    bsp = w(i)*bsp;
    x = x + px(i)*bsp;
    y = y + py(i)*bsp;
    sum = sum + bsp; 
end
x = x./sum;
y = y./sum;
data = [x',y'];


function y = bspline_basis(i,p,t,x)
N = zeros(size(t));
nt = length(t);
for j = 1:nt-1
    if(t(j+1) < t(end))
        if(x >= t(j)) && (x < t(j+1))
            N(j) = 1.0;
            break;
        end
    elseif(t(j+1) == t(end))
        N(j) = 1.0;
        break;
    end
end

for r = 1:p
    for j = 1:nt-1-r
        dn0 = x - t(j);
        dd0 = t(j+r) - t(j);
        dn1 = t(j+r+1) - x;
        dd1 = t(j+r+1) - t(j+1);
        if(dd0 ~= 0)
            N(j) = (dn0/dd0)*N(j);
        end
        if(dd1 ~= 0)
            N(j) = N(j) + (dn1/dd1)*N(j+1);
        end
    end
end
y = N(i);
