function data = nurbsClosed(profile)

np = profile.np;
p = profile.p;
n = profile.nc;
px = profile.xc;
py = profile.yc;
w = profile.w;

px = [px; px(1); px(2); px(3)];
py = [py; py(1); py(2); py(3)];
w = [w; w(1); w(2); w(3)];
n = n + 3;

u = linspace(0,1,n-p+1 + 2*p); % knots vector
t = linspace(u(4),u(end-3),np);

x = zeros(size(t));
y = zeros(size(t));
sum = zeros(size(t));
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

function t = mylinspace(a,b,np)
t = linspace(-1,1,np);
t = atan(3*t);
t = t - t(1);
t = t/max(t);
t = a + (b-a)*t;







