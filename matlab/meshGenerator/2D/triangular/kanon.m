function data = kanon

delka = 8;
sirka = 3;
polomer = 0.5;
polomer2 = 0.3;
tlouska = 0.1;
a = 0.5;
b = 1;

% 0 -sirka/2;
% delka -sirka/2;

body = [delka 0;
        delka sirka/2;
        0 sirka/2;
        0 polomer+tlouska;
        a+tlouska polomer+tlouska;
        a+tlouska polomer2;
        a polomer2;
        a polomer;
        -b polomer;
        -b 0];
pocet = size(body,1);
    
data = zeros(2*pocet-1, 2);
data(pocet:end,:) = body(:,:);
for i = 1 : size(body,1)
    data(i,1) = body(end-i+1,1);
    data(i,2) = -body(end-i+1,2);
end