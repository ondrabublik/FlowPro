function scramjet
% okrajove a pocatecni podminky: M = 5, r = 1, u = 1, v = 0

x = [0 3.9899 3.9899 0];
y = [-0.865 -0.3986 0.3986 0.865];

data = [x',y'];

save 'sj1' data

x = [1.088 2.1071 2.1998 3.394 2.941];
y = [0.346 0.1211 0.1211 0.2811 0.346];

data = [x',y'];

save 'sj2' data

y = -y;

data = [x',y'];

save 'sj3' data