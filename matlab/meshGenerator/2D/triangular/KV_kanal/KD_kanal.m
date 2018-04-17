function KD_kanal
% konvergentni-divergentni kanal
y1 = 0;
y2 = 0.3;
y3 = 1;
x1 = 1;
x2 = 2;
x3 = 5;

x = [0,x1,x2,x3,x3,0];
y = [y1,y1,y2,y2,y3,y3];

data = [x',y'];

save 'KD_kanal' data