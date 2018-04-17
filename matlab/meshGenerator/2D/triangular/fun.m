function h = fun(x,y)

n = length(x);
h = 10*ones(size(x));

% pro obtekani profilu
r0 = 2;    % radius jemne oblasti
r1 = 0.6;
x0 = 0.5;  % x-ova souradnice stredu
y0 = 0;  % y-ova souradnice stredu
h0 = 0.15;
h1 = 0.05;
for i = 1 : n
    radius = sqrt((x(i)-x0)^2 + (y(i)-y0)^2);
    if radius < r0
        h(i) = h0;
    end
    if radius < r1
        h(i) = h1;
    end
end

% for i = 1 : n
%     radius = abs(x(i)-0.5);
%     if radius < 0.3
%         h(i) = 0.005;
%     end
%     if radius < 0.15
%         h(i) = 0.001;
%     end
% end

% % %   % % %  for NACA0012 aerofoil, angle of attack = 2 [deg], inviscid  % % %   % % %
% r2 = 1;
% r3 = 2;
% x0 = 0.5;
% y0 = 0.25;
% for i = 1 : n
%     radius = sqrt((x(i)-x0)^2 + (y(i)-y0)^2);
%     if x(i) > 0.65 && x(i) < 0.9 && y(i) > 0 && y(i) < 1.2
%         h(i) = 0.015;
%     elseif radius < r2
%         h(i) = 0.05;
%     elseif radius < r3
%         h(i) = 0.1;
%     else
%         h(i) = 5;
%     end
% end

% % %   % % %  for NACA0012 aerofoil, angle of attack = 0, viscous  % % %   % % %
% d0 = 2;
% h1 = 0.02;
% h2 = 0.3;
% for i = 1 : n
%     distance = d0;
%     if x(i) < 0
%         distance = sqrt(x(i)^2 + y(i)^2);
%     elseif x(i) < 3
%         distance = abs(y(i));
%     end
%     
%     if distance < d0
%         h(i) = (h2-h1)/d0 * distance + h1;
%     else
%         h(i) = 0.5;
%     end
% end

% % %   % % %  for NACA0012 aerofoil, angle of attack = 0.12, inviscid  % % %   % % %
% d0 = 2;
% h1 = 0.02;
% h2 = 0.3;
% for i = 1 : n
%     if x(i) < 0
%         distance = sqrt(x(i)^2 + y(i)^2);
%     elseif x(i) < 1
%         distance = abs(y(i));
%     else
%         distance = sqrt((x(i)-1)^2 + y(i)^2);
%     end
%     
% %     if x(i) > 0.4 && x(i) < 0.65 && y(i) > 0 && y(i) < 1
% %         h(i) = 0.015;
% %     elseif x(i) > 0.2 && x(i) < 0.4 && y(i) > -0.2 && y(i) < 0
% %         h(i) = 0.015;
%     if distance < d0
%         h(i) = (h2-h1)/d0 * distance + h1;
%     else
%         h(i) = 1;
%     end
% end




% r2 = 1;
% r3 = 2;
% x0 = 0;
% y0 = 0;
% for i = 1 : n
%     radius = sqrt((x(i)-x0)^2 + (y(i)-y0)^2);
%     if radius <= r2 || (x(i) >= 0 && x(i) <= 3 && y(i) >= -r2 && y(i) <= r2)
%         h(i) = 0.1;
%     if radius <= r2 || (x(i) >= 0 && x(i) <= 3 && y(i) >= -r2 && y(i) <= r2)
%         h(i) = 0.1;
%     elseif radius <= r3 || (x(i) >= 0 && x(i) <= 5 && y(i) >= -r3 && y(i) <= r3)
%         h(i) = 0.25;
%     else
%         h(i) = 0.5;
%     end
% end

% n = length(x);
% h = 0.1*ones(n,1);
% I = 1:n;
% a = 1;
% b = 1;
% x0 = 1.5;
% y0 = 0;
% Hp = a*(x-x0).^2 + b*(y-y0).^2;
% 
% h(I) = 0.01 + Hp(I)/20;
