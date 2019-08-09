% Constants
T = 60;
d = 6;
num_modes = 1;

% Variables
t = msspoly('t', 1);
x = cell(num_modes, 1);
F = cell(num_modes, 1);
hX = cell(num_modes, 1);
hXT = cell(num_modes, 1);
sX = cell(num_modes, num_modes);
R = cell(num_modes, num_modes);

% Set the boundary to be a box constraint 
hB = cell(num_modes, 1);

% Box constraints are used here
outerRadius = 2;
options.outerRadius = outerRadius;

% Dynamics
mu = 1;
for i = 1:num_modes
    x{i} = msspoly('x', 2);
end
y = x{1};
F{1} = T*[-y(2); -mu*(1- y(1)^2)*y(2) + y(1)];
% y = x{2};
% F{2} = T*[-y(2); -mu*(1- y(1)^2)*y(2) + y(1)];

% Domains
% Mode 1
y = x{1};
ny = size(y, 1);
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{1} = boxLimit;
hX{1} = [hB{1}];
hXT{1} = [0.1 - y'*y; hX{1}];
% sX{1, 2} = [y(1); -y(1)];
% R{1, 2} = y;
% 
% % Mode 2
% y = x{2};
% ny = size(y, 1);
% boxLimit = [];
% for j = 1:ny
%     boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
% end
% hB{2} = boxLimit;
% hX{2} = [-y(1); hB{1}];
% hXT{2} = [0.1 - y'*y; hX{1}];
% sX{2, 1} = [y(1); -y(1)];
% R{2, 1} = y;

options.freeFinalTime = 0;

sol = outerApproximation(t, x, F, hX, sX, R, hXT, d, options);
[t, y] = ode45(@(t, x) vanderPol(1, t, x), [0; 30], [0; 0.1]);
plot(y(:, 1), y(:, 2))