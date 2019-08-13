clear; close all;
% Constants
T = 40; options.Tspan = T;
d = 12;
M = 4;

num_modes = 2;

% Variables
t = msspoly('t', 1);
x = cell(num_modes, 1);
z = cell(num_modes, 1);

F = cell(num_modes, 1);
hX = cell(num_modes, 1);
hXT = cell(num_modes, 1);
sX = cell(num_modes, num_modes);
R = cell(num_modes, num_modes);

% Set the boundary to be a box constraint 
hB = cell(num_modes, 1);

% Box constraints are used here
outerRadius = M;
options.outerRadius = outerRadius;

% Dynamics
mu = 1;
for i = 1:num_modes
    x{i} = msspoly('x', 2);
    z{i} = msspoly('z', 2);
    y = z{i};
    F{i} = T*[-y(2); -mu*(1- y(1)^2)*y(2) + y(1)];

    F{i} = subs(F{i}, y, M*x{i})/M; 
end

% Gobal box constraint
for i = 1:num_modes
y = z{i};
ny = size(y, 1);

% Construct hX
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{i} = boxLimit;
hB{i} = subs(boxLimit, y, M*x{i})/M;
hX{i} = [hB{i}];
end

% Mode 1
y = z{1};
ny = size(y, 1);
exprX = [ y(2)];
exprX = subs(exprX, y, M*x{1});
hX{1} = [hX{1}; exprX];

exprXT = [0.1 - y'*y];
exprXT = subs(exprXT, y, M*x{1});
hXT{1} = [exprXT; hX{1}];


% Mode 2
y = z{2};
ny = size(y, 1);
exprX = [ -y(2)];
exprX = subs(exprX, y, M*x{2});
hX{2} = [hX{2}; exprX];

exprXT = [0.1 - y'*y];
exprXT = subs(exprXT, y, M*x{2});
hXT{2} = [exprXT; hX{2}];

% Guard and reset
y = z{1};
sX{1,2} = [ y(2); -y(2)];
sX{1,2} = subs(sX{1,2}, y, M*x{1});

% R{1,2} = y;
% R{1,2} = subs(R{1,2}, y, M*x{1});
R{1,2} = x{1};

y = z{2};
sX{2,1} = [ y(2); -y(2)];
sX{2,1} = subs(sX{2,1}, y, M*x{2});

% R{2,1} = y;
% R{2,1} = subs(R{2,1}, y, M*x{2});
R{2,1} = x{2};

options.freeFinalTime = 0;

[sol, Wsol] = outerApproximationHybrid(t, x, F, hX, sX, R, hXT, d, options);

sol.status

left = -3*M;
right = 3*M;

for i = 1 : num_modes
figure;
nx = size(x{i});
% M = options.outerRadius;
plotW = Wsol.w_monom{i}'*Wsol.w_coeff{i}; plotW = subs(plotW, x{i}, z{i}/M);
plotX = z{i};
contourSpotless(plotW, plotX(1), plotX(2),[left right],[left right],[], [], 1);

hold on;

Tspan = options.Tspan; 
[t, y] = ode45(@(t, x) vanderPol(1, t, x), [0; Tspan], [0; 0.1]);
plot(y(:, 1), y(:, 2));

plot(left:0.1:right, zeros(size(left:0.1:right)));

% xlim([-3 3]); ylim([-3 3]);
title(['M = ', num2str(M), ' d = ', num2str(d)]);
end

figure;
for i = 1 : num_modes
nx = size(x{i});
% M = options.outerRadius;
plotW = Wsol.w_monom{i}'*Wsol.w_coeff{i}; plotW = subs(plotW, x{i}, z{i}/M);
plotX = z{i};
contourSpotless(plotW, plotX(1), plotX(2),[left right],[left right],[], [], 1);

hold on;

Tspan = options.Tspan; 
[t, y] = ode45(@(t, x) vanderPol(1, t, x), [0; Tspan], [0; 0.1]);
plot(y(:, 1), y(:, 2));

plot(left:0.1:right, zeros(size(left:0.1:right)));

% xlim([-3 3]); ylim([-3 3]);
title(['M = ', num2str(M), ' d = ', num2str(d)]);
end