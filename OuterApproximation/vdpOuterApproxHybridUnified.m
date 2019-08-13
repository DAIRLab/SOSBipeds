clear; close all;
% Constants
T = 40; options.Tspan = T;
d = 14;
M = 4;

num_modes = 2;

% Variables
t = msspoly('t', 1);
x = msspoly('x', 2);
z = msspoly('z', 2);

F = cell(num_modes, 1);
hX = cell(num_modes, 1);

% Box constraints are used here
outerRadius = M;
options.outerRadius = outerRadius;

% Dynamics
mu = 1;
for i = 1:num_modes
    y = z;
    F{i} = T*[-y(2); -mu*(1- y(1)^2)*y(2) + y(1)];

    F{i} = subs(F{i}, y, M*x)/M; 
end

% Gobal box constraint
y = z;
ny = size(y, 1);

% Construct hX
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB = boxLimit;
hB = subs(boxLimit, y, M*x)/M;

y = z;
exprXT = [0.1 - y'*y];
exprXT = subs(exprXT, y, M*x);
hXT = [hB; exprXT];

% Mode 1
ny = size(y, 1);
exprX = [ y(2)];
exprX = subs(exprX, y, M*x);
hX{1} = [hB; exprX];

% Mode 2
exprX = [-y(2)];
exprX = subs(exprX, y, M*x);
hX{2} = [hB; exprX];

options.freeFinalTime = 0;

[sol, Wsol] = outerApproximationUnified(t, x, F, hX, hXT, hB, d, options);

sol.status

left = -M;
right = M;

figure;
nx = length(x);
% M = options.outerRadius;
plotW = Wsol.w_monom'*Wsol.w_coeff; plotW = subs(plotW, x, z/M);
plotX = z;
contourSpotless(plotW, plotX(1), plotX(2),[left right],[left right],[], [], 1);

hold on;

Tspan = options.Tspan; 
[t, y] = ode45(@(t, x) vanderPol(1, t, x), [0; Tspan], [0; 0.1]);
plot(y(:, 1), y(:, 2));

plot(left:0.1:right, zeros(size(left:0.1:right)));

% xlim([-3 3]); ylim([-3 3]);
title(['M = ', num2str(M), ' d = ', num2str(d)]);
