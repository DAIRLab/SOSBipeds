clear;
% Create model
g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = .05; % set to 0 to get point foot model with no continuous inputs
u_wrt_com = 0;
model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_com);

omega_0 = sqrt(model.g/model.z_nom);

% Box constraints are used here
outerRadius = 10;
options.outerRadius = outerRadius;


% Constants
T = 10;
d = 6;
num_modes = 4;

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

% Dynamics
for i = 1:num_modes
    x{i} = msspoly('x', model.num_states);
end
[f, g] = model.controlAffineDynamics(t, x{1});
F{1} = T*(f + g*[1; 1]);
[f, g] = model.controlAffineDynamics(t, x{2});
F{2} = T*(f + g*[1; -1]);
[f, g] = model.controlAffineDynamics(t, x{3});
F{3} = T*(f + g*[-1; 1]);
[f, g] = model.controlAffineDynamics(t, x{4});
F{4} = T*(f + g*[-1; -1]);

load('V0_inner_LIPMSwingLeg.mat', 'V');
% x_0 = msspoly('x', model.num_states);
W = find_W(t, x{1}, model.reset(t, x{1}, []), V);%subs(V, x_0, x{1})

% Domains
y = x{1};
ny = size(y, 1);
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{1} = boxLimit;
hX{1} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB{1}];
hXT{1} = [1 - W; hX{1}];
sX{1, 2} = [y(1) + y(2)/omega_0 - y(3); -1*(y(1) + y(2)/omega_0 - y(3))];
sX{1, 3} = [y(1) + y(2)/omega_0; -1*(y(1) + y(2)/omega_0)];
sX{1, 4} = [...
    y(1) + y(2)/omega_0; 
    y(1) + y(2)/omega_0 - y(3); 
    -1*(y(1) + y(2)/omega_0); 
    -1*(y(1) + y(2)/omega_0 - y(3)) ];
R{1, 2} = y;
R{1, 3} = y;
R{1, 4} = y;

y = x{2};
ny = size(y, 1);
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{2} = boxLimit;
hX{2} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB{2}];
hXT{2} = [1 - W; hX{2}];
sX{2, 1} = [y(1) + y(2)/omega_0 - y(3); -1*(y(1) + y(2)/omega_0 - y(3))];
sX{2, 3} = [...
    y(1) + y(2)/omega_0; 
    y(1) + y(2)/omega_0 - y(3); 
    -1*(y(1) + y(2)/omega_0); 
    -1*(y(1) + y(2)/omega_0 - y(3))];
sX{2, 4} = [y(1) + y(2)/omega_0; -1*(y(1) + y(2)/omega_0)];
R{2, 1} = y;
R{2, 3} = y;
R{2, 4} = y;

y = x{3};
ny = size(y, 1);
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{3} = boxLimit;
hX{3} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB{3}];
hXT{3} = [1 - W; hX{3}];
sX{3, 1} = [y(1) + y(2)/omega_0; -1*(y(1) + y(2)/omega_0)];
sX{3, 2} = [...
    y(1) + y(2)/omega_0; 
    y(1) + y(2)/omega_0 - y(3); 
    -1*(y(1) + y(2)/omega_0); 
    -1*(y(1) + y(2)/omega_0 - y(3))];
sX{3, 4} = [y(1) + y(2)/omega_0 - y(3); -1*(y(1) + y(2)/omega_0 - y(3))];
R{3, 1} = y;
R{3, 2} = y;
R{3, 4} = y;

y = x{4};
ny = size(y, 1);
boxLimit = [];
for j = 1:ny
    boxLimit = [boxLimit; options.outerRadius  + y(j); options.outerRadius - y(j)];
end
hB{4} = boxLimit;
hX{4} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB{4}];
hXT{4} = [1 - W; hX{4}];
sX{4, 1} = [...
    y(1) + y(2)/omega_0; 
    y(1) + y(2)/omega_0 - y(3); 
    -1*(y(1) + y(2)/omega_0); 
    -1*(y(1) + y(2)/omega_0 - y(3))];
sX{4, 2} = [y(1) + y(2)/omega_0; -1*(y(1) + y(2)/omega_0)];
sX{4, 3} = [y(1) + y(2)/omega_0 - y(3); -1*(y(1) + y(2)/omega_0 - y(3))];
R{4, 1} = y;
R{4, 2} = y;
R{4, 3} = y;

options.freeFinalTime = 0;

outerApproximation(t, x, F, hX, sX, R, hXT, d, options);