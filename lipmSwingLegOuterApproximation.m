% Create model
g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = .05; % set to 0 to get point foot model with no continuous inputs
u_wrt_com = 0;
model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_com);

omega_0 = sqrt(model.g/model.z_nom);

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

% Dynamics
for i = 1:num_modes
    x{i} = msspoly('x', model.num_states);
end
[f, g] = model.controlAffineDynamics(t, x{1});
F{1} = f + g*[1; 1];
[f, g] = model.controlAffineDynamics(t, x{2});
F{2} = f + g*[1; -1];
[f, g] = model.controlAffineDynamics(t, x{3});
F{3} = f + g*[-1; 1];
[f, g] = model.controlAffineDynamics(t, x{4});
F{4} = f + g*[-1; -1];

% Domains
y = x{1};
hB = 1 - y'*y;
hX{1} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB];
load('V0_inner_LIPMSwingLeg.mat', 'V');
x_0 = msspoly('x', model.num_states);
W = find_W(t, x{1}, model.reset(t, x{1}, []), subs(V, x_0, x{1}));
hXT{1} = [1 - W; hX{1}];
% sX{1, 2} = y(1) + y(2)/omega_0 - y(3);

y = x{2};
hB = 1 - y'*y;
hX{2} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB];
x_0 = msspoly('x', model.num_states);
W = find_W(t, x{1}, model.reset(t, x{2}, []), subs(V, x_0, x{2}));
hXT{2} = [1 - W; hX{2}];

y = x{3};
hB = 1 - y'*y;
hX{3} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB];
x_0 = msspoly('x', model.num_states);
W = find_W(t, x{3}, model.reset(t, x{3}, []), subs(V, x_0, x{3}));
hXT{3} = [1 - W; hX{3}];

y = x{4};
hB = 1 - y'*y;
hX{4} = [ y(1) + y(2)/omega_0; y(1) + y(2)/omega_0 - y(3); hB];
x_0 = msspoly('x', model.num_states);
W = find_W(t, x{4}, model.reset(t, x{4}, []), subs(V, x_0, x{4}));
hXT{4} = [1 - W; hX{4}];

options.freeFinalTime = 0;

outerApproximation(t, x, F, hX, sX, R, hXT, d, options);