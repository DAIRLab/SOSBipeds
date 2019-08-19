clear;
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
d = 12;
num_modes = 5;

% Variables
t = msspoly('t', 1);
x = msspoly('x', model.num_states);
F = cell(num_modes, 1);
hX = cell(num_modes, 1);

% Dynamics
[f, g] = model.controlAffineDynamics(t, x);
F{1} = T*(f + g*[1; 1]);
F{2} = T*(f + g*[1; -1]);
F{3} = T*(f + g*[-1; 1]);
F{4} = T*(f + g*[-1; -1]);
F{5} = [0; 0; 0];

% Target set
load('V0_inner_LIPMSwingLeg.mat', 'V');
W = find_W(t, x, model.reset(t, x, []), V);
hXT = 1 - W;

% Domain of each mode
% Mode 1
hX{1} = [ x(1) + x(2)/omega_0; x(1) + x(2)/omega_0 - x(3); W - 1];

% Mode 2
hX{2} = [ x(1) + x(2)/omega_0; - x(1) - x(2)/omega_0 + x(3); W - 1];

% Mode 3
hX{3} = [ - x(1) - x(2)/omega_0; x(1) + x(2)/omega_0 - x(3); W - 1];

% Mode 4
hX{4} = [ - x(1) - x(2)/omega_0; - x(1) - x(2)/omega_0 + x(3); W - 1];

% Mode 5
hx{5} = [1- W];

% Boundary
hB = [1 - x'*x];

% Parameters
T = 1;
d = 6;

% Integration options
options.sphereVars = x;
options.Adiag = [1, 1, 1];
options.boxVars = [];
options.boxLims = [];

% Inner approximation
innerApproximationNew(t, x, F, hX, hXT, hB, T, d, options);


