num_modes = 2;

% Variables
t = msspoly('t', 1);
x = msspoly('x', 2);

F = cell(num_modes, 1);
hX = cell(num_modes, 1);

% Dynamics
F{1} = [-2*x(2); 0.8*x(1) + 10*(x(1)^2 - 0.21)*x(2)];
F{2} = [-2*x(2); 0.8*x(1) + 10*(x(1)^2 - 0.21)*x(2)];

% Target region
hXT = [0.3 - x'*x];

% Boundary
hB = [1.3 - x'*x];

% Mode boundary
hX{1} = [x(1)];
hx{2} = [-x(1)];

% Parameters
T = 1;
d = 12;

% Integration options
options.sphereVars = x;
options.Adiag = [1.3, 1.3];
options.boxVars = [];
options.boxLims = [];

% Inner Approximation
innerApproximationNew(t, x, F, hX, hXT, hB, T, d, options);
