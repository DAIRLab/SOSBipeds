function lipmSwingLegNStepCapturability(n)
% Outer approximation

if nargin < 1
  n = 0;
end

% Constants
g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = 0.05;
u_wrt_cm = 0;

model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm);

% Ellipsoid within which the states should lie
R_diag = [1 1 1];

% Increased step time for n=0
if n > 0
  T = step_time;
else
  T = 5;
end

%TODO: Remove unwanted options
options.degree = 6;
options.scale = 0.5;
options.control_design = false;
options.korda_control_design = false;
options.beta = 0;
options.infinite_time = false;
options.free_final_time = false;
options.swing_leg = true;

% Goal region for n=0
target = [];
% If goal region is more than just the origin use the following two lines
% goal_radius = 0.1;
% target = @(x) goal_radius^2 - x'*x;

[Vsol, Wsol] = nStepCapturabilitySOSwLeg(model, T, R_diag, target, n, options);

end
