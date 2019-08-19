function [sol, Wsol, Vsol] = innerApproximationNew(t, x, F, hX, hXT, hB, T, d, options)

num_modes = size(F, 1);

% Time constraints
hT = t*(T - t);

% Complement of target region
hXTc = -hXT;

%%  SOS Program and variables setup
prog = spotsosprog;
prog = prog.withIndeterminate(t);
prog = prog.withIndeterminate(x);

V_vars = [t; x];
[prog, V] = prog.newFreePoly(monomials(V_vars, 0:d));

W_vars = x;
[prog, W] = prog.newSOSPoly(monomials(W_vars, 0:d));

%% SOS Constraints

for m=1:num_modes
    % Vdot <= 0
    dVdt = diff(V, t);
    dVdx = diff(V, x);
    Vdot = dVdt + dVdx * F{m};
    [prog, VdotSOS] = spotless_add_sprocedure(prog, -Vdot, hT, V_vars, 0, d); 
    for i=hB
        [prog, VdotSOS] = spotless_add_sprocedure(prog, VdotSOS, i, V_vars, 0, d); 
    end
    for i=hX{m}
        [prog, VdotSOS] = spotless_add_sprocedure(prog, VdotSOS, i, V_vars, 0, d);
    end
    prog = prog.withSOS(VdotSOS);
end

% V(t, x) >= 0 on the boundary
Vsos = V;
for i=hB
    [prog, Vsos] = spotless_add_eq_sprocedure(prog, Vsos, i, V_vars, d);
end
[prog, Vsos] = spotless_add_sprocedure(prog, Vsos, hT, V_vars, 0, d);
prog = prog.withSOS(Vsos);

% V(T, x) >= 0
V_TSOS = subs(V, t, T);
for i=hB
    [prog, V_TSOS] = spotless_add_sprocedure(prog, V_TSOS, i, V_vars, 0, d);
end
[prog, V_TSOS] = spotless_add_sprocedure(prog, V_TSOS, hXTc, V_vars, 0, d);
prog = prog.withSOS(V_TSOS);

% W >= V(0, x) + 1
V_0 = subs(V, t, 0);
Wsos = W - V_0 - 1;
for i=hB
    [prog, Wsos] = spotless_add_sprocedure(prog, Wsos, i, W_vars, 0, d);
end
prog = prog.withSOS(Wsos);

%% Objective function
if ~isfield(options, 'sphereVars')
    options.sphereVars = [];
end

if ~isfield(options, 'Adiag')
    options.Adiag = [];
end

if ~isfield(options, 'boxVars')
    options.boxVars = [];
end

if ~isfield(options, 'boxLims')
    options.boxLims = [];
end

obj = spotlessIntegral(prog, W, options.sphereVars, options.Adiag, options.boxVars, options.boxLims);

%% Solve
spot_options = spotprog.defaultOptions;
spot_options.verbose = true;
spot_options.do_fr = false;
solver = @spot_mosek;

sol = prog.minimize(obj,solver,spot_options);

V_sol = sol.eval(V);
W_sol = sol.eval(W);

%% Plotting
figure();
contourSpotless(V_sol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0, {'r'});
hold on;
for i=hB
    contourSpotless(i, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0, {'k'});
end
% [tt, yy] = ode45(@(t, x) vanderPol(1, t, x), [0, 20], [0; 0.1]);
% plot(yy(:, 1), yy(:, 2), 'b');
contourSpotless(W_sol, x(1), x(2), [-1, 1], [-1, 1], [x(3)], [0], 1, {'g'});
contourSpotless(hXT, x(1), x(2), [-1, 1], [-1, 1], [x(3)], [0], 0, {'m'});
keyboard;

end