load('V1_inner_swing_leg.mat', 'W', 'model');
load('V0_inner_LIPMSwingLeg.mat', 'V');
V0 = V;

x = msspoly('x', model.num_states);
t = msspoly('t', 1);

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, rho] = prog.newSOSPoly(monomials(x, 0:25), 1);

[prog, rho_sos] = spotless_add_sprocedure(prog, rho - 1, 1 - W, x, 0, 4);
% [prog, rho_sos] = spotless_add_sprocedure(prog, rho_sos, 1 - V0, x, 0, 4);
[prog, rho_sos] = spotless_add_sprocedure(prog, rho_sos, 1 - x'*x, x, 0, 4);
prog = prog.withSOS(rho_sos);

objective = spotlessIntegral(prog, rho, x, [1, 1, 1], [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
spot_options.scale_monomials = true;
solver = @spot_mosek;
sol = prog.minimize(objective, solver, spot_options);

display_array = [0, 0, 0, 0; 0, 0.033, 0.066, 1];
close all;
plot_figures(1, [V0, W, 1-x'*x], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1, 1, 0], ['r', 'g', 'k']);
plot_figures(1, [sol.eval(rho)], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1], ['b']);

keyboard
% I_w = sol.eval(rho);
% save('../data/W_indicator_function.mat', 'I_w');
