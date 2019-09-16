load('V0_inner_LIPMSwingLeg.mat', 'V', 'model');
load('V1_inner_safe.mat', 'W');
load('W_indicator_function.mat', 'I_w');
load('V_intersection_W_indicator_fn.mat', 'I_vw');
V0 = V;

x = msspoly('x', model.num_states);
t = msspoly('t', 1);

% contourSpotless(V0, x(1), x(2), [-1 1], [-1 1], [t; x(3)], [0; 0.1], 1, {'r'});
% hold on;
% contourSpotless(W, x(1), x(2), [-1 1], [-1 1], [t; x(3)], [0; 0.1], 1, {'g'});

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, V] = prog.newSOSPoly(monomials(x, 2:10), 1);

% V0 >= 1 && W >= 1 => V >= 1
[prog, V_sos] = spotless_add_sprocedure(prog, V - 1, V0 - 1, x, 0);
[prog, V_sos] = spotless_add_sprocedure(prog, V_sos, W - 1, x, 0);
[prog, V_sos] = spotless_add_sprocedure(prog, V_sos, 1-x'*x, x, 0);
prog = prog.withSOS(V_sos);

cost = spotlessIntegral(prog, V, x, [1 1 1], [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
spot_options.scale_monomials = true;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

display_array = [0, 0, 0, 0; 0, 0.05, 0.1, 0.2];
plot_figures(1, [V0, W, 1-x'*x], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1, 1, 0], ['r', 'g', 'k']);
plot_figures(1, sol.eval(V), [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, 1, 'b');



