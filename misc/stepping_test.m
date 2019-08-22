load('V0_inner_LIPMSwingLeg.mat', 'V', 'model')

x = msspoly('x', model.num_states);
t = msspoly('t', 1);
u = msspoly('u', model.num_inputs);

r = model.reset(t, x, []);

V_r = subs(V, x, r);

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, W] = prog.newSOSPoly(monomials(x, 2:2), 1);

% V_sos = (1 + x'*x)*(1 - W);
% [prog, V_sos, mult, ~] = spotless_add_sprocedure(prog, V_sos, 1 - V, x, 0, 4);
% [prog, V_sos, mult, ~] = spotless_add_sprocedure(prog, V_sos, 1 - x'*x, x, 0, 4);
% prog = prog.withSOS(V_sos);

Vr_sos = (1 + x'*x)*(W - 1);
[prog, Vr_sos, mult, ~] = spotless_add_eq_sprocedure(prog, Vr_sos, V_r - 1, x, 4);
[prog, Vr_sos, mult, ~] = spotless_add_sprocedure(prog, Vr_sos, 1 - x'*x, x, 0, 4);
prog = prog.withSOS(Vr_sos);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
solver = @spot_mosek;
%   solver = @spot_sedumi;
sol = prog.minimize(spotlessIntegral(prog, W, x, ones(1, model.num_states), [], []), solver, spot_options);

contourSpotless(sol.eval(W), x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'r'});
hold on;
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'b'});
load('V1_LIPMSwingLeg.mat', 'Vsol')
contourSpotless(Vsol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 0, {'m'});