load('V0_inner_LIPMSwingLeg.mat', 'V', 'model');
x = msspoly('x', model.num_states);
t = msspoly('t', 1);
u = msspoly('u', model.num_inputs);
r = model.reset(t, x, []);

V_r = subs(V, x, r);

prog = spotsosprog;
prog = prog.withIndeterminate(x);

[prog, W] = prog.newSOSPoly(monomials(x, 2:2), 1);

W_sos = 1 - W;
[prog, W_sos] = spotless_add_eq_sprocedure(prog, W_sos, 1 - V, x, 4);
[prog, W_sos] = spotless_add_eq_sprocedure(prog, W_sos, 1 - V_r, x, 4);
[prog, W_sos] = spotless_add_sprocedure(prog, W_sos, 1 - x'*x, x, 0, 4);
prog = prog.withSOS(W_sos);

% W_sos = W - 1;
% [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, V - 1, x, 0, 4);
% [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, 1 - x'*x, x, 0, 4);
% prog = prog.withSOS(W_sos);
 
% W_sos = W - 1;
% [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, V_r - 1, x, 0, 4);
% [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, 1 - x'*x, x, 0, 4);
% prog = prog.withSOS(W_sos);

cost = spotlessIntegral(prog, W, x, ones(1, size(x, 1)), [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(cost, solver, spot_options);

W = sol.eval(W);
keyboard
