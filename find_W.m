function [W] = find_W(t, x, r, V)
% Lyapunov function mapped through reset
V_r = subs(V, x, r);

% Setup SOS Program and variables
prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, W] = prog.newFreePoly(monomials(x, 2:2), 1); % Should this be a SOS polynomial?

% SOS polynomial
% (V_r >= 1, x'*x <= 1) => (W >= 1)
Vr_sos = (1 + x'*x)*(W - 1);
[prog, Vr_sos, ~, ~] = spotless_add_sprocedure(prog, Vr_sos, V_r - 1, x, 0, 4);
[prog, Vr_sos, ~, ~] = spotless_add_sprocedure(prog, Vr_sos, 1 - x'*x, x, 0, 4);
prog = prog.withSOS(Vr_sos);

% Cost
% int W(x) dx over an ellipse
cost = spotlessIntegral(prog, W, x, ones(1, size(x, 1)), [], []);

% Solve the optimization
spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
solver = @spot_mosek;
sol = prog.minimize(cost, solver, spot_options);

W = sol.eval(W);
end


% % Lyapunov function mapped through reset
% V_r = subs(V, x, r);
% 
% % Setup SOS Program and variables
% prog = spotsosprog;
% prog = prog.withIndeterminate(x);
% [prog, W] = prog.newFreePoly(monomials(x, 2:2), 1); % Should this be a SOS polynomial?
% 
% % SOS polynomial
% % (V_r >= 1, x'*x <= 1) => (W >= 1)
% Vr_sos = (1 + x'*x)*(W - 1);
% [prog, Vr_sos, ~, ~] = spotless_add_sprocedure(prog, Vr_sos, V_r - 1, x, 0, 4);
% [prog, Vr_sos, ~, ~] = spotless_add_sprocedure(prog, Vr_sos, 1 - x'*x, x, 0, 4);
% prog = prog.withSOS(Vr_sos);
% 
% % Cost
% % int W(x) dx over an ellipse
% cost = spotlessIntegral(prog, W, x, ones(1, model.num_states), [], []);
% 
% % Solve the optimization
% spot_options = spotprog.defaultOptions;
% spot_options.verbose = 1;
% spot_options.sos_slack = -1e-6;
% spot_options.clean_primal = false;
% solver = @spot_mosek;
% sol = prog.minimize(cost, solver, spot_options);
% 
% W = sol.eval(W);