function [sol, Wsol] = outerApproximationUnified(t, x, F, hX, hXT, hB, d, options)

% time variables
T = 1;
hT = t*(T - t);

num_modes = size(x, 1);

% setup SOS program
prog = spotsosprog;
prog = prog.withIndeterminate(t);
prog = prog.withIndeterminate(x);    


% create v and w
% Change made here
vmonom = monolist( [ t; x ], d );
[ prog, v, coeff_v ] = prog.newFreePoly( vmonom );
wmonom = monolist( [ x ], d );
[ prog, w, coeff_w] = prog.newFreePoly( wmonom );

v0 = subs(v, t, 0);
vT = subs(v, t, T);
dvdt = diff(v, t);
dvdx = diff(v, x);
    
for i = 1:num_modes
    Lv{i} = dvdt + dvdx*F{i};
end

for i = 1:num_modes
    % Lv_i(t, x) <= 0
    prog = sosOnK(prog, -Lv{i}, [t; x], [hT; hX{i}], d); 
end
% v_i(T, x) >= 0
prog = sosOnK(prog, vT, x, hXT, d);

% w_i(x) > v_i(0, x) + 1
prog = sosOnK(prog, w - v0 - 1, x, hB, d);

% w_i(x) >= 0
prog = sosOnK(prog, w, x, hB, d);

% Objective function
% because X is always set to be a unit box
nx = length(x);
boxlimit =  [-ones(1, nx); ones(1, nx)];
ld = getLebesgueMoments(d, boxlimit, 1);

obj = coeff_w'*ld; 


spot_options = spot_sdp_default_options();
spot_options.verbose = 1;

if isfield(options, 'solver_options')
    spot_options.solver_options = options.solver_options;
end

sol = prog.minimize(obj, @spot_mosek, spot_options);

Wsol = struct();

w_coeff =  sol.eval(coeff_w);
w_monom =  wmonom;

Wsol.w_coeff = w_coeff;
Wsol.w_monom = w_monom;
end
