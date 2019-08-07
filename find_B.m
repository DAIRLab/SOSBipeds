function [B] = find_B(t, x, u, f, g, r, V, S)
V_r = subs(V, x, r);

nU = size(g,2);
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, gamma] = prog.newPos(1);

[prog, B] = prog.newFreePoly(monomials(x, 0:4), 1);

% SOS 1
% V = 1 and V.R >= 1 => B <= -gamma
B_sos = (1 + (x'*x))*(-B - gamma);
[prog, B_sos] = spotless_add_eq_sprocedure(prog, B_sos, 1 - V, x, 4);
[prog, B_sos] = spotless_add_sprocedure(prog, B_sos, V_r - 1, x, 0, 4);
[prog, B_sos] = spotless_add_sprocedure(prog, B_sos, 1 - x'*x, x, 0, 4);
prog = prog.withSOS(B_sos);

% SOS 2
% V = 1 and \dot{V} >= 0 => B >= gamma
% divided into 2^nU parts
for j=1:2^nU
    B_sos_vdot = (1 + (x'*x))*(B - gamma);
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    [prog, B_sos_vdot] = spotless_add_eq_sprocedure(prog, B_sos_vdot, 1 - V, x, 4);
    [prog, B_sos_vdot] = spotless_add_sprocedure(prog, B_sos_vdot, Vdot, x, 0, 4);
    for k=1:nU
      [prog, B_sos_vdot] = spotless_add_sprocedure(prog, B_sos_vdot, umat(j,k)*S(k),x,0,4);
    end
    [prog, B_sos_vdot] = spotless_add_sprocedure(prog, B_sos_vdot, 1-x'*x,x,0,4);
    prog = prog.withSOS(B_sos_vdot);
end
 
spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = true;
solver = @spot_mosek;
% solver = @spot_sedumi;
sol = prog.minimize(gamma, solver, spot_options);

if sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    B = sol.eval(B);
else
    disp("Unsuccessful find_B");
    keyboard
end
end