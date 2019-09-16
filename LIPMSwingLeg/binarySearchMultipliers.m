function [mult_opt,bmult_opt,smult_opt,mult_r_opt,bmult_r_opt] = binarySearchMultipliers(x, f, g, r, umat, V, S, B, outer_radius, inner_radius, T, max_iter)
if nargin < 12
    max_iter = 1;
end
nU = size(g,2);

% for i=1:length(S)
%   [~,~,coeff]=decomp(S(i));
%   mult_scale = norm(coeff,inf);
%   S(i) = S(i)/mult_scale;
% end
% 
% for i=1:length(B)
%   [~,~,coeff]=decomp(B(i));
%   mult_scale = norm(coeff,inf);
%   B(i) = B(i)/mult_scale;
% end

for i=1:max_iter
  %%
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newFree(1);
  
  Vr_sos = 1 - subs(V, x, r);
  [prog, Vr_sos, mult_r, ~] = spotless_add_eq_sprocedure(prog, Vr_sos, 1 - V, x, 4);
  [prog, Vr_sos, bmult_r, ~] = spotless_add_sprocedure(prog, Vr_sos, B, x, 0, 4);
  [prog, Vr_sos] = spotless_add_sprocedure(prog, Vr_sos, outer_radius-x'*x, x, 0, 4);
  prog = prog.withSOS(Vr_sos);
  
  mult = cell(1, 2^nU);
  bmult = cell(1, 2^nU);
  smult = cell(2^nU, nU);
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    [prog, Vdot_sos, mult{j}, ~] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V,x,4);
    [prog, Vdot_sos, bmult{j}, ~] = spotless_add_sprocedure(prog, Vdot_sos, -B, x, 0, 4);
    for k=1:nU
      [prog, Vdot_sos, smult{j}{k}, ~] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*S(k),x,0,4);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,0,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = 1;
  spot_options.sos_slack = -1e-6;
  spot_options.clean_primal = false;
  solver = @spot_mosek;
%   solver = @spot_sedumi;
  sol = prog.minimize(gamma,solver,spot_options);
  
  if double(sol.eval(gamma)) < -1e-6 % && sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    for j=1:2^nU
      mult_opt{j} = sol.eval(mult{j});
      bmult_opt{j} = sol.eval(bmult{j});
      for k=1:nU
        smult_opt{j}{k} = sol.eval(smult{j}{k});
      end
    end
    mult_r_opt = sol.eval(mult_r);
    bmult_r_opt = sol.eval(bmult_r);
  else
    keyboard
  end
end
end