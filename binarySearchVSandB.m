function [V, S, B] = binarySearchVSandB(x, f, g, r, umat, V_init, Q_init, ...
    mult, bmult, smult, mult_r, bmult_r, outer_radius, inner_radius, T, degree)
% cost_option
% 1 - determinant
% 2 - integral
cost_option = 2;
max_iter = 1;
is_fail = true;

nX = length(x);
nU = size(g,2);

for i=1:max_iter
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newFree(1);
  if degree == 2
    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;
  else
    % Does newSOSPoly work better then free poly and then imposing
    % constraint?
    [prog,V] = prog.newSOSPoly(monomials(x,2:2));
    
    % Is Q PSD already?
    Q = subs(diff(diff(V,x)',x)/2,x,zeros(nX,1));
  end  
  
  [prog, S] = prog.newFreePoly(monomials(x,0:degree), nU);
  [prog, B] = prog.newFreePoly(monomials(x,0:degree), 1);
  
  Vr_sos = 1 - subs(V, x, r);
  Vr_sos = Vr_sos - mult_r*(1 - V);
  Vr_sos = Vr_sos - bmult_r*B;
  [prog, Vr_sos] = spotless_add_sprocedure(prog, Vr_sos, outer_radius - x'*x, x, 0, 4);
  prog = prog.withSOS(Vr_sos);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    Vdot_sos = -Vdot - mult{j}*(1-V);
    Vdot_sos = Vdot_sos - bmult{j}*(-B);
    for k=1:nU
      Vdot_sos = Vdot_sos - smult{j}{k}*umat(j,k)*S(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,0,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*x,x,0,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
 
  scale_mat = T;
  if i==1
    % initialize cost
    if cost_option == 1
      [~,cost_const] = calc_cost(Q,Q_init,scale_mat);
      
      % based on the fact that cost(V0) = 0
      % gets around the fact that subs(Q) doesn't work with Q PSD
      cost_init = -cost_const;
    else
      [~,cost_const] = calc_integral_cost(prog,x,V_init,outer_radius,scale_mat);
      cost_init = cost_const;
    end
    
    if sign(cost_init) > 0
      cost_mult = 1/1.03;
    else
      cost_mult = 1.03;
    end
    
    cost_min = -inf;
    cost_max = cost_init;
    cost_val = cost_init;
  end
  
  if cost_option == 1
    cost = calc_cost(Q,Q_init,scale_mat);
  else
    cost = calc_integral_cost(prog,x,V,outer_radius,scale_mat);
  end
  
  prog = prog.withPos(cost_val - cost);
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = 1;
  spot_options.sos_slack = -1e-6;
  spot_options.clean_primal = false;
  solver = @spot_mosek;
%   solver = @spot_sedumi;
  sol = prog.minimize(gamma,solver,spot_options);
  
  if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
%           keyboard
  end
  
  if sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE && double(sol.eval(gamma)) < -1e-6 
    cost_max = cost_val;
    V_opt = sol.eval(V);
    S_opt = sol.eval(S);
    B_opt = sol.eval(B);
    Q_opt = sol.eval(Q);
    is_fail = false;
  else
    cost_min = cost_val;
    if is_fail
      keyboard
      cost_val = cost_val*1.05;
    end
  end
  
  if ~is_fail
    if ~isinf(cost_min)
      cost_val = (cost_max + cost_min)/2;
    else
      cost_val = cost_val*cost_mult;
    end
  end
end

if is_fail
  keyboard %uh oh
end

Q_init_det = det(scale_mat*Q_init*scale_mat');
Q_det = det(scale_mat*Q_opt*scale_mat');
fprintf('Determinant from %f to %f, percent change %f\n',Q_init_det,Q_det,100 - 100*Q_det/Q_init_det);
V = V_opt;
S = S_opt;
B = B_opt;
end

function [cost,cost_const] = calc_cost(Q,Q_init,scale_mat)
Q_trans = scale_mat*Q*scale_mat';
Q_init_trans = scale_mat*Q_init*scale_mat';

cost_coeffs = det(Q_init_trans)*inv(Q_init_trans);

% add cost on Q(1,1)
cost_coeffs = cost_coeffs + 5*scale_mat(:,1)*scale_mat(1,:);

cost_coeffs = cost_coeffs/norm(cost_coeffs(:),inf);
cost = sum(sum((Q_trans - Q_init_trans).*cost_coeffs));


coeffs = cost.coeff;
pows = cost.pow;
assert(pows(1) == 0)

cost_const = coeffs(1);
cost = cost - cost_const;
end

function [cost,cost_const] = calc_integral_cost(prog,x,V,radius,scale_mat)
nX = length(x);
A_diag = ones(1,nX)*radius;
cost = spotlessIntegral(prog,V,x,A_diag,[],[]);

% add cost in x-direction
% the "1" isn't quite right,
% cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,1)),x(1),radius/norm(scale_mat(:,1)),[],[]);
% cost = cost + 5000*cost_line;
% 
% cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,3)),x(1),radius/norm(scale_mat(:,3)),[],[]);
% cost = cost + 1000*cost_line;

%
if isnumeric(cost)
  cost_const = cost;
  cost = 0;
else
  coeffs = cost.coeff;
  pows = cost.pow;
  
  if any(pows == 0)
    assert(pows(1) == 0)
    cost_const = coeffs(1);
    cost = cost - cost_const;
  else
    cost_const = 0;
  end
end
end