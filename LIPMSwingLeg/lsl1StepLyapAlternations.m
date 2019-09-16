function [V, S] = lsl1StepLyapAlternations(x,f,g,V0,S0,W,R)
if nargin < 7
  R = [];
end

degree = 2;

N = 1;
nX = length(x);
nU = size(g,2);

Q_init = double(subs(diff(diff(V0,x)',x)/2,x,zeros(nX,1)));

V = V0;
S = S0;
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end


if nargin < 5
  R = [];
end

% scale problem data
% y = T*x
% x = T^-1 y
y = msspoly('y',nX);
T = Q_init^(.5);
x_y = inv(T)*y;
f_y = T*subs(f,x,x_y);
g_y = T*subs(g,x,x_y);
V_y = subs(V,x,x_y);
W_y = subs(W, x, x_y);
S_y = subs(S,x,x_y);
R_y = cell(length(R),1);
for i=1:length(R)
  R_y{i} = inv(T)'*R{i}*inv(T);
end
Q_init_y = eye(nX);

%% iter 1
outer_radius = 1;
inner_radius = .2;
rho0 = 1;
% [rho,mult,bmult] = binarySearchRho(x,f,g,umat,V,B,outer_radius,inner_radius,rho0);
[rho_y,mult_y,bmult_y,B_scale] = binarySearchRho(y,f_y,g_y,umat,V_y,S_y,W_y,outer_radius,inner_radius,rho0,T);
V_y = V_y/rho_y;

%% iter 2
% [V,B] = binarySearchVandB(x,f,g,umat,mult,bmult,outer_radius,inner_radius,Q_init,R,eye(nX));
[V_y,S_y] = binarySearchVandS(y,f_y,g_y,umat,mult_y,bmult_y,outer_radius,inner_radius,V_y,Q_init_y,W_y,R_y,T,degree,B_scale);
V = subs(V_y,y,T*x);
S = subs(S_y,y,T*x);
end

function [rho_opt,mult_opt,bmult_opt,S] = binarySearchRho(x,f,g,umat,V,S,W,outer_radius,inner_radius,rho0,T)
max_iter = 1;
rho = rho0;
rho_min = -inf;
rho_max = inf;
nU = size(g,2);

% for i=1:length(S),
%   [~,~,coeff]=decomp(S(i));
%   mult_scale = norm(coeff,inf);
%   S(i) = S(i)/mult_scale;
% end

for i=1:max_iter
  %%
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  [prog,gamma] = prog.newFree(1);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, -Vdot, rho-V,x,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
    for k=1:nU
      [prog, Vdot_sos,bmult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*S(k),x,0,4);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*inv(T)'*inv(T)*x,x,0,4); %#ok<*MINV>
%     [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*x,x,0,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
  
  spot_options = spotprog.defaultOptions;
  spot_options.verbose = 1;
  spot_options.sos_slack = -1e-6;
  spot_options.clean_primal = false;
  solver = @spot_mosek;
%   solver = @spot_sedumi;
  sol = prog.minimize(gamma,solver,spot_options);
  
  disp(double(sol.eval(gamma)));
  %%
  if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    %       keyboard
  end
  
  if double(sol.eval(gamma)) < -1e-6% && sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    rho_min = rho;
    for j=1:2^nU
      mult_opt{j} = sol.eval(mult{j});
      for k=1:nU,
        bmult_opt{j}{k} = sol.eval(bmult{j}{k});
      end
    end
  else
    rho_max = rho;
  end
  
  if isinf(rho_max)
    rho = 1.01*rho;
  elseif isinf(rho_min)
    rho = .99*rho;
  else
    rho = (rho_max + rho_min)/2;
  end
end
if isinf(rho_min)
  keyboard % uh oh
else
  rho_opt = rho_min;
end
rho_opt = min(rho_opt,1); % reset rho, but use best multipliers
end


function [V,S] = binarySearchVandS(x,f,g,umat,mult,bmult,outer_radius,inner_radius,V_init,Q_init,W,R,T,degree,B_init)
%%
% cost_option
% 1 - determinant
% 2 - integral
cost_option = 2;
max_iter = 8;
is_fail = true;

rho = 1;

% for i=1:length(bmult),
%   for j=1:length(bmult{i}),
%     [~,~,coeff]=decomp(bmult{i}{j});
%     mult_scale = max(1,norm(coeff,inf));
%     bmult{i}{j} = bmult{i}{j}/mult_scale;
%   end
% end
nX = length(x);
nU = size(g,2);

% scale_mat = eye(nX);
scale_mat = T';

% c(Q) = f(Q) + g
% c'(Q) = f(Q) = c(Q) - g
% c(Q) < 0 <==> c'(Q) < -g

for i=1:max_iter
  %%
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  rho = 1;
  [prog,gamma] = prog.newFree(1);
  if degree == 2
    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;
  else
    [prog,V] = prog.newFreePoly(monomials(x,2:2));
    % could add s-procedure here, if needed
    prog = prog.withSOS(V);
    
    Q = subs(diff(diff(V,x)',x)/2,x,zeros(nX,1));
  end  
  
  [prog,S] = prog.newFreePoly(monomials(x,1:degree),nU);
  %   B = B0;
  %   R = [];
  if ~isempty(R)
    if ~iscell(R)
      R = {R};
    end
    for j=1:length(R),
      if degree == 2
        prog = prog.withPSD(Q-R{j});
      else
        [prog, R_sos] = spotless_add_sprocedure(prog, V - 1, x'*R{j}*x-1,x,0,degree-2);
        prog = prog.withSOS(R_sos);
      end
    end
  end
  
%     V = V_init;
%   B = B_init;
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    Vdot_sos = -Vdot - mult{j}*(rho-V);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
    for k=1:nU,
      Vdot_sos = Vdot_sos - bmult{j}{k}*umat(j,k)*S(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*inv(T)'*inv(T)*x,x,0,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*inv(T)'*inv(T)*x,x,0,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
  
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
      cost_mult = 1/1.02;
    else
      cost_mult = 1.02;
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
  disp(double(sol.eval(gamma)));
  
  if sol.status ~= spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    %       keyboard
  end
  
  if double(sol.eval(gamma)) < -1e-6% && sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
    cost_max = cost_val;
    V_opt = sol.eval(V);
    S_opt = sol.eval(S);
    Q_opt = sol.eval(Q);
    is_fail = false;
  else
    cost_min = cost_val;
    if is_fail,
      keyboard
      cost_val = cost_val*cost_mult;
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
cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,1)),x(1),radius/norm(scale_mat(:,1)),[],[]);
cost = cost + 5000*cost_line;

cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,3)),x(1),radius/norm(scale_mat(:,3)),[],[]);
cost = cost + 1000*cost_line;

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