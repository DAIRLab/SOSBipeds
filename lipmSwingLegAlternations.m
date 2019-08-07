%% 0 STEP (balancing) ROA
% Constants
g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = 0.05;
u_wrt_cm = 0;
x_leg = 0;

model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm);

% Get an initial quadratic Lyapunov candidate
x = msspoly('x',model.num_states);
t = msspoly('t',1);
u = msspoly('u',model.num_inputs);
system_dynamics = model.dynamics(t,x,u);
[f,g] = model.controlAffineDynamics(t,x);

A = double(subs(diff(system_dynamics,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
B = double(subs(diff(system_dynamics,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

Q = eye(model.num_states);
R = eye(model.num_inputs);
[K,Q] = lqr(A,B,Q,R);

if exist(solutionFilename(model, 0), 'file') && 0
    load('V0_inner_LIPMSwingLeg.mat', 'V', 'model', 'S');
%     contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; x_leg],1,{'r'});
    disp("V0 and Bu loaded.");
else
    %
    A_state = {};
    A_state{1} = diag([0, 0, 0]);%diag([1/.25^2;1/.25^2;1/.25^2]);
    V0 = x'*Q*x;
    B0 = -diff(V0,x)*B;
    [ V,Bu ] = lipmSwingLegLyapunovAlternations(x,f,g,100*V0,B0,A_state);

    figure(1)
    hold off
    contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t;x(3)],[0; x_leg],1,{'r'});
    hold on;
    xlim([-1 1]);
    ylim([-1 1]);
    filename_suffix = class(model);
    filename = sprintf(['V%d_' filename_suffix '.mat'], 0);
    data = load(filename);
    V0 = data.Vsol;
    contourSpotless(V0, x(1), x(2), [-1 1], [-1 1], [t; x(3)], [0; x_leg], [0 0], {'g'});
    %%
    for i=1:30
      [ V,Bu ] = lipmSwingLegLyapunovAlternations(x,f,g,V,Bu,A_state);


    %   figure(1)
      hold on
      if mod(i,2) == 0
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; x_leg],1,{'r'});
      else
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; x_leg],1,{'b'});
      end
      sqrt(1./diag(double(subs(diff(diff(subs(V,t,0),x)',x)/2,x,x*0))))
    end
    % keyboard

    S = Bu;
    filename = solutionFilename(model, 0);
    save(filename,'V','model','S');
    return
end

V0 = V;
S0 = S;
load('boundary_test.mat', 'B_plot');
B0 = B_plot;
r = model.reset(t, x, []);

for iterations=1:1
    degree = 2;

    N = 1;
    nX = length(x);
    nU = size(g,2);

    Q_init = double(subs(diff(diff(V0,x)',x)/2,x,zeros(nX,1)));

    V = V0;
    B = B0;
    S = S0;
    ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
    [ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
    umat = zeros(2^nU,nU);
    for i=1:nU
      umat(:,i) = ugrid{i}(:);
    end
    
    % scale problem data
    % y = T*x
    % x = T^-1 y
%     y = msspoly('y',nX);
%     T = Q_init^(.5);
%     x_y = inv(T)*y;
%     r_y = subs(r,x,x_y);
%     f_y = T*subs(f,x,x_y);
%     g_y = T*subs(g,x,x_y);
%     V_y = subs(V,x,x_y);
%     B_y = subs(B,x,x_y);
%     S_y = subs(S,x,x_y);
%     Q_init_y = eye(nX);
    
    outer_radius = 2;
    inner_radius = .1;
    
    for j=1:1
%         [mult, smult, bmult] = binarySearchMultipliers(V_y, S_y, B_y, umat, f_y, g_y, r_y, y, outer_radius);
        [mult, smult, bmult] = binarySearchMultipliers(V, S, B, umat, f, g, r, x, outer_radius);
        [V_y, ~, S_y] = binarySearchVBandS(mult, bmult, smult, umat, f_y, g_y, r_y, y, V_y, Q_init_y, outer_radius, inner_radius, degree);
        contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0)
    end
    keyboard
end

function [mult_opt, smult_opt, bmult_opt] = binarySearchMultipliers(V, S, B, umat, f, g, r, x, outer_radius)
    nU = size(g,2);
%     for i=1:length(S)
%       [~,~,coeff]=decomp(S(i));
%       mult_scale = norm(coeff,inf);
%       S(i) = S(i)/mult_scale;
%     end
%     
%     for i=1:length(B)
%       [~,~,coeff]=decomp(B(i));
%       mult_scale = norm(coeff,inf);
%       B(i) = B(i)/mult_scale;
%     end
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newPos(1);
%     [prog, gamma] = prog.newFree(1);

    V_sos = (1 + x'*x)*(1 - subs(V, x, r));
    [prog, V_sos, mult1, coeff] = spotless_add_eq_sprocedure(prog, V_sos, 1 - V, x, 4);
    [prog, V_sos, bmult1, coeff] = spotless_add_sprocedure(prog, V_sos, B, x, 0, 4);
    [prog, V_sos] = spotless_add_sprocedure(prog, V_sos, 1 - x'*x, x, 0, 4);
    prog = prog.withSOS(V_sos);
    
    bmult_opt = cell(1);
    for j=1:2^nU
        Vdot = diff(V,x)*(f + g*umat(j,:)');
        [prog, Vdot_sos,mult{j},coeff] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V, x, 4);
        [prog, Vdot_sos,bmult{j},coeff] = spotless_add_sprocedure(prog, Vdot_sos, -B, x, 0, 4);
        for k=1:nU
          [prog, Vdot_sos,smult{j}{k},coeff] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*S(k),x,0,4);
        end
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,0,4);
        prog = prog.withSOS(Vdot_sos+gamma);
    end
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-2;
    spot_options.clean_primal = false;
    solver = @spot_mosek;
%     solver = @spot_sedumi;
    sol = prog.minimize(gamma,solver,spot_options);
    
    if sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE && double(sol.eval(gamma)) < 1e-6
        for j=1:2^nU
          mult_opt{j} = sol.eval(mult{j});
          for k=1:nU,
            smult_opt{j}{k} = sol.eval(smult{j}{k});
          end
        end 
    else
        keyboard
        disp("Unsuccessful");
    end
end

function [V, B, S] = binarySearchVBandS(mult, bmult, smult, umat, f, g, r, x, V_init, Q_init, outer_radius, inner_radius, degree)
%%
% cost_option
% 1 - determinant
% 2 - integral
cost_option = 2;
max_iter = 6;
is_fail = true;

nX = length(x);
nU = size(g,2);

for i=1:max_iter
  prog = spotsosprog;
  prog = prog.withIndeterminate(x);
  rho = 1;
  [prog,gamma] = prog.newFree(1);
  if degree == 2
    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;
  else
    [prog,V] = prog.newFreePoly(monomials(x,2:2));
    prog = prog.withSOS(V);
    
    Q = subs(diff(diff(V,x)',x)/2,x,zeros(nX,1));
  end  
  
  [prog,S] = prog.newFreePoly(monomials(x,1:degree),nU);
  
  prog = prog.withPSD(Q);
  
%   V_sos = 1 - subs(V, x, r);
%   [prog, V_sos, mult1, coeff] = spotless_add_eq_sprocedure(prog, V_sos, 1 - V, x, 4);
%   [prog, V_sos, bmult1, coeff] = spotless_add_sprocedure(prog, V_sos, B, x, 0, 4);
%   prog = prog.withSOS(V_sos);
  
  for j=1:2^nU
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    Vdot_sos = -Vdot - mult{j}*(rho-V);
    for k=1:nU
      Vdot_sos = Vdot_sos - smult{j}{k}*umat(j,k)*S(k);
    end
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, outer_radius-x'*x,x,0,4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, -inner_radius+x'*x,x,0,4);
    prog = prog.withSOS(Vdot_sos+gamma);
  end
 
  scale_mat = eye(nX);
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
  %%
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
display(sprintf('Determinant from %f to %f, percent change %f',Q_init_det,Q_det,100 - 100*Q_det/Q_init_det));
V = V_opt;
S = S_opt;
B = 0;
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

function [filename] = solutionFilename(model, n)
    filename_suffix = class(model);
    filename = sprintf(['V%d_inner_' filename_suffix '.mat'], n);
end