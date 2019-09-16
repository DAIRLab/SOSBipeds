load('V0_inner_LIPMSwingLeg.mat', 'V', 'model', 'S');
global I_w
load('W_indicator_function.mat', 'I_w');
 
nX = model.num_states;
nU = model.num_inputs;

x = msspoly('x', nX);
t = msspoly('t', 1);
u = msspoly('u', nU);

r = model.reset(t, x, []);
 
W = find_W(t, x, r, V);

% Load this to evaluate the multipliers when the alternations fail (last
% iteration of alternations)
load('V1_inner_safe.mat', 'V', 'S', 'model', 'W');

display_array = [0, 0, 0, 0, 0, 0; 0, 0.1, 0.2, 0.25, 0.27, 0.3];
close all;
% Display for W <= 1
plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1, 1], ['r', 'g']);
load('V1_LIPMSwingLeg.mat', 'Vsol');
plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [0], ['m']);
plot_figures(1, x'*[1, 0, 0; 0, 1, 0; 0, 0, 4]*x, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1], ['k']);

radius.outer = 1;
radius.inner = 0.1;

Q_init = double(subs(diff(diff(V,x)',x)/2,x,zeros(nX,1)));
T = eye(model.num_states);

% Controller initialization
S(1) = x(1) + x(2)/sqrt(model.g/model.z_nom);
S(2) = x(1) - x(3) +  x(2)/sqrt(model.g/model.z_nom);
rho = 1;
for i=1:50
    num_iters = 0;
    [old_multipliers] = find_multipliers(t, x, model, V, W, S, radius, T, 1);
    while ~isstruct(old_multipliers) && num_iters < 10
        [old_multipliers] = find_multipliers(t, x, model, (1.02^(num_iters+1))*V, W, S, radius, T, 1);
        num_iters = num_iters + 1;
    end
    if num_iters == 10
        keyboard
    end
    multipliers = old_multipliers;
    [V, S] = find_functions(t, x, model, V, W, multipliers, radius, T, 6);
    
    if mod(i, 20) == 0
        close all;
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [1, 1], ['r', 'g']);
        plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [0], ['m']);
        plot_figures(1, x'*[1, 0, 0; 0, 1, 0; 0, 0, 4]*x, [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [1], ['k']);
    elseif mod(i, 2) == 0
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [1, 1], ['r', 'g']); 
    else
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [1, 1], ['b', 'g']); 
    end
end
keyboard;

function [multipliers] = find_multipliers(t, x, model, V, W, S, radius, scaling, iterations)
nX = model.num_states;
nU = model.num_inputs;

[f, g] = model.controlAffineDynamics(t, x);

% Switching controller possibilities
ndgrid_arg = mat2cell(repmat([-1;1],1,nU), 2, ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

mult_deg = 6;

for i=1:iterations
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newFree(1);

    mult = cell(1, 2^nU);
    wmult = cell(1, 2^nU);
    smult = cell(2^nU, nU);
    for j=1:2^nU
        Vdot = (diff(V,x)*(f + g*umat(j,:)'));
        [prog, Vdot_sos, mult{j}, ~] = spotless_add_eq_sprocedure(prog, -Vdot, ...
            1-V,x,mult_deg);
        [prog, Vdot_sos, wmult{j}, ~] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, ...
            x, 0, mult_deg);
        for k=1:nU
          [prog, Vdot_sos, smult{j}{k}, ~] = spotless_add_sprocedure(prog, Vdot_sos, ...
              umat(j,k)*S(k),x,0,mult_deg);
        end
        ellipsoid_matrix = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, ...
            radius.outer-x'*ellipsoid_matrix*x,x,0,mult_deg);
        prog = prog.withSOS(Vdot_sos+gamma);
    end
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    solver = @spot_mosek;
    sol = prog.minimize(gamma,solver,spot_options);
    
    disp('In Multipliers:')
    disp('status: ');
    disp(sol.status);
    disp('gamma: ');
    disp(double(sol.eval(gamma)));
    if sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE && double(sol.eval(gamma)) < 1e-2 % modified from -1e-6
        for j=1:2^nU
            mult_opt{j} = sol.eval(mult{j});
            wmult_opt{j} = sol.eval(wmult{j});
            for k=1:nU
                smult_opt{j}{k} = sol.eval(smult{j}{k});
            end
        end
        multipliers.mult = mult_opt;
        multipliers.wmult = wmult_opt;
        multipliers.smult = smult_opt;
    elseif sol.status == spotsolstatus.STATUS_NUMERICAL_PROBLEMS && double(sol.eval(gamma)) < 1e-2 % modified from -1e-6
        for j=1:2^nU
            mult_opt{j} = sol.eval(mult{j});
            wmult_opt{j} = sol.eval(wmult{j});
            for k=1:nU
                smult_opt{j}{k} = sol.eval(smult{j}{k});
            end
        end
        multipliers.mult = mult_opt;
        multipliers.wmult = wmult_opt;
        multipliers.smult = smult_opt;
    else
        multipliers = false;
    end
end
end

function [V, S] = find_functions(t, x, model, V_init, W, multipliers, radius, scaling, iterations)
cost_option = 2;

nX = model.num_states;
nU = model.num_inputs;
is_fail = true;

[f, g] = model.controlAffineDynamics(t, x);

% Switching controller possibilities
ndgrid_arg = mat2cell(repmat([-1;1],1,nU), 2, ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end 

Q_init = double(subs(diff(diff(V_init, x)', x)/2, x, zeros(nX, 1))); %#ok<UDIM>

mult_deg = 6;

for i=1:iterations
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newFree(1);

    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;
%     [prog, V] = prog.newSOSPoly(monomials(x, 2:2), 1);
    [prog, S] = prog.newFreePoly(monomials(x,0:2), nU);

    mult = multipliers.mult;
    wmult = multipliers.wmult;
    smult = multipliers.smult;
    for j=1:2^nU
        Vdot = diff(V,x)*(f + g*umat(j,:)');
        Vdot_sos = -Vdot - mult{j}*(1 - V);
        [prog, Vdot_sos, wmult{j}, ~] = spotless_add_sprocedure(prog, ...
            Vdot_sos, W - 1, x, 0, mult_deg);
        for k=1:nU
          Vdot_sos = Vdot_sos - smult{j}{k}*umat(j,k)*S(k);
        end
        ellipsoid_matrix = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, ...
            radius.outer-x'*ellipsoid_matrix*x,x,0,mult_deg);
%         [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, x'*x-radius.inner,x,0,4);
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
          [~,cost_const] = calc_integral_cost(prog,x,V_init,radius.outer,scale_mat);
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
        cost = calc_integral_cost(prog,x,V,radius.outer,scale_mat);
    end

    prog = prog.withPos(cost_val - cost);
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    solver = @spot_mosek;
    sol = prog.minimize(gamma,solver,spot_options);
    
    disp(sol.status)
    disp(double(sol.eval(gamma)))
    
    if double(sol.eval(gamma)) < -1e-3
        cost_max = cost_val;
        V_opt = sol.eval(V);
        S_opt = sol.eval(S);
        is_fail = false;
    else
        cost_min = cost_val;
        if is_fail
          keyboard
          cost_val = cost_val*1.02;
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
global I_w;
nX = length(x);
A_diag = ones(1,nX)*radius;
cost = spotlessIntegral(prog, V*I_w, x, A_diag, [], []);

% When V is not known this line throws an error. Debug this later.
% Q_init = double(subs(diff(diff(V, x)', x)/2, x, zeros(nX, 1))); %#ok<UDIM>
% cost = spotlessIntegral(prog, V, x, diag(Q_init)', [], []);

% add cost in x-direction
% the "1" isn't quite right,
% cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,1)),x(1),radius/norm(scale_mat(:,1)),[],[]);
% cost = cost + 5000*cost_line;
% 
% cost_line = spotlessIntegral(prog,subs(V,x,x(1)*scale_mat(:,3)),x(1),radius/norm(scale_mat(:,3)),[],[]);
% cost = cost + 1000*cost_line;

% cost_line = spotlessIntegral(prog, subs(V, [x(1)], [0]), [x(2); x(3)], ones(1, 2), [], []);
% cost = cost + 500*cost_line;

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