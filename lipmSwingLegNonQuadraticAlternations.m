load('V0_inner_LIPMSwingLeg.mat', 'V', 'model', 'S');
global I_w
load('W_indicator_function.mat', 'I_w');
 
nX = model.num_states;
nU = model.num_inputs;

x = msspoly('x', nX);
t = msspoly('t', 1);
u = msspoly('u', nU);

r = model.reset(t, x, []);

load('V1_inner_safe.mat', 'V', 'S');
load('V1_inner_safe.mat', 'W');

display_array = [0, 0, 0, 0, 0, 0; 0, 0.1, 0.2, 0.25, 0.27, 0.3];
close all;
plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1, 1], ['r', 'g']);
load('V1_LIPMSwingLeg.mat', 'Vsol');
plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [0], ['m']);
plot_figures(1, x'*x, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1], ['k']);

radius.outer = 1;
radius.inner = 0.1;

T = eye(model.num_states);

% Controller initialization
S(1) = x(1) + x(2)/sqrt(model.g/model.z_nom);
S(2) = x(1) - x(3) +  x(2)/sqrt(model.g/model.z_nom);
rho = 1;
for i=1:50
    num_iters = 0;
    [I_v] = find_indicator_function(V, x, 10);
    [multipliers] = find_multipliers(t, x, model, V, W, S, radius, T, 1, I_v);
    while ~isstruct(multipliers) && num_iters < 10
        [multipliers] = find_multipliers(t, x, model, (1.02^(num_iters+1))*V, W, S, radius, T, 1);
        num_iters = num_iters + 1;
    end
    if num_iters == 10
        keyboard
    end
    [V, S] = find_functions(t, x, model, V, W, multipliers, radius, T, 6);
    
    if mod(i, 20) == 0
        close all;
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [1, 1], ['r', 'g']);
        plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], display_array, [0], ['m']);
        plot_figures(1, x'*x, [x(1), x(2)], [-1, 1], [-1, 1], ...
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

function [I_v] = find_indicator_function(V, x, deg)
prog = spotsosprog;
prog = prog.withIndeterminate(x);

[prog, rho] = prog.newSOSPoly(monomials(x, 0:deg), 1);
[prog, rho_sos] = spotless_add_eq_sprocedure(prog, (1 + x'*x)*(rho - 1), 1 - V, x, 6);
[prog, rho_sos] = spotless_add_sprocedure(prog, rho_sos, 1 - x'*x, x, 0, 6);
prog = prog.withSOS(rho_sos);

cost = spotlessIntegral(prog, rho, x, [1, 1, 1], [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
spot_options.scale_monomials = true;
solver = @spot_mosek;

sol = prog.minimize(cost, solver, spot_options);

I_v = sol.eval(rho);
end

function [multipliers] = find_multipliers(t, x, model, V, W, S, radius, scaling, iterations, I_v)
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
            x, 0);
        for k=1:nU
          [prog, Vdot_sos, smult{j}{k}, ~] = spotless_add_sprocedure(prog, Vdot_sos, ...
              umat(j,k)*S(k),x,0);
        end
        ellipsoid_matrix = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, ...
            radius.outer-x'*ellipsoid_matrix*x,x,0);
        prog = prog.withSOS(Vdot_sos+gamma);
    end
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    solver = @spot_mosek;
    sol = prog.minimize(gamma, solver, spot_options);
    
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
mult_deg = 6;

for i=1:iterations
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newFree(1);

    [prog, V] = prog.newSOSPoly(monomials(x, 2:8), 1);
    [prog, S] = prog.newFreePoly(monomials(x,0:4), nU);

    mult = multipliers.mult;
    wmult = multipliers.wmult;
    smult = multipliers.smult;
    for j=1:2^nU
        Vdot = diff(V,x)*(f + g*umat(j,:)');
        Vdot_sos = -Vdot - mult{j}*(1 - V);
        [prog, Vdot_sos, wmult{j}, ~] = spotless_add_sprocedure(prog, ...
            Vdot_sos, W - 1, x, 0);
        for k=1:nU
          Vdot_sos = Vdot_sos - smult{j}{k}*umat(j,k)*S(k);
        end
        ellipsoid_matrix = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, ...
            radius.outer-x'*ellipsoid_matrix*x,x,0);
%         [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, x'*x-radius.inner,x,0,4);
        prog = prog.withSOS(Vdot_sos+gamma);
    end
    
    scale_mat = eye(nX);
    if i==1
        % initialize cost
        [~,cost_const] = calc_integral_cost(prog,x,V_init,radius.outer,scale_mat);
        cost_init = cost_const;

        if sign(cost_init) > 0
          cost_mult = 1/1.02;
        else
          cost_mult = 1.02;
        end

        cost_min = -inf;
        cost_max = cost_init;
        cost_val = cost_init;
    end

    cost = calc_integral_cost(prog,x,V,radius.outer,scale_mat);
    
    prog = prog.withPos(cost_val - cost);
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    solver = @spot_mosek;
    sol = prog.minimize(gamma,solver,spot_options);
    
    dbtxt = sprintf('In functions : \n Iterations: %d/%d \n', i, iterations);
    dbtxt = dbtxt + sprintf("Status: %s \n", sol.status);
    dbtxt = dbtxt + sprintf("Gamma: %f \n", sol.eval(gamma));
    disp(dbtxt);
    
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

function [cost,cost_const] = calc_integral_cost(prog,x,V,radius,scale_mat)
global I_w;
nX = length(x);
A_diag = ones(1,nX)*radius;
cost = spotlessIntegral(prog, V*I_w, x, A_diag, [], []);

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