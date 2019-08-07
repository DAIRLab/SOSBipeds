load('V1_inner_safe.mat', 'model', 'V', 'S', 'W');

nX = model.num_states;
nU = model.num_inputs;

x = msspoly('x', nX);
t = msspoly('t', 1);
u = msspoly('u', nU);

% Switching controller possibilities
ndgrid_arg = mat2cell(repmat([-1;1],1,nU), 2, ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

r = model.reset(t, x, []);
[f, g] = model.controlAffineDynamics(t, x);

% W = find_W(t, x, r, V);

% Plot the 0-step lyapunov function along with the reset map
plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], [0, 0, 0, 0; 0, 0.1, 0.2, 0.3], [1, 1], ['r', 'g']);
load('V1_LIPMSwingLeg.mat', 'Vsol');
plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], [0, 0, 0, 0; 0, 0.1, 0.2, 0.3], [0], ['m']);

iterations = 0;
has_converged = false;

radius.outer = 1;
radius.inner = 0.2;

% Initialize multipliers
[multiplier_v, multiplier_s] = find_multipliers(t, x, f, g, umat, V, S, W, radius);
S(2) = x(1) + x(2)/sqrt(model.g) - x(3);
while ~has_converged && iterations < 30
    iterations = iterations + 1;
    [multiplier_v, S, is_valid] = find_controller(t, x, f, g, umat, V, W, multiplier_s, radius);
    if ~is_valid
        has_converged = true;
        break;
    else
        [multiplier_s, V] = find_lyapunov_function(t, x, f, g, umat, V, S, W, multiplier_v, radius);
    end
    
    if mod(iterations, 2)
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], [0, 0, 0, 0; 0, 0.1, 0.2, 0.3], [1, 1], ['r', 'g']); 
    else
        plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
            [t; x(3)], [0, 0, 0, 0; 0, 0.1, 0.2, 0.3], [1, 1], ['b', 'g']); 
    end
end

function [multiplier_v, S, is_valid] = find_controller(t, x, f, g, umat, V, W, multiplier_s, radius)
    nU = size(g, 2);
    
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog, gamma] = prog.newFree(1);

    [prog, S] = prog.newFreePoly(monomials(x, 0:2), nU);

    multiplier_v = cell(1, 2^nU);
    for j=1:2^nU
        Vdot = diff(V, x)*(f + g*umat(j, :)');
        [prog, Vdot_sos, multiplier_v{j}, ~] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V, x, 4);
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
        for k=1:nU
            Vdot_sos = Vdot_sos - multiplier_s{j}{k}*umat(j, k)*S(k);
        end
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, radius.outer-x'*x, x, 0, 4);
        prog = prog.withSOS(Vdot_sos+gamma);
    end

    prog = prog.withPos(-gamma);

    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    sol = prog.minimize(gamma, @spot_mosek, spot_options);

    fprintf("In find_controller: \nStatus: %s \nGamma: %f\n", sol.status, double(sol.eval(gamma)));
    if double(sol.eval(gamma)) < -1e-6
        S = sol.eval(S);
        for j=1:2^nU
            multiplier_v{j} = sol.eval(multiplier_v{j});
        end
        is_valid = true;
    else
        keyboard
        % Set S and multiplier_v before returning
        is_valid = false;
    end
end

function [multiplier_s, V] = find_lyapunov_function(t, x, f, g, umat, V_init, S, W, multiplier_v, radius)
nX = size(x, 1);
nU = size(g, 2);

iterations = 6;
for i=1:iterations
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newFree(1);

    [prog,Q] = prog.newPSD(nX);
    V = x'*Q*x;

    multiplier_s = cell(2^nU, nU);
    for j=1:2^nU
        Vdot = diff(V,x)*(f + g*umat(j,:)');
        Vdot_sos = -Vdot - multiplier_v{j}*(1 - V);
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
        for k=1:nU
          [prog, Vdot_sos, multiplier_s{j}{k}, ~] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*S(k), x, 0, 4);
        end
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, radius.outer-x'*x,x,0,4);
        prog = prog.withSOS(Vdot_sos+gamma);
    end

    cost = spotlessIntegral(prog, V, x, ones(1, nX), [], []);
    
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
    
    disp(sol.status)
    disp(double(sol.eval(gamma)))
    
    if double(sol.eval(gamma)) < -1e-6 % && sol.status == spotsolstatus.STATUS_PRIMAL_AND_DUAL_FEASIBLE
        cost_max = cost_val;
        V_opt = sol.eval(V);
        multiplier_s_opt = cell(2^nU, nU);
        for j=1:2^nU
            for k=1:nU
                multiplier_s_opt{j}{k} = sol.eval(multiplier_s{j}{k});
            end
        end
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
multiplier_s = multiplier_s_opt;
end

function [multiplier_v, multiplier_s] = find_multipliers(t, x, f, g, umat, V, S, W, radius)
    nU = size(g, 2);

    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog, gamma] = prog.newFree(1);

    multiplier_v = cell(1, 2^nU);
    multiplier_s = cell(2^nU, nU);
    % (W >= 1 and V = 1) => \dot{V} <= 0
    for j=1:2^nU
        Vdot = diff(V, x)*(f + g*umat(j, :)');
        [prog, Vdot_sos, multiplier_v{j}] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V, x, 4);
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
        for k=1:nU
          [prog, Vdot_sos, multiplier_s{j}{k}] = spotless_add_sprocedure(prog, Vdot_sos, umat(j,k)*S(k), x, 0, 4);
        end
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, radius.outer-x'*x, x, 0, 4);
        prog = prog.withSOS(Vdot_sos+gamma);
    end

    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    sol = prog.minimize(gamma, @spot_mosek, spot_options);

    fprintf("In multipliers: \nStatus: %s \nGamma: %f\n", sol.status, double(sol.eval(gamma)));
    if double(sol.eval(gamma)) < 1e-6
        for j=1:2^nU
            multiplier_v{j} = sol.eval(multiplier_v{j});
            for k=1:nU
                multiplier_s{j}{k} = sol.eval(multiplier_s{j}{k});
            end
        end
    else 
        keyboard
        for j=1:2^nU
            multiplier_v{j} = NaN;
            for k=1:nU
                multiplier_s{j}{k} = NaN;
            end
        end
    end
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