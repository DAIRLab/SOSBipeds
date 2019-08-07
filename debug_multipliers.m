function [] = debug_multipliers(model, V, S, W, multipliers)
    nX = model.num_states;
    nU = model.num_inputs;
    
    x = msspoly('x', nX);
    t = msspoly('t', 1);
    
    [f, g] = model.controlAffineDynamics(t, x);
    
    % Switching controller possibilities
    ndgrid_arg = mat2cell(repmat([-1;1],1,nU), 2, ones(1,nU)');
    [ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
    umat = zeros(2^nU,nU);
    for i=1:nU
      umat(:,i) = ugrid{i}(:);
    end

    mult_deg = 6;

    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog,gamma] = prog.newFree(1);

    mult = multipliers.mult;
    smult = multipliers.smult;
    for j=1:2^nU
        Vdot = diff(V,x)*(f + g*umat(j,:)');
        Vdot_sos = -Vdot - mult{j}*(1 - V);
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, ...
            Vdot_sos, W - 1, x, 0, mult_deg);
        for k=1:nU
          Vdot_sos = Vdot_sos - smult{j}{k}*umat(j,k)*S(k);
        end
        ellipsoid_matrix = [1, 0, 0; 0, 1, 0; 0, 0, 1];
        [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, ...
            1-x'*ellipsoid_matrix*x,x,0,mult_deg);
        prog = prog.withSOS(Vdot_sos+gamma);
    end
    
    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 0;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    solver = @spot_mosek;
    sol = prog.minimize(gamma,solver,spot_options);
    
    double(sol.eval(gamma))
end