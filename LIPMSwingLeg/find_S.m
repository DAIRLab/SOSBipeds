function [S] = find_S(x, f, g, V, W)
    prog = spotsosprog;
    prog = prog.withIndeterminate(x);
    [prog, gamma] = prog.newFree(1);
    
    nU = 2;
    [prog, S] = prog.newFreePoly(monomials(x, 0:4), nU);

    Vdot = diff(V, x)*(f + g*S);
    [prog, Vdot_sos] = spotless_add_eq_sprocedure(prog, -Vdot, 1-V, x, 4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, 1-x'*x, x, 0, 4);
    prog = prog.withSOS(Vdot_sos+gamma);
    
    Vdot = diff(V, x)*(f + g*S) + 100;
    [prog, Vdot_sos] = spotless_add_eq_sprocedure(prog, Vdot, 1-V, x, 4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, W - 1, x, 0, 4);
    [prog, Vdot_sos] = spotless_add_sprocedure(prog, Vdot_sos, 1-x'*x, x, 0, 4);
    prog = prog.withSOS(Vdot_sos);

    prog = prog.withPos(-gamma);

    spot_options = spotprog.defaultOptions;
    spot_options.verbose = 1;
    spot_options.sos_slack = -1e-6;
    spot_options.clean_primal = false;
    spot_options.scale_monomials = true;
    sol = prog.minimize(gamma, @spot_mosek, spot_options);
    
    sol.eval(gamma)
    keyboard;
end