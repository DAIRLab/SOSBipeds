function [] = outerApproximation(t, x, F, hX, sX, R, hXT, d, options)

% time variables
T = 1;
hT = t*(T - t);

% empty R => Identityt map
if isempty(R)
    R = cell(nmodes,nmodes);
    for i = 1 : nmodes
    for j = 1 : nmodes
        if ~isempty(sX{i,j}), R{i,j} = x{i}; end
    end
    end
end

num_modes = size(x, 1);

% setup SOS program
prog = spotsosprog;
prog = prog.withIndeterminate(t);

for i = 1:num_modes
    prog = prog.withIndeterminate(x{i});
    
    % create v and w
    vmonom{ i } = monomials( [ t; x{ i } ], 0:d );
    [ prog, v{ i }, ~ ] = prog.newFreePoly( vmonom{ i } );
    wmonom{ i } = monomials( [ x{ i } ], 0:d );
    [ prog, w{ i }, ~ ] = prog.newFreePoly( wmonom{ i } );
    
    % constraints variables
    v0{i} = subs(v{i}, t, 0);
    vT{i} = subs(v{i}, t, T);
    dvdt{i} = diff(v{i}, t);
    dvdx{i} =  diff(v{i}, x{i});
    Lv{i} = dvdt{i} + dvdx{i}*F{i};
end

obj = 0;
for i = 1:num_modes
    % Lv_i(t, x) <= 0
    prog = sosOnK(prog, -Lv{i}, [t; x{i}], [hT; hX{i}], d); 
    
    % v_i(T, x) >= 0
    prog = sosOnK(prog, vT{i}, x{i}, hXT{i}, d);
    
    % w_i(x) > v_i(0, x) + 1
    prog = sosOnK(prog, w{i} - v0{i} - 1, x{i}, hX{i}, d);
    
    % w_i(x) >= 0
    prog = sosOnK(prog, w{i}, x{i}, hX{i}, d);
    
    % v_i'(t, R_ii'(x)) >= v_i(t, x)
    for j = 1:num_modes
       if(~isempty(sX{i, j}))
           vj_helper = subs( v{j}, x{j}, R{i, j});
           prog = sonOnK(prog, vj_helper - v{i}, [t; x{i}], [hT; sX{i, j}], d);
       end
    end
    
    % Objective function
    shape.sphere_vars = x{i};
    shape.A_diag = 1*eye(size(x{i}, 1));
    shape.box_vars = [];
    shape.box_lims = [];
    ld = generateMoments(prog, size(x{i}, 1), d, shape);
    obj = obj + wmonom{i}' * ld;
end

spot_options = spot_sdp_default_options();
spot_options.verbose = 1;

if isfield(options, 'solver_options')
    spot_options.solver_options = options.solver_options;
end

sol = prog.minimize(obj, @spot_mosek, spot_options);
keyboard
end