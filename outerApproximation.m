function [sol] = outerApproximation(t, x, F, hX, sX, R, hXT, d, options)

% time variables
T = 1;
hT = t*(T - t);

num_modes = size(x, 1);

% empty R => Identityt map
if isempty(R)
    R = cell(num_modes,num_modes);
    for i = 1 : num_modes
    for j = 1 : num_modes
        if ~isempty(sX{i,j}), R{i,j} = x{i}; end
    end
    end
end

% setup SOS program
prog = spotsosprog;
prog = prog.withIndeterminate(t);

for i = 1:num_modes
    prog = prog.withIndeterminate(x{i});
    
    % create v and w
% Change made here
%   vmonom{ i } = monomials( [ t; x{ i } ], 0:d );
    vmonom{ i } = monolist( [ t; x{ i } ], d );

    [ prog, v{ i }, coeff_v{i} ] = prog.newFreePoly( vmonom{ i } );
%     wmonom{ i } = monomials( [ x{ i } ], 0:d );
    wmonom{ i } = monolist( [ x{ i } ], d );

    [ prog, w{ i }, coeff_w{i}] = prog.newFreePoly( wmonom{ i } );
    
    % constraints variables
    v0{i} = subs(v{i}, t, 0);
    vT{i} = subs(v{i}, t, T);
    dvdt{i} = diff(v{i}, t);
    dvdx{i} =  diff(v{i}, x{i});
    Lv{i} = dvdt{i} + dvdx{i}*F{i};
end

obj = 0;
for i = 1:num_modes
    nX = size(x{i},1);
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
           prog = sosOnK(prog, vj_helper - v{i}, [t; x{i}], [hT; sX{i, j}], d);
       end
    end
    
    % Objective function
    shape.sphere_vars = x{i};
    shape.A_diag = 1/options.outerRadius^2*eye(size(x{i}, 1));
    shape.box_vars = x{i};
    shape.box_lims = [-options.outerRadius*ones(1,nX);  options.outerRadius*ones(1,nX)];
    
% ld = generateMoments(prog, size(x{i}, 1), d, shape);
% Change made here
    ld = getLebesgueMoments(d,  shape.box_lims, 1);
    obj = obj + coeff_w{i}' * ld; 
end

spot_options = spot_sdp_default_options();
spot_options.verbose = 1;

if isfield(options, 'solver_options')
    spot_options.solver_options = options.solver_options;
end


sol = prog.minimize(obj, @spot_mosek, spot_options);
keyboard;
end