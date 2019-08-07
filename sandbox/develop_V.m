load('V1_inner_swig_leg.mat', 'V', 'W', 'S', 'model');

x = msspoly('x', 3);
t = msspoly('t', 1);

[f, g] = model.controlAffineDynamics(t, x);

nU = size(g,2);
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

prog = spotsosprog;
prog = prog.withIndeterminate(x);
[prog, gamma] = prog.newPos(1);

[prog, W1] = prog.newSOSPoly(monomials(x, 2:2), 1);
sos_list = [];
for j=1:2^nU
    W_sos = (1 + x'*x)*(1 - W1);
    Vdot = diff(V,x)*(f + g*umat(j,:)');
    [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, -Vdot, x, 0, 4);
    [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, W - 1, x, 0, 4);
    for k=1:nU
      [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, umat(j,k)*S(k),x,0,4);
    end
    [prog, W_sos] = spotless_add_sprocedure(prog, W_sos, 1-x'*x,x,0,4);
    prog = prog.withSOS(W_sos - gamma);
    sos_list = [sos_list, W_sos - gamma]; %#ok<AGROW>
end

% cost = spotlessIntegral(prog, W1, x, [1, 1, 1], [], []);

spot_options = spotprog.defaultOptions;
spot_options.verbose = 1;
spot_options.sos_slack = -1e-6;
spot_options.clean_primal = false;
solver = @spot_mosek;
%   solver = @spot_sedumi;
sol = prog.minimize(gamma,solver,spot_options);

double(sol.eval(gamma))

for i=sos_list
    disp(subs(sol.eval(i), x, [0; 0.28; -0.2]))
end

contourSpotless(sol.eval(W1), x(1), x(2), [-10, 10], [-10, 10], [t; x(3)], [0; 0.1], 1, {'r'});
hold on;
contourSpotless(V, x(1), x(2), [-10, 10], [-10, 10], [t; x(3)], [0; 0.1], 1, {'b'});
contourSpotless(W, x(1), x(2), [-10, 10], [-10, 10], [t; x(3)], [0; 0.1], 1, {'g'});
keyboard;