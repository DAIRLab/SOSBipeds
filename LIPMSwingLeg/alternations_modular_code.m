if exist("V0_inner_LIPMSwingLeg.mat", 'file')
    load('V0_inner_LIPMSwingLeg.mat', 'V', 'model', 'S');
    disp('V and S for 0-step capturability are loaded');
    
    x = msspoly('x', model.num_states);
    t = msspoly('t', 1);
    u = msspoly('u', model.num_inputs);
    
    [f,g] = model.controlAffineDynamics(t,x);
    r = model.reset(t, x, []);
else
    g = 10;
    z_nom = 1;
    step_time = 0.3;
    cop_max = 0.05;
    u_wrt_cm = 0;
    x_leg = 0;

    model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm); 
    
    x = msspoly('x',model.num_states);
    t = msspoly('t',1);
    u = msspoly('u',model.num_inputs);
    system_dynamics = model.dynamics(t,x,u);
    [f,g] = model.controlAffineDynamics(t,x);
    r = model.reset(t, x, []);

    A = double(subs(diff(system_dynamics,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
    B = double(subs(diff(system_dynamics,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

    Q = eye(model.num_states);
    R = eye(model.num_inputs);
    [K,Q] = lqr(A,B,Q,R);
    
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
    
    figure(1);
    for i=1:30
      [ V,Bu ] = lipmSwingLegLyapunovAlternations(x,f,g,V,Bu,A_state);
      hold on
      if mod(i,2) == 0
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; x_leg],1,{'r'});
      else
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; x_leg],1,{'b'});
      end
    end

    S = Bu;
    save('V0_inner_LIPMSwingLeg.mat','V','model','S');
    return
end

% Load boundary function from external SOS program
% load('boundary_test.mat', 'B_plot');
% B = B_plot;
load('1step_boundary_test.mat', 'W');
B = 1 - W;

% Degree of Lyapunov function
degree = 2;

nX = length(x);
nU = length(u);

% Equation: V = x'*Q_init*x
Q_init = double(subs(diff(diff(V, x)', x)/2, x, zeros(nX, 1))); %#ok<UDIM>

% Switching controller possibilities
ndgrid_arg = mat2cell(repmat([-1;1],1,nU), 2, ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

% TODO: Scaling the function
% y = T*x
% => x = T^{-1}*y
% y = msspoly('y', nX);
% T = Q_init^(0.5);
% x_y = inv(T)*y; %#ok<MINV>
% r_y = subs(r, x, x_y);
% f_y = T*subs(f, x, x_y);
% g_y = T*subs(g, x, x_y);
% V_y = subs(V, x, x_y);
% B_y = subs(B, x, x_y);
% S_y = subs(S, x, x_y);
% Q_init_y = eye(nX);

outer_radius = 1;
inner_radius = 0.1;

figure(1);
subplot(2, 2, 1);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 1, {'b'});
hold on;
contourSpotless(B, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0, {'g'});

subplot(2, 2, 2);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'b'});
hold on;
contourSpotless(B, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 0, {'g'});

subplot(2, 2, 3);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 1, {'b'});
hold on;
contourSpotless(B, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 0, {'g'});

subplot(2, 2, 4);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 1, {'b'});
hold on;
contourSpotless(B, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 0, {'g'});

for j=1:50
    if j > 5
        B = find_B(t, x, u, f, g, r, V, S);
    end
    
    [mult, bmult, smult, mult_r, bmult_r] = ...
        binarySearchMultipliers(x, f, g, r, umat, V, S, B, ...
        outer_radius, inner_radius, eye(nX), 1);
%     [V, S, B] = ...
%         binarySearchVSandB(x, f, g, r, umat, V, Q_init, ...
%         mult, bmult, smult, mult_r, bmult_r, ...
%         outer_radius, inner_radius, eye(nX), degree);
    [V, S] = ...
        binarySearchVandS(x, f, g, r, umat, V, Q_init, B, ...
        mult, bmult, smult, mult_r, bmult_r, ...
        outer_radius, inner_radius, eye(nX), degree);
    
%     [mult, bmult, smult, mult_r, bmult_r] = ...
%         binarySearchMultipliers(y, f_y, g_y, r_y, umat, V_y, S_y, B_y, ...
%         outer_radius, inner_radius, eye(nX), 1);
%     [V, S, B] = ...
%         binarySearchVSandB(y, f_y, g_y, r_y, umat, V_y, Q_init_y, ...
%         mult, bmult, smult, mult_r, bmult_r, ...
%         outer_radius, inner_radius, degree);
    
    if mod(i, 2) == 0
        figure(1);
        subplot(2, 2, 1);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0],1,{'b'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0], 0, {'g'});
        subplot(2, 2, 2);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.1],1,{'b'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.1], 0, {'g'});
        subplot(2, 2, 3);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.2],1,{'b'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.2], 0, {'g'});
        subplot(2, 2, 4);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.5],1,{'b'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.5], 0, {'g'});
    else
        figure(1);
        subplot(2, 2, 1);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0],1,{'r'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0], 0, {'g'});
        subplot(2, 2, 2);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.1],1,{'r'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.1], 0, {'g'});
        subplot(2, 2, 3);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.2],1,{'r'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.2], 0, {'g'});
        subplot(2, 2, 4);
        contourSpotless(V,x(1),x(2),[-1 1],[-1 1],[t; x(3)],[0; 0.5],1,{'r'});
        contourSpotless(B, x(1), x(2), [-1,1], [-1,1], [t; x(3)], [0; 0.5], 0, {'g'});
    end
end