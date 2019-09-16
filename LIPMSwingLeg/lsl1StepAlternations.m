% Load 0-step capturability solution
load('V0_inner_LIPMSwingLeg.mat', 'V', 'S', 'model');
load('W1_inner.mat', 'W');

% Load 1-step capturability solution
% load('V1_inner_safe.mat', 'V', 'S', 'model', 'W');

x = msspoly('x', model.num_states);
t = msspoly('t', 1);
u = msspoly('u', model.num_inputs);

[f, g] = model.controlAffineDynamics(t, x);
r = model.reset(t, x, []);

V_r = subs(V, x, r);

display_array = [0, 0, 0, 0; 0, 0.1, 0.2, 0.3];
% Display for W <= 1
plot_figures(1, [V, W], [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1, 1], ['r', 'g']);
load('V1_LIPMSwingLeg.mat', 'Vsol');
plot_figures(1, Vsol, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [0], ['m']);
plot_figures(1, x'*[1, 0, 0; 0, 1, 0; 0, 0, 4]*x, [x(1), x(2)], [-1, 1], [-1, 1], ...
    [t; x(3)], display_array, [1], ['k']);

A_state = {};
A_state{1} = diag([0, 0, 0]);
V0 = V;
S0 = S;
S0(2) = x(1) - x(3) +  x(2)/sqrt(model.g);
[ V,S ] = lsl1StepLyapAlternations(x,f,g,V0,S0,W);

contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t;x(3)],[0; 0.01],1,{'r'});

for i=1:60
    S(2) = x(1) + (x(2)/sqrt(model.g)) - x(3);
    [ V,S ] = lsl1StepLyapAlternations(x,f,g,V,S,W);

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
    sqrt(1./diag(double(subs(diff(diff(subs(V,t,0),x)',x)/2,x,x*0))))
end

keyboard;