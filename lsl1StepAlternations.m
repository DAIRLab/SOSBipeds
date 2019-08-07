% load('V0_inner_LIPMSwingLeg.mat', 'V', 'S', 'model');
% load('W1_inner.mat', 'W');
load('V1_inner_safe.mat', 'V', 'S', 'model', 'W');

x = msspoly('x', model.num_states);
t = msspoly('t', 1);
u = msspoly('u', model.num_inputs);

[f, g] = model.controlAffineDynamics(t, x);
r = model.reset(t, x, []);

V_r = subs(V, x, r);

figure(1);
subplot(2, 2, 1);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 1, {'b'});
hold on;
contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 1, {'g'});
load('V1_LIPMSwingLeg.mat', 'Vsol')
contourSpotless(Vsol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0, {'m'});
contourSpotless(1 - x'*x, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 0, {'k'});

subplot(2, 2, 2);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'b'});
hold on;
contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'g'});
contourSpotless(Vsol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 0, {'m'});
contourSpotless(1 - x'*x, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 0, {'k'});

subplot(2, 2, 3);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 1, {'b'});
hold on;
contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 1, {'g'});
contourSpotless(Vsol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 0, {'m'});
contourSpotless(1 - x'*x, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 0, {'k'});

subplot(2, 2, 4);
contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 1, {'b'});
hold on;
contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 1, {'g'});
contourSpotless(Vsol, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 0, {'m'});
contourSpotless(1 - x'*x, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 0, {'k'});

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

    %   figure(1)
    if mod(i, 2) == 0
        figure(1);
        subplot(2, 2, 1);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0],1,{'b'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 1, {'g'});
        subplot(2, 2, 2);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.1],1,{'b'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'g'});
        subplot(2, 2, 3);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.2],1,{'b'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 1, {'g'});
        subplot(2, 2, 4);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.5],1,{'b'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 1, {'g'});
    else
        figure(1);
        subplot(2, 2, 1);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0],1,{'r'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0], 1, {'g'});
        subplot(2, 2, 2);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.1],1,{'r'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.1], 1, {'g'});
        subplot(2, 2, 3);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.2],1,{'r'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.2], 1, {'g'});
        subplot(2, 2, 4);
        contourSpotless(V,x(1),x(2),[-1, 1],[-1, 1],[t; x(3)],[0; 0.5],1,{'r'});
        contourSpotless(W, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.5], 1, {'g'});
    end
    sqrt(1./diag(double(subs(diff(diff(subs(V,t,0),x)',x)/2,x,x*0))))
end

keyboard;