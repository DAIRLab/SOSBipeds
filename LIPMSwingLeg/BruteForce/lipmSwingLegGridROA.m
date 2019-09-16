function [] = lipmSwingLegGridROA(inner_approx)
if nargin < 1
    inner_approx = 0;
end
if inner_approx==2
    if isfile('gridded_ROA.mat')
        load('gridded_ROA.mat', 'X1', 'X2', 'X3')
        scatter3(X1, X2, X3, 'b');
        hold on;
    end
    if isfile('gridded_1step_goal.mat')
        load('gridded_1step_goal.mat', 'W1', 'W2', 'W3')
        scatter3(W1, W2, W3, 'g');
    end
    if isfile('all_points_inside_V.mat')
        load('all_points_inside_V.mat', 'X1', 'X2', 'X3')
        scatter3(X1, X2, X3, 'r');
    end
    if isfile('points_inside_V.mat')
        load('points_inside_V.mat', 'X1', 'X2', 'X3')
        scatter3(X1, X2, X3, 'r');
    end
    xlabel('$x_{cm}$', 'Interpreter', 'latex');
    ylabel('$v_{cm}$', 'Interpreter', 'latex');
    zlabel('$x_{leg}$', 'Interpreter', 'latex');
    legend('1-step capturability', '1-step goal region', 'Inner Approximation');
    return
end
load('V0_inner_LIPMSwingLeg.mat', 'V');
V0 = V;
load('V1_inner_swing_leg.mat', 'V', 'S', 'W', 'model');
x = msspoly('x', model.num_states);
S(1) = x(1) + x(2)/sqrt(model.g/model.z_nom);
S(2) = x(1) - x(3) +  x(2)/sqrt(model.g/model.z_nom);

X1 = []; X2 = []; X3 = [];
W1 = []; W2 = []; W3 = [];
for i = -1:0.1:1
    for j = -1:0.1:1
        for k = -1:0.1:1
            x_s = [i; j; k]
            if double(subs(W, x, x_s)) <= 1
                W1 = [W1; i];
                W2 = [W2; j];
                W3 = [W3; k];
            end
            if inner_approx && double(subs(V, x, x_s)) > 1
                continue
            end
            success = false;
            for step = 1:1000
                u = double(subs(S, x, x_s));
                u = (u > 0)*2 - 1;
                dx = model.dynamics([], x_s, u);
                x_s = x_s + dx*0.01;
                if(x_s'*x_s > 1)
                    break
                end
                if double(subs(W, x, x_s)) <= 1 || double(subs(V, x, x_s)) <= 1
                    success = true;
                    break
                end
            end
            if success == true
                X1 = [X1; i];
                X2 = [X2; j];
                X3 = [X3; k];
            end
        end
    end
end
if inner_approx == 0
    scatter3(X1, X2, X3, 'b')
elseif inner_approx == 1
    scatter3(X1, X2, X3, 'r')
end
hold on;
scatter3(W1, W2, W3, 'g')
keyboard
end