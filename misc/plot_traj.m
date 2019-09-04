function [] = plot_traj(x)
load('V0_inner_LIPMSwingLeg.mat', 'V');
V0 = V;
load('V1_inner_swing_leg.mat', 'V', 'S', 'model', 'W');
x_mss = msspoly('x', model.num_states);
S(1) = x_mss(1) + x_mss(2)/sqrt(model.g/model.z_nom);
S(2) = x_mss(1) - x_mss(3) +  x_mss(2)/sqrt(model.g/model.z_nom);

t = msspoly('t', 1);
x_s = msspoly('x', model.num_states);

dt = 0.01;
x_hist = x;
has_stepped = false;
for i=1:1000
    u = double(subs(S, x_s, x));
    u = (u > 0)*2 - 1;
    dx = model.dynamics(t, x, u);
    x = x + dx*dt;
    x_hist = [x_hist, x]; %#ok<AGROW>
    if double(subs(W, x_s, x)) < 1 && ~has_stepped
        has_stepped = true;
        fprintf('Step\n');
        fprintf('%d\n%d\n%d\n',x(1), x(2), x(3));
        disp(double(subs(W, msspoly('x', 3), x)));
        br = x;
        x = model.reset([], x, []);
        ar = x;
        x_hist = [x_hist, ar];
        load('V0_inner_LIPMSwingLeg', 'V', 'S', 'model');
    elseif has_stepped && abs(x(1)) < 0.01 && abs(x(2)) < 0.01
        fprintf('Balance');
        i
        break
    end
    if x(1)^2 + x(2)^2 + x(3)^2 > 1
        fprintf("Outside the ball\n");
    end
end

figure(1);
plot3(x_hist(1, :), x_hist(2, :), x_hist(3, :), 'm', 'LineWidth', 3);
hold on;
scatter3(br(1), br(2), br(3), 'g');
scatter3(ar(1), ar(2), ar(3), 'r');
grid on;
xlabel('x');
ylabel('x_dot');
zlabel('x_leg');
end