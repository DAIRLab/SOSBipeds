g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = 0.05;
u_wrt_cm = 0;
x_leg = 0;

model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm);

x = msspoly('x', 3);
t = msspoly('t', 1);
u = msspoly('u', 2);
r = model.reset(t, x, []);

load('V0_inner_LIPMSwingLeg.mat', 'V');
V_r = subs(V, x, r);

X = []; Y = []; Z = [];
X_r = []; Y_r = []; Z_r = [];

for i=linspace(-1, 1, 100)
    for j=linspace(-1, 1, 100)
        for k=linspace(-1, 1, 100)
            v_subs = subs(V, x, [i; j; k]);
            if double(v_subs) < 1
                vr_subs = subs(V_r, x, [i; j; k]);
                if double(vr_subs) < 1
                    r = model.reset([], [i; j; k], []);
                    X = [X, i]; Y = [Y, j]; Z = [Z, k];
                    X_r = [X_r, r(1)]; Y_r = [Y_r, r(2)]; Z_r = [Z_r, r(3)];
                end
            end
        end
    end
end

figure;
scatter3(X, Y, Z);