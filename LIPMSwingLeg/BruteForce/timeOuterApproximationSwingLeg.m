function timeOuterApproximationSwingLeg(n)
% Constants
g = 10;
z_nom = 1;
step_time = 0.3;
cop_max = 0.05;
u_wrt_cm = 1;

% Model of the robot
model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm);

% Constant variable
x_leg = 0.5;

% Load or generate Lyapunov function for the current step
filename = solutionFileName(model, n);
if ~exist(filename, 'file')
  fprintf('Solving for %d-step first\n',n)
  lipmSwingLegNStepCapturability(n)
end
data = load(filename);
V1 = data.Vsol;

% Load or generate data for previous step
filename = solutionFileName(model, n - 1);
data = load(filename);
V0 = data.Vsol;

x = msspoly('x', 3);
t = msspoly('t', 1);

% Grid the space and look for times where stepping will be successful
x1 = linspace(-1, 1, 15);
x2 = linspace(-1, 1, 15);
T = cell(length(x1), length(x2));
for i=1:length(x1)
    for j=1:length(x2)
        T{i,j} = [];
    end
end

for i=1:length(x1)
    for j=1:length(x2)
        if norm([x1(i); x2(j); x_leg]) > 1
            continue
        end
%         xp = model.reset(t, [x1(i); x2(j); x_leg], []);
%         V0_i = subs(V0, [t; x], [0; xp]);

%         if double(V0_i) >= 0
        V1_i = subs(V1, x, [x1(i); x2(j); x_leg]);
        for k=linspace(0, step_time, 20)
            if double(subs(V1_i, t, k)) >= 0
                T{i,j} = [T{i,j}, k];
            end
        end
%         end
    end
end

% Plotting
for i=1:length(x1)
    for j=1:length(x2)
        for k=1:length(T{i,j})
            scatter3(x1(i), x2(j), T{i,j}(k), 'b');
            hold on;
        end
    end
end

h0 = contourSpotless(V1, x(1), x(2), [-1 1], [-1 1], [x(3); t], [x_leg; 0], [0 0], {'r'});
hT = contourSpotless(V1, x(1), x(2), [-1 1], [-1 1], [x(3); t], [x_leg; step_time], [0 0], {'g'});
hT{1}.ContourZLevel = step_time;
xlim([-1, 1])
ylim([-1 1])
xlabel('x_{cm}')
ylabel('v_{cm}')
zlabel('Time')
end

function filename = solutionFileName(model, n)
    filename_suffix = class(model);
    filename = sprintf(['V%d_' filename_suffix '.mat'], n);
end