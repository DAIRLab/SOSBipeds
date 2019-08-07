lipmSLA = LIPMSwingLegAutonomous(10, 1, 1, 0.1, 0);
[t, y] = ode45(@(t, x) dynamics(lipmSLA, t, x), [0, 1], [0.1; -0.1; 0.1]);
plot(t, y)
