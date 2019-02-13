g = 10;
z_nom = 1;
step_max = .7;
step_time = 0.3;
cop_max = .1; % set to 0 to get point foot model with no continuous inputs

model = LIPM2D(g, z_nom, step_max, step_time, cop_max);


% Get an initial quadratic Lyapunov candidate
x = msspoly('x',model.num_states);
t = msspoly('t',1);
u = msspoly('u',model.num_inputs);
f = model.dynamics(t,x,u);

A = double(subs(diff(f,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
B = double(subs(diff(f,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

Q = 100*eye(model.num_states);
R = eye(model.num_inputs);
[K,Q] = lqr(A,B,Q,R);

%
V0 = x'*Q*x;
[V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V0*200)
figure(1)
hold off
  contourSpotless(V,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),1);

%% 0 step (balancing) alternations
for i=1:40,
  [V,u_fn] = quadraticControlLyapunovAlternations(x,u,f,V);
%   figure(1)
  hold off
  contourSpotless(V,x(1),x(2),[-1 1],[-3 3],[t;x(3:end)],zeros(model.num_states-1,1),1);
end;

%%
V_0step = V;
s = msspoly('s',model.num_reset_inputs);
r = model.reset(t,x,s);

% get inner approximation of pre-reset
% s = (c'SA)^-1 * c'Sx
c=double(diff(r,s));
Q_V0=double(diff(diff(V,x)',x))/2;

a = c'*Q_V0*c;
b = 2*x'*Q_V0*c;
d = x'*Q_V0*x - 1;
%% init
rho2 = msspoly(1) + 20*t;
% V2 = x'*Q*x*500;
V2 = V_0step*20;
% [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u,f,V2,step_time,rho2,a,b,d);


%% 1 step alternations
for i=1:25,
  [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u,f,V2,step_time,rho2,a,b,d,[],false);
% [V2,u2,rho2] =  quadraticControlAlternationsWithResetNoGThreeSteps(x,u_alt,f_alt,V2,step_time,rho2,a,b,d,constraint_alt);
  figure(2)

  % green: 0 step viable capture region
  % black: 1 step viable capture region (t=0)
  % red: 1-step states for (t=T)
  % blue: max-boundary on the red (pre-reset states that map into the
  % green)
  
  hold off
  contourSpotless(V,x(1),x(2),[-2 2],[-3 3],t,zeros(model.num_states-1,1),1,{'g'});
  hold on
  contourSpotless(V2,x(1),x(2),[-2 2],[-3 3],t,zeros(model.num_states-1,1),dmsubs(rho2,t,0),{'k'});
  contourSpotless(V2,x(1),x(2),[-2 2],[-3 3],t,[step_time;zeros(model.num_states-2,1)],dmsubs(rho2,t,step_time),{'r'});
  contourSpotless(b^2-4*a*d,x(1),x(2),[-2 2],[-3 3],t,[step_time;zeros(model.num_states-2,1)],0,{'y'});
  contourSpotless(b,x(1),x(2),[-2 2],[-3 3],t,[step_time;zeros(model.num_states-2,1)],0,{'b'});
  contourSpotless(2*a-b,x(1),x(2),[-2 2],[-3 3],t,[step_time;zeros(model.num_states-2,1)],0,{'b'});

  hold off
  
end
%%