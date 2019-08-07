load('V0_inner_LIPMSwingLeg.mat', 'V', 'S', 'model');

x = msspoly('x', model.num_states);
t = msspoly('t', 1);
u = msspoly('u', model.num_inputs);

r = model.reset(t, x, []);

nU = model.num_inputs;
ndgrid_arg = mat2cell(repmat([-1;1],1,nU),2,ones(1,nU)');
[ugrid{1:nU}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nU,nU);
for i=1:nU
  umat(:,i) = ugrid{i}(:);
end

% contourSpotless3D(V, x, 1, [1; 1; 1])
% hold on;
% contourSpotless3D(subs(V, x, r), x, 1, [1; 1; 1]);

contourSpotless(V, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.02], 1, {'b'});
hold on;
contourSpotless(subs(V, x, r), x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.02], 1, {'r'});

% [f, g] = model.controlAffineDynamics(t, x);
 
% colors = ['r', 'g', 'k', 'm'];
% for i = 1:2^nU
%     V_dot = diff(V, x)*(f + g*umat(i, :)');
%     contourSpotless(V_dot, x(1), x(2), [-1, 1], [-1, 1], [t; x(3)], [0; 0.02], 0, {colors(i)});
% end