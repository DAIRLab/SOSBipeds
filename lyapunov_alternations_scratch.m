function [V, B] = lyapunov_alternations_scratch(x, f, g, V0, S0, R)
if nargin < 6
    R = [];
end

nx = length(x);
nu = size(g, 2);

P_init = double(subs(diff(diff(V0, x)', x), x, zeros(nx, 1)));

ndgrid_arg = mat2cell(repmat([-1;1],1,nu),2,ones(1,nu)');
[ugrid{1:nu}] = ndgrid(ndgrid_arg{:});
umat = zeros(2^nu,nu);
for i=1:nu
  umat(:,i) = ugrid{i}(:);
end

y = msspoly('y',nX);
T = Q_init^(.5);
x_y = inv(T)*y;
f_y = T*subs(f,x,x_y);
g_y = T*subs(g,x,x_y);
V_y = subs(V,x,x_y);
B_y = subs(B,x,x_y);
R_y = cell(length(R),1);
for i=1:length(R),
  R_y{i} = inv(T)'*R{i}*inv(T);
end
Q_init_y = eye(nX);