function weights = spotlessIntegral(prog,poly,sphere_vars,A_diag,box_vars,box_lims)
% weights = spotlessIntegral(poly,sphere_vars,A,box_vars,box_lims)
%  Computes the integral of polynomial over a given region
%    If the coefficients of poly are `coeff`, and then the result is
%  
%    <weights, coeff> = int(poly)
%
%   @param sphere_vars the indeterminates corresponding to an ellipsoidal
%   region
%   @param box_vars the indetermintes corresponding to a box region
%
%   Where the region is specified in two ways: variables either correspond
%   to an ellipsoid with diagonal matrix
%            {sphere_vars'*A_diag*sphere_vars <= 1}
%   or a box regbion
%          box_lims(i,1) <= box_vars(i) <= box_lims(i,2)
% 

[vars,alphas,coeff] = decomp(poly, [prog.freeVar;prog.coneVar]);

% find sphere_vars
n_sphere = length(sphere_vars);

sphere_inds = [];
sphere_var_inds = [];
sphere_var_inds_c = [];
for i=1:length(sphere_vars)
  is_var = false;
  for j=1:length(vars)
    if isequal(vars(j),sphere_vars(i))
      sphere_inds(end+1) = j;
      sphere_var_inds(end+1) = i;
      is_var = true;
      break;
    end
  end
  if ~is_var
    sphere_var_inds_c(end+1) = i;
  end
end

% find box
n_box = length(box_vars);

box_inds = [];
box_var_inds = [];
box_var_inds_c = [];
for i=1:length(box_vars)
  is_var = false;
  for j=1:length(vars)
    if isequal(vars(j),box_vars(i))
      box_inds(end+1) = j;
      box_var_inds(end+1) = i;
      is_var = true;
      break;
    end
  end
  if ~is_var
    box_var_inds_c(end+1) = i;
  end
end

assert(length([sphere_inds box_inds]) == length(unique([sphere_inds box_inds])));
assert(length([sphere_inds box_inds]) == length(vars))

sphere_alphas = zeros(size(alphas,1),length(sphere_vars));
sphere_alphas(:,sphere_var_inds) = alphas(:,sphere_inds);

box_alphas = zeros(size(alphas,1),length(box_vars));
box_alphas(:,box_var_inds) = alphas(:,box_inds);

% eliminate box variables first
% /int x^alpha = 1/(alpha+1)*(x_max^(alpha+1) - x_min^(alpha+1))
for i=1:length(box_vars)
  pows = box_alphas(:,i) + 1;
  vals = 1./pows .* (box_lims(i,2).^pows - box_lims(i,1).^pows);
  coeff = coeff.*vals';
end

if n_sphere > 0
  if isempty(alphas)
    sphere_alphas = zeros(size(A_diag));
  end
  
  sphere_betas = 0.5*(sphere_alphas + 1);
  Ra = (1.^(sum(sphere_alphas,2) + n_sphere))./(sum(sphere_alphas,2) + n_sphere);
  IS = 2*prod(gamma(sphere_betas),2)./(gamma(sum(sphere_betas,2)));
  l = Ra.*IS;
  alphaszero = (mod(sphere_alphas,2) ~= 0);
  alphaszero = any(alphaszero,2);
  l(alphaszero) = 0;
  
  l = l.*prod(repmat(A_diag,size(sphere_alphas,1),1).^(sphere_alphas+1),2);
  
  weights = coeff*l;
else
  weights = sum(coeff);
end
end