function [prog, eqn, mult, coefmult] = spotless_add_sprocedure(prog, eqn, h, vars, min_degree, degree)
%SPOTLESS_ADD_SPROCEDURE adds S-procedure multiplier for an inequality
% constraint. Uses:
% [prog, eqn, mult, coefmult] = spotless_add_sprocedure(prog, eqn, h, vars)
% [prog, eqn, mult, coefmult] = spotless_add_sprocedure(prog, eqn, h, vars, min_degree)
% [prog, eqn, mult, coefmult] = spotless_add_sprocedure(prog, eqn, h, vars, min_degree, degree)
%
% @param prog the spot program
% @param eqn the current SOS equation
% @param h the equality constraint(s), g(vars) = 0
%    can be specified as a scalar constraint or as a vector
% @param vars the indeterminates to use in the multiplier
% @param min_degree the minimum degree to use for polynomials (default 0)
% @param degree the degree of the multiplier(s). If not provided, will
%   automatically determine the degree based on the total degree of eqn
%   and of the mulitpliers.
%   Can be specified as a scalar, or as a vector (one per-constraint)
% @return prog the modified program
% @return eqn the modified SOS equation
% @param mult the multiplier polynomial (must be SOS)
% @param coefmult the coefficientgs of mult, as decision variables

% eqn_deg = full(deg(eqn,vars));
% if ~even(eqn_deg)
%   eqn_deg = eqn_deg + 1;
% end
% degree = max(2,min(degree, 2*floor((eqn_deg - deg(h))/2)));



if nargin < 5
  min_degree = 0;
end

if nargin < 6
  degree = [];
end

original_deg = even_degree(eqn,vars);


if isempty(degree)
  degree = zeros(length(h),1);
  for i = 1:length(h)
    degree(i) = original_deg - even_degree(h(i),vars);
  end
else
  if length(degree) == 1
    degree = repmat(degree,length(h),1);
  end  
end

mult = msspoly;
coefmult = msspoly;
for i = 1 : length(h)
  [prog,mult_i,coefmult_i] = prog.newSOSPoly(monomials(vars,min_degree:degree(i)));
  
  eqn = eqn - h(i) * mult_i;
  mult = [mult; mult_i]; %#ok<AGROW>
  coefmult = [coefmult; coefmult_i]; %#ok<AGROW>
  if isnumeric(h(i))
    deg_h = 0;
  else
    deg_h = deg(h(i));
  end
  if original_deg ~= even_degree(h(i),vars) + degree(i);
    warning('S-procedure degree mismatch')
  end
end
end

