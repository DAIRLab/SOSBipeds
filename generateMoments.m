function [moments] = generateMoments(prog, nX, deg, shape)

% maybe we can extract x from prog.
x = msspoly('x', nX);

% N: the number of coefficients
N = nchoosek(nX + deg, deg);
moments = zeros(N, 1);
monomialList = monomials(x, 0:deg);

sphere_vars = shape.sphere_vars;
A_diag = shape.A_diag;
box_vars = shape.box_vars;
box_lims = shape.box_lims;

for i = 1: N 
	weight = spotlessIntegral(prog, monomialList(i), x, A_diag, box_vars, box_lims);
	monomialList(i) = weight;
end

end