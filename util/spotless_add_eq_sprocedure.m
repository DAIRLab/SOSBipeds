function [prog, eqn, mult, coefmult ] = spotless_add_eq_sprocedure(prog, eqn, g, vars, degree)
%SPOTLESS_ADD_EQ_SPROCEDURE adds S-procedure multiplier for an equality
% constraint
% @param prog the spot program 
% @param eqn the current SOS equation
% @param g the equality constraint, g(vars) = 0
% @param vars the indeterminates to use in the multiplier
% @param degree the degree of the multiplier [0,degree]
%
% @return prog the modified program
% @return eqn the modified SOS equation
% @param mult the multiplier polynomial
% @param coefmult the coefficientgs of mult, as decision variables 

[prog,mult,coefmult] = prog.newFreePoly(monomials(vars,0:degree));
eqn = eqn - g*mult;

% display(sprintf('S-proc eq. SOS deg: %d, g deg: %d, mult deg: %d',full(deg(eqn,vars)), deg(g), degree))

end

