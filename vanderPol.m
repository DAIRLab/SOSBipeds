function xdot = vanderPol(mu, t, x)
   xdot = [0; 0];
%    xdot(1) = x(2);
%    xdot(2) = mu*(1 - x(1)^2)*x(2) - x(1);
   xdot(1) = 2*x(2);
   xdot(2) = - 0.8*x(1) + 2.1*x(2) - 10*x(2)*x(1)^2;
end
