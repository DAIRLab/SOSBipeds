classdef LIPMSwingLegAutonomous < NStepROASystem

  properties
    g; % gravitational acceleration
    z_nom; % nominal center of mass height
    step_max; % max step distance
    T; % step time
    cop_max; % max distance between foot and CoP
    u_wrt_com;
  end

  methods
    function obj = LIPMSwingLegAutonomous(g, z_nom, step_time, cop_max, u_wrt_com)
      if cop_max > 0
          num_inputs = 2;
      else
          num_inputs = 1;
      end
      obj@NStepROASystem(3, num_inputs, 0);
      obj.g = g;
      obj.z_nom = z_nom;
      obj.T = step_time;
      obj.cop_max = cop_max;
      assert(u_wrt_com == 1 || u_wrt_com == 0);
      obj.u_wrt_com = u_wrt_com;
    end

    function xdot = dynamics(obj, t, x)
      [f,g] = controlAffineDynamics(obj, t, x);
      xdot = f + g*controller(obj, t, x);
    end

    function u = controller(obj, t, x)
        u = zeros(obj.num_inputs, 1);
        omega_0 = sqrt(obj.g/obj.z_nom);
%         u(1) = (2 + 1/omega_0)*(x(1) + x(2)/omega_0);
%         u(2) = (2 + 1/omega_0)*(x(1) + x(2)/omega_0 - x(3));
        u(1) = sign(x(1) + x(2)/omega_0);
        u(2) = sign(x(1) + x(2)/omega_0 - x(3));
    end

    function [f,g] = controlAffineDynamics(obj,t,x)
      qcm = x(1);
      vcm = x(2);
      ql = x(3);
      f = [vcm; qcm*obj.g/obj.z_nom; obj.u_wrt_com*vcm];
      if obj.cop_max > 0
          g = [0, 0; -obj.cop_max*obj.g/obj.z_nom, 0; 0, 1];
      else
          g = [0; 0; 1];
      end
    end

    function xp = reset(obj, t, xm, s)
      % control input changes q only
      % qp = qm - u
      qcm_m = xm(1);
      vcm_m = xm(2);
      ql_m = xm(3);
      xp = [qcm_m - ql_m; vcm_m; -ql_m];
    end

    function ret = inputLimits(obj, u, x)
      ret = 1 - u'*u;
    end

    function[umin, umax, A] = simpleInputLimits(obj, x)
      umin = -1;
      umax = 1;
      A = [];
    end

    % Transform u
    % y = A*u + b s.t. limits on y are |y_i| <= 1
    % u = C*y + d inverse transform
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = 1;
      b = 0;
      C = 1;
      d = 0;
    end


    function ret = resetInputLimits(obj, s)
      ret = [];
    end

    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u_sol, sub_val)
      if nargin <= 9
          sub_val = [0; 0];
      end
      qcm = x(1);
      vcm = x(2);
      ql = x(3);

      sub_vars = [t; ql(1)];
      plot_vars = [qcm(1);vcm(1)];

      figure(1)
      contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)], sub_vars, sub_val, [1 0], {'b','r'});
      hold on
      grid on
      xlabel('q_1')
      ylabel('v_1')
      title('W(x)')

      % similar to Koolen et. al IJRR
      figure(n*10+2)
      omega = sqrt(obj.g/obj.z_nom);
      r_ic = qcm + vcm/omega;
      [umin, umax, ~] = obj.simpleInputLimits([]);
      if n == 1 && obj.cop_max == 0
        x = linspace(-R_diag(1), R_diag(1), 25);
        y = linspace(-R_diag(2), R_diag(2), 25);
        [X, Y] = meshgrid(x, y);
        scatter(reshape(X, 25*25, 1), reshape(Y, 25*25, 1), 'k');
        hold on;
        [X, Y] = lipmSwingLegCaptureLimit(x, y, sub_val, obj);
        scatter(X, Y, 'r');
        conservative_limits = (r_ic - sub_val(2)*exp(-omega*obj.T))'*(r_ic - sub_val(2)*exp(-omega*obj.T));
        dN = (umax*obj.T)*exp(-omega*obj.T);
      elseif n == 0
        step_max = umax*obj.T;
        conservative_limits = r_ic'*r_ic;
        dN = lipmCaptureLimit(obj.T, obj.cop_max, step_max, obj.z_nom, obj.g, n); % theoretical max ICP distance
      end

      if n==0
          contourSpotless([Vsol; h_X; conservative_limits],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0 dN^2],{'b','r','g'});
      else
          contourSpotless([Vsol; h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'b','r','g'});
      end
      hold on
      grid on
      xlabel('q_1')
      ylabel('v_1')
      title(['V(', num2str(sub_val(1)), ',x) with x_3 = ', num2str(sub_val(2))])

      figure(n*10+4)
      hold off
      h=contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'k','r'});
      hold on
      grid on
      h=contourSpotless(u_sol,plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val);
      xlabel('q_1')
      ylabel('v_1')
      title('u(0,x)')
    end

    function draw(obj,t,x)
      x_stance = x(end-2:end-1);
      % draw line from origin to COM
      h=line([x_stance(1);x(1)],[x_stance(2);0]);
      set(h,'LineWidth',3,'Color','red')
      radius = .1;
      rectangle('Position',[[x(1);0]+x_stance-radius/2;radius;radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-2 2])
      ylim([-2 2])
      axis off
    end

  end
end
