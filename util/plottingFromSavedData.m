function [] = plottingFromSavedData(obj, n, Vsol, h_X, R_diag, t, x)
    qcm = x(1);
    vcm = x(2);
    ql = x(3);

    sub_vars = [t; ql(1)];
    sub_val = [0; 0];
    plot_vars = [qcm(1);vcm(1)];

    figure(1)
    omega = sqrt(obj.g/obj.z_nom);
    r_ic = qcm + vcm/omega;
    cop_max = 0;
    [umin, umax, ~] = obj.simpleInputLimits([]);
    if n == 1
      x = linspace(-R_diag(1), R_diag(1), 25);
      y = linspace(-R_diag(2), R_diag(2), 25);
      options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
      for i = x
          for j = y
              icp = i + j/omega;
              if icp*exp(omega*obj.T) - sub_val(2) > 0 && icp - sub_val(2) > 0
                  [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - sub_val(2) - umax*t, obj.T, options);
              elseif icp*exp(omega*obj.T) - sub_val(2) <= 0 && icp - sub_val(2) <= 0
                  [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - sub_val(2) - umin*t, obj.T, options);
              else
                  [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - sub_val(2) - umax*t, obj.T, options);
                  if exitflag < 0 || t < 0 || t > obj.T
                      [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - sub_val(2) - umin*t, obj.T, options);
                  end
              end
              if exitflag > 0 && t > 0 && t < obj.T
                  scatter(i, j, 'r');
                  hold on;
              else
                  scatter(i, j, 'k');
                  hold on;
              end
          end
      end
      dN = (umax*obj.T)*exp(-omega*obj.T);
    else
      step_max = umax*obj.T;
     dN = lipmCaptureLimit(obj.T, cop_max, step_max, obj.z_nom, obj.g, n); % theoretical max ICP distance
    end

    contourSpotless([Vsol; h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(2) R_diag(2)],sub_vars,sub_val,[0 0],{'b','r'});
    xlabel('q_1')
    ylabel('v_1')
    title('V(0,x)')
end
