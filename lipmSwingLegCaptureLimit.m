function [X, Y] = lipmSwingLegCaptureLimit(qcm, vcm, sub_val, obj)
    t_lim = sub_val(1);
    ql = sub_val(2);
    omega = sqrt(obj.g/obj.z_nom);
    [umin, umax, ~] = obj.simpleInputLimits([]);
    X = [];
    Y = [];
    options = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    for i = qcm
      for j = vcm
        icp = i + j/omega;
%         if icp*exp(omega*obj.T) - ql > 0 && icp - ql > 0
%           [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umax*t + obj.u_wrt_com*(qcm - qcm*exp(obj.g*t/obj.z_nom)), obj.T, options);
%         elseif icp*exp(omega*obj.T) - ql <= 0 && icp - ql <= 0
%           [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umin*t + obj.u_wrt_com*(qcm - qcm*exp(obj.g*t/obj.z_nom)), obj.T, options);
%         else
%           [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umax*t + obj.u_wrt_com*(qcm - qcm*exp(obj.g*t/obj.z_nom)), obj.T, options);
%           if exitflag < 0 || t < 0 || t > obj.T
%               [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umin*t + obj.u_wrt_com*(qcm - qcm*exp(obj.g*t/obj.z_nom)), obj.T, options);
%           end
%         end
        c1 = 0.5*(i + j/omega);
        c2 = 0.5*(i - j/omega);
        [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umax*t - obj.u_wrt_com*(c1*exp(omega*t) + c2*exp(-omega*t) - c1 - c2), obj.T - t_lim, options);
        if exitflag <= 0 || t < 0 || t >= obj.T - t_lim
            [t, ~, exitflag, ~] = fsolve(@(t) icp*exp(omega*t) - ql - umin*t - obj.u_wrt_com*(c1*exp(omega*t) + c2*exp(-omega*t) - c1 - c2), obj.T - t_lim, options);
        end
        if exitflag > 0 && t > 0 && t <= obj.T 
            X = [X, i];
            Y = [Y, j];
        end
      end
    end
end