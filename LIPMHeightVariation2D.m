classdef LIPMHeightVariation2D < NStepCapturabilitySOSSystem
  % 2D model with variable height and a foot
  
  properties
    gravity; % gravitational acceleration
    z_nom; % nominal height
    step_max; % max step distance
    T; % step time
    foot_radius
    uz_bnd
    z_inv_degree
  end
  
  methods
    function obj = LIPMHeightVariation2D(g, z_nom, step_max, step_time, uz_bnd, foot_radius,z_inv_degree)
      %  obj = LIPMHeightVariation2D(g, z_nom, step_max, step_time, uz_bnd, foot_radius,z_inv_degree)
      %   @param g gravity magnitude
      %   @param z_nom nominal height (m)
      %   @param step_max max step length (m)
      %   @param step_time stepping period (s)
      %   @param uz_bnd the limit on vertical acceleration (in g's)
      %   @param foot_radius (m)
      %   @param z_inv_degree the degree of the Taylor expansion of 1/z
      obj@NStepCapturabilitySOSSystem(4, 2, 1);
      obj.gravity = g;
      obj.z_nom = z_nom;
      obj.step_max = step_max;
      obj.T = step_time;      obj.foot_radius = foot_radius;
      obj.uz_bnd = uz_bnd;
      obj.z_inv_degree = z_inv_degree;
    end
    
    % x = [q;v]
    % q = [x_com;z_com;theta]
    function xdot = dynamics(obj, t, x, u)
      [f,g] = controlAffineDynamics(obj, t, x);
      xdot = f + g*u;
    end
    
    function [f,g] = controlAffineDynamics(obj,t,x)
      q = x(1 : 2);
      q(2) = q(2) + obj.z_nom;
      v = x(3 : 4);
      
      z_inv = obj.zInverse(q(2));
      
      f_x_ddot = obj.gravity*z_inv*q(1);
      f_z_ddot = 0;      
      g_x_ddot = obj.gravity*z_inv*[-1 q(1)];
      g_z_ddot = [0 obj.gravity];
      
      f = [v;f_x_ddot;f_z_ddot];
      g = [zeros(2); g_x_ddot; g_z_ddot];
      g = g*diag([obj.foot_radius;obj.uz_bnd]);
    end
    
    % taylor approximation of 1/z
    function h = zInverse(obj,z)
      h = 1/obj.z_nom;
      for i=1:obj.z_inv_degree,
        h = h + (-1)^i*obj.z_nom^(-i-1)*(z-obj.z_nom)^i;
      end
    end
    
    function xp = reset(obj, t, xm, s)
      % control input changes q(1) only
      % qp = qm - u
      qm = xm(1 : 2);
      vm = xm(3 : 4);
      xp = [qm - [obj.step_max; 0] * s; vm];
    end
    
    function ret = inputLimits(obj, u, x)
      ret = 1 - u.*u;
    end
    
    function[umin,umax,A] = simpleInputLimits(obj,x)
      umin = -ones(2,1);
      umax = ones(2,1);
      A = [];
    end
    
    function [A,b,C,d] = unitBoxInputTransform(obj)
      A = 1;
      b = 0;
      C = A;
      d = b;
    end       
    
    function ret = resetInputLimits(obj, s)
      ret = 1 - s'*s;
    end

    function [smin,smax] = simpleResetInputLimits(obj,x)
      smin = -1;
      smax = 1;
    end

    function rp = stanceUpdate(obj,x,r,s)
      rp = r + obj.step_max*[s;0];
    end
    
    function plotfun(obj, n, Vsol, Wsol, h_X, R_diag, t, x, u)
      q = x(1:2);
      v = x(3:4);
      
      sub_vars = [q(2);v(2);t];
      sub_val = [0;0;0];
      plot_vars = [q(1);v(1)];
      
      if ~isempty(Wsol)
        figure(1)
        contourSpotless([Wsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[1 0],{'b','r'});
        xlabel('q_1')
        ylabel('v_1')
        title('W(x)')
      end
      
      figure(n*10+2)
      contourSpotless([Vsol;h_X],plot_vars(1),plot_vars(2),[-R_diag(1) R_diag(1)],[-R_diag(3) R_diag(3)],sub_vars,sub_val,[0 0],{'b','r'});
      xlabel('q_1')
      ylabel('v_1')
      title('V(0,x)')
      
%       % 3d plot for t = 0, zdot = 0
%       contour_inds = [1 3 2];
%       sub_inds = setdiff(1:6,contour_inds);
%       hFig = figure(n * 10 + 3);
%       clf;
%       contourSpotless3D(subs(Vsol,  [t;x(sub_inds)], [0;0;0]), x(contour_inds), 0, R_diag(contour_inds));
%       xlabel('x_c_m'); ylabel('xdot_c_m'); zlabel('z_c_m');

      % video of rotating ROA
      create_video = false;
      if create_video
        createRotatingVideo(hFig,[class(obj) '_V' num2str(n)]);
      end
    end
    
    function draw(obj,t,x)
      % plotting function
      x_stance = x(end-1);
      % draw line from origin to COM
      h=line(x_stance+[0;x(1)],[0;x(2) + obj.z_nom]);
      set(h,'LineWidth',3,'Color','red')
      axis_major = .3;
      axis_minor = .2;
      theta = linspace(0,2*pi,100);
      x_ellipse = axis_minor*cos(theta);
      z_ellipse = axis_major*sin(theta);
      x_body = x_stance + x(1) + x_ellipse*cos(x(3)) + z_ellipse*sin(x(3));
      z_body = obj.z_nom + x(2) - x_ellipse*sin(x(3)) + z_ellipse*cos(x(3));
      patch(x_body,z_body,'k')
      
      h=line([-10 10],[0 0]);
      set(h,'LineWidth',5,'Color','black')
      radius = .1;
      rectangle('Position',[x_stance+x(1)-radius/2,x(2)+obj.z_nom-radius/2,radius,radius],'Curvature',[1,1], 'FaceColor','k')
      xlim([-.5 1.5])
      ylim([-.5 1.5])
      axis off
    end
  end
end
