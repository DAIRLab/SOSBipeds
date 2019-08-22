function resetVisualization()
    % Constants
    g = 10;
    z_nom = 1;
    step_time = 0.3;
    cop_max = 0.05;
    u_wrt_cm = 0;
  
    model = LIPMSwingLeg(g, z_nom, step_time, cop_max, u_wrt_cm);

    x = msspoly('x', 3);
    t = msspoly('t', 1);
    u = msspoly('u', 2);
    f = model.dynamics(t,x,u);
    
    A = double(subs(diff(f,x),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));
    B = double(subs(diff(f,u),[t;x;u],zeros(1+model.num_states+model.num_inputs,1)));

    Q = eye(model.num_states);
    R = eye(model.num_inputs);
    [K,Q] = lqr(A,B,Q,R);
    
    x = sym('x', [3 1], 'real');
    r = model.reset(t, x, []);
    V = x'*Q*x;
    figure;
    fimplicit3(V - 1);
    hold on;
    
    V_r = subs(V, x, r);
    fimplicit3(V_r - 1);
    
%     p = Polyhedron('Ae', [0 0 1], 'be', 0);
%     p.plot();
%     l = Polyhedron('Ae', [2 0 -1], 'be', 0);
%     l.plot();
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    return;
    X = []; Y = []; Z = [];
    X_r = []; Y_r = []; Z_r = [];
    for i=linspace(-1, 1, 20)
        for j=linspace(-1, 1, 20)
           V_subs = subs(V, x(1:2), [i;j]);
           zs = vpa(solve(V_subs == 1, x(3)));
           if ~isempty(zs)
               for z=zs'
                   X = [X, i];
                   Y = [Y, j];
                   Z = [Z, z];
                   X_r = [X_r, i-z];
                   Y_r = [Y_r, j];
                   Z_r = [Z_r, -z];
               end
           end
%           V_r_subs = subs(V_r, x, [i;j;z]);
        end
    end
    figure;
    scatter3(X, Y, Z, 'r');
    hold on;
    scatter3(X_r, Y_r, Z_r, 'b');
    return
    
%     r = 1;
    
%     N = 20;
%     thetavec = linspace(0,pi,N);
%     phivec = linspace(0,2*pi,2*N);
%     [th, ph] = meshgrid(thetavec,phivec);
%     R = r*ones(size(th));

%     x = reshape(R.*sin(th).*cos(ph), 1, []);
%     y = reshape(R.*sin(th).*sin(ph), 1, []);
%     z = reshape(R.*cos(th), 1, []);
    X = [x; y; z];
    
    xp = x-z;
    yp = y;
    zp = -1*z;
    Xp = [xp; yp; zp];
    
    normXp = sum(Xp.^2, 1);
    notX = X(:, normXp >= r^2);
    X(:, normXp >= r^2) = [];
    Xp(:, normXp >= r^2) = [];
    
    figure;
    scatter3(X(1, :), X(2, :), X(3, :), 'b');
%     hold on;
%     scatter3(Xp(1, :), Xp(2, :), Xp(3, :), 'r');
    axis vis3d
    
    X = X';
    notX = notX';
    Y = [];
    t = ones(size(X, 1), 1);
    t = [t; 2*ones(size(notX, 1), 1)];
    
%     X = [1, 2, 3]
%     notX = [4 5 6]
%     for k=0:1
%         for j=0:1
%             for i=0:1
%                 Z = [(X(:,1).^i).*(X(:,2).^j).*(X(:,3).^k); (notX(:,1).^i).*(notX(:,2).^j).*(notX(:,3).^k)];
%                 Y = [Y, Z];
%             end
%         end
%     end
%     Y
   
%     [B, dev, stats] = mnrfit([X; notX], t);
%     B
hold on;
p = Polyhedron('Ae', [0 0 1], 'be', 0);
p.plot();
l = Polyhedron('Ae', [2 0 -1], 'be', 0);
l.plot();
xlabel('x')
ylabel('y')
zlabel('z')

all(([0 0 1]*X').*([2 0 -1]*X') > 0)
all(([0 0 1]*notX').*([2 0 -1]*notX') < 0)
end
