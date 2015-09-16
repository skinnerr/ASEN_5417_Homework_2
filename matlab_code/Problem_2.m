function [] = Problem_2()

    %%%%%%
    % Solves the stream function for a 2D jet in self-similar form, using RK2 and ode45.
    %   Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();

    % Initialize times at which to evaluate solution.
    eta0 = 0;
    etaf = 4;
    deta = 0.05;
    rk2.eta  = (eta0:deta:etaf)';
    
    % Initialize velocities and positions.
    rk2.f   = zeros(length(rk2.eta),1);
    rk2.fp  = zeros(length(rk2.eta),1);
    rk2.fpp = zeros(length(rk2.eta),1);
    rk2.f(1)   = 0;
    rk2.fp(1)  = 1;
    rk2.fpp(1) = 0;
    
    % Perform integration using an RK2 method.
    for n = 1:(length(rk2.eta)-1)
        f   = rk2.f(n);    % k
        fp  = rk2.fp(n);   % l
        fpp = rk2.fpp(n);  % m
        k1 = deta * fp;
        l1 = deta * fpp;
        m1 = deta * -(f * fp);
        k2 = deta * (fp + l1/2);
        l2 = deta * (fpp + m1/2);
        m2 = deta * -(f + k1/2) * (fp + l1/2);
        rk2.f(n+1)   = f + k2;
        rk2.fp(n+1)  = fp + l2;
        rk2.fpp(n+1) = fpp + m2;
    end
    rk2.uou0 = rk2.fp;
    rk2.vou0 = rk2.eta .* rk2.fp - 0.5 * rk2.f;
    
    % Repeat integration using a higher-order method (ode45).
    [T, Y] = ode45(@stream, [0,4], [0,1,0]);
    o45.eta = T;
    o45.f   = Y(:,1);
    o45.fp  = Y(:,2);
    o45.fpp = Y(:,3);
    o45.uou0 = o45.fp;
    o45.vou0 = o45.eta .* o45.fp - 0.5 * o45.f;
    
    % Calculate solution using best-fit expression.
    bst.eta = rk2.eta;
    bst.f   = 1.07313 * erf(0.825833 * bst.eta);
    bst.fp  =                     exp(-0.682*bst.eta.^2);
    bst.fpp = -0.682 * bst.eta .* exp(-0.682*bst.eta.^2);
    bst.uou0 = bst.fp;
    bst.vou0 = bst.eta .* bst.fp - 0.5 * bst.f;
    
    % Plot fields.
    figure();
    hold on;
    plot(rk2.eta, rk2.uou0,       'DisplayName', 'U/U_0 (RK2)');
    plot(o45.eta, o45.uou0, '--', 'DisplayName', 'U/U_0 (ode45)');
    plot(bst.eta, bst.uou0, '-.', 'DisplayName', 'U/U_0 (Best Fit)');
    plot(rk2.eta, rk2.vou0,       'DisplayName', 'V/U_0 (RK2)');
    plot(o45.eta, o45.vou0, '--', 'DisplayName', 'V/U_0 (ode45)');
    plot(bst.eta, bst.vou0, '-.', 'DisplayName', 'V/U_0 (Best Fit)');
    plot(rk2.eta, rk2.f,          'DisplayName', 'f (RK2)');
    plot(o45.eta, o45.f,    '--', 'DisplayName', 'f (ode45)');
%     plot(bst.eta, bst.f,    '-.', 'DisplayName', 'f (Best Fit)');
    hleg = legend('show');
    set(hleg, 'location', 'southwest');
    xlim([eta0,etaf]);
    ylim([-1,1.5]);
    xlabel('\eta');
    
end

function dy = stream(~, y)
    dy = zeros(3,1);
    dy(1) = y(2);
    dy(2) = y(3);
    dy(3) = -y(1)*y(2);
end