function [] = Problem_2()

    %%%%%%
    % Solves the stream function for a 2D jet in self-similar form, using RK2 and RK4.
    %   Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();

    % Initialize times at which to evaluate solution.
    eta0 = 0;
    etaf = 4;
    deta = 0.05;
    rk2.eta  = [eta0:deta:etaf]';
    
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
    
    % Repeat integration using a higher-order method (here, ode45).
    [T, Y] = ode45(@stream, [0,4], [0,1,0]);
    o45.eta = T;
    o45.f   = Y(:,1);
    o45.fp  = Y(:,2);
    o45.fpp = Y(:,3);
    o45.uou0 = o45.fp;
    o45.vou0 = o45.eta .* o45.fp - 0.5 * o45.f;
    
    % Plot fields.
    figure();
    hold on;
    plot(rk2.eta,rk2.uou0,      'DisplayName','RK2: U/U_0');
    plot(o45.eta,o45.uou0, '--','DisplayName','O45: U/U_0');
    plot(rk2.eta,rk2.vou0,      'DisplayName','RK2: V/U_0');
    plot(o45.eta,o45.vou0, '--','DisplayName','O45: V/U_0');
    plot(rk2.eta,rk2.f   ,      'DisplayName','RK2: f');
    plot(o45.eta,o45.f   ,'k--','DisplayName','O45: f');
    plot(rk2.eta,exp(-0.692*rk2.eta.^2),'DisplayName','BEST FIT (f)');
    hleg = legend('show');
    set(hleg, 'location', 'southwest');
    xlim([eta0,etaf]);
    xlabel('\eta');
    
    
%     % Calculate actual solution using best-fit expression.
%     best.eta = rk2.eta;
%     best.f   =                      exp(-0.682*rk2.eta.^2);
%     best.fp  = -0.682   * eta    .* exp(-0.682*rk2.eta.^2);
%     best.fpp =  0.682^2 * eta.^2 .* exp(-0.682*rk2.eta.^2);
%     figure();
%     hold on;
%     plot(rk2.eta,best.f,'DisplayName','f');
%     plot(rk2.eta,best.fp,'DisplayName','f''');
%     plot(rk2.eta,best.fpp,'DisplayName','f''''');
%     legend('show');
    
end

function dy = stream(~, y)
    dy = zeros(3,1);
    dy(1) = y(2);
    dy(2) = y(3);
    dy(3) = -y(1)*y(2);
end





