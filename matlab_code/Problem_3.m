function [] = Problem_3()

    %%%%%%
    % Solves a stiff ODE using an RK4 method and two different step sizes.
    %   Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();

    % Initialize times at which to evaluate solution (with dt = hmax*2).
    t0 = 0;
    tf = 15;
    dt = 8/19;
    t  = (t0:dt:tf)';
    
    % Initialize velocities and positions.
    y   = zeros(length(t),1);
    yp  = zeros(length(t),1);
    ypp = zeros(length(t),1);
    y(1)   = -9;
    yp(1)  = 0;
    ypp(1) = 42.75;
    
    % Perform integration using an RK4 method.
    for n = 1:(length(t)-1)
        yn   = y(n);    % k
        ypn  = yp(n);   % l
        yppn = ypp(n);  % m
        k1 = dt * ypn;
        l1 = dt * yppn;
        m1 = dt * -(19/4)*yn - 10*ypn;
        k2 = dt * (ypn + l1/2);
        l2 = dt * (yppn + m1/2);
        m2 = dt * -(19/4)*(yn + k1/2) - 10*(ypn + l1/2);
        k3 = dt * (ypn + l2/2);
        l3 = dt * (yppn + m2/2);
        m3 = dt * -(19/4)*(yn + k2/2) - 10*(ypn + l2/2);
        k4 = dt * (ypn + l3);
        l4 = dt * (yppn + m3);
        m4 = dt * -(19/4)*(yn + k3) - 10*(ypn + l3);
        y(n+1)   = yn   + k1/6 + k2/3 + k3/3 + k4/6;
        yp(n+1)  = ypn  + l1/6 + l2/3 + l3/3 + l4/6;
        ypp(n+1) = yppn + m1/6 + m2/3 + m3/3 + m4/6;
    end
    hmt2.t = t;
    hmt2.y = y;
    hmt2.y10 = interp1(hmt2.t,hmt2.y,10);

    % Initialize times at which to evaluate solution (with dt = hmax/2).
    t0 = 0;
    tf = 15;
    dt = 2/19;
    t  = (t0:dt:tf)';
    
    % Initialize velocities and positions.
    y   = zeros(length(t),1);
    yp  = zeros(length(t),1);
    ypp = zeros(length(t),1);
    y(1)   = -9;
    yp(1)  = 0;
    ypp(1) = 42.75;
    
    % Perform integration using an RK4 method.
    for n = 1:(length(t)-1)
        yn   = y(n);    % k
        ypn  = yp(n);   % l
        yppn = ypp(n);  % m
        k1 = dt * ypn;
        l1 = dt * yppn;
        m1 = dt * -(19/4)*yn - 10*ypn;
        k2 = dt * (ypn + l1/2);
        l2 = dt * (yppn + m1/2);
        m2 = dt * -(19/4)*(yn + k1/2) - 10*(ypn + l1/2);
        k3 = dt * (ypn + l2/2);
        l3 = dt * (yppn + m2/2);
        m3 = dt * -(19/4)*(yn + k2/2) - 10*(ypn + l2/2);
        k4 = dt * (ypn + l3);
        l4 = dt * (yppn + m3);
        m4 = dt * -(19/4)*(yn + k3) - 10*(ypn + l3);
        y(n+1)   = yn   + k1/6 + k2/3 + k3/3 + k4/6;
        yp(n+1)  = ypn  + l1/6 + l2/3 + l3/3 + l4/6;
        ypp(n+1) = yppn + m1/6 + m2/3 + m3/3 + m4/6;
    end
    hmd2.t = t;
    hmd2.y = y;
    hmd2.y10 = interp1(hmd2.t,hmd2.y,10);
    
    % Calculate solution analytically.
    ya = (1/2) * exp(-19*hmd2.t/2) - (19/2) * exp(-hmd2.t/2);
    ya10 = (1/2) * exp(-19*10/2) - (19/2) * exp(-10/2);
    
    % Plot fields.
    figure();
    hold on;
    plot(hmd2.t,    ya,     'DisplayName','analytic');
    plot(hmt2.t,hmt2.y,'--','DisplayName','RK4 (hmax*2)');
    plot(hmd2.t,hmd2.y,'-.','DisplayName','RK4 (hmax/2)');
    hleg = legend('show');
    set(hleg, 'location', 'eastoutside');
    xlim([0,10]);
    xlabel('\eta');
    ylabel('y');
    
    % Print statistics.
    fprintf('Values of y(10) for each method:\n');
    fprintf('RK4 (hmax/2): %.2f\n',hmd2.y10);
    fprintf('RK4 (hmax*2): %.2f\n',hmt2.y10);
    fprintf('Actual      : %.2f\n',ya10);
    
end





