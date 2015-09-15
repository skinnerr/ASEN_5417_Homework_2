function [] = Problem_3()

    %%%%%%
    % Solves a stiff ODE using an RK4 method and two different step sizes.
    %   Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();

    % Initialize times at which to evaluate solution.
    t0 = 0;
    tf = 10;
    dt = 0.5;
    t  = (t0:dt:tf)';
    
    % Calculate real values analytically.
    ya   = (1/2) * exp(-19*t/2) - (19/2) * exp(-t/2);
    ypa  = (19/4) * (exp(-t/2) - exp(-19*t/2));
    yppa = (19/4) * ((-1/2)*exp(-t/2) - (-19/2)*exp(-19*t/2));
    
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
    
    % Plot fields.
    figure();
    hold on;
%     plot(t,yppa,'DisplayName','y'''' analytic');
%     plot(t,ypa ,'DisplayName','y'' analytic');
    plot(t,ya  ,'DisplayName','y analytic');
%     plot(t,ypp ,'DisplayName','y''''');
%     plot(t,yp  ,'DisplayName','y''');
    plot(t,y   ,'DisplayName','y');
    hleg = legend('show');
    set(hleg, 'location', 'northeast');
    xlim([t0,tf]);
    xlabel('\eta');
    
end