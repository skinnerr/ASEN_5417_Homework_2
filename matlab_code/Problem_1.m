function [] = Problem_1()

    %%%%%%
    % Solves the rocket equation using a second-order Runge-Kutta method.
    %   Ryan Skinner, September 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    % Define constants.
    c.mc    = 51.02;    % Rocket casing mass
    c.g     = 9.8;      % Gravitational accel
    c.rho   = 1.23;     % Air density
    c.A     = 0.1;      % Max cross-sectional area
    c.ve    = 360;      % Exhaust speed
    c.CD    = 0.15;     % Drag coefficient
    c.mp0   = 102.04;   % Initial propellant mass
    c.v0    = 0;        % Initial velocity
    c.z0    = 0;        % Initial altitude

    % Initialize times at which to evaluate solution.
    t0 = 0;
    tf = 60;
    dt = 0.1;
    t = t0:dt:tf;
    
    % Initialize velocities and positions.
    v = zeros(length(t),1);
    z = zeros(length(t),1);
    v(1) = c.v0;
    z(1) = c.z0;
    
    % Perform integration using an RK2 method.
    for n = 1:(length(t)-1)
        tn = t(n);
        vn = v(n);
        zn = z(n);
        k1 = dt * (alpha(c,tn)      + beta(c,tn)      * vn   * abs(vn));
        l1 = dt * vn;
        k2 = dt * (alpha(c,tn + dt/2) + beta(c,tn + dt/2) * (vn + k1/2) * abs(vn + k1/2));
        l2 = dt * (vn + l1/2);
        v(n+1) = vn + k2;
        z(n+1) = zn + l2;
    end
    
    % Perform validation using Matlab's ODE45.
    t_span = [t0, tf];
    initials = [v(1), z(1)];
    [t_ml, sol_ml] = ode45(@(t,y) rocket(t,y,c), t_span, initials);
    v_ml = sol_ml(:,1);
    z_ml = sol_ml(:,2);
    
    % Plot velocity.
    figure();
    hold on;
    plot(t_ml,v_ml,'DisplayName','Matlab ODE45');
    plot(t,v,'--','DisplayName','Custom RK2');
    legend('show');
    plot(t_span,[0,0],'k:');
    xlim(t_span);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    
    % Plot altitude.
    figure();
    hold on;
    plot(t_ml,z_ml,'DisplayName','Matlab ODE45');
    plot(t,z,'--','DisplayName','Custom RK2');
    hleg = legend('show');
    set(hleg,'location','south');
    plot(t_span,[0,0],'k:');
    xlim(t_span);
    xlabel('Time (s)');
    ylabel('Altitude (m)');
    
    % Stats on rocket flight (Custom RK2)
    i = find(v == max(v));
    fprintf('RK2  : Max velocity is %.2f at time %.2f and height %.2f.\n',max(v),t(i),z(i));
    i = find(z == max(z));
    fprintf('RK2  : Max altitude is %.2f at time %.2f.\n',max(z),t(i));
    ignore = 100;
    t_crash = interp1(z(ignore:end),t(ignore:end),0);
    fprintf('RK2  : Crash occurs at time %.2f with velocity %.2f\n',t_crash,interp1(t,v,t_crash));
    
    % Stats on rocket flight (Matlab's ODE45)
    i = find(v_ml == max(v_ml));
    fprintf('ODE45: Max velocity is %.2f at time %.2f and height %.2f.\n',max(v_ml),t_ml(i),z_ml(i));
    i = find(z_ml == max(z_ml));
    fprintf('ODE45: Max altitude is %.2f at time %.2f.\n',max(z_ml),t_ml(i));
    ignore = 50;
    t_crash = interp1(z_ml(ignore:end),t_ml(ignore:end),0);
    fprintf('ODE45: Crash occurs at time %.2f with velocity %.2f\n',t_crash,interp1(t_ml,v_ml,t_crash));
    

end

function [ mp_dot ] = mp_dot ( c, t )
% Calculates time derivative of propellant mass.
% This is the exact value of mp(t) given in the problem statement.
    if 0 <= t && t < 1
        mp_dot = t;
    elseif 1 <= t && t < 4
        mp_dot = 1;
    elseif 4 <= t && t < 5
        mp_dot = 5-t;
    else
        mp_dot = 0;
    end
    mp_dot = mp_dot * c.mp0 / 4;
end

function [ mp ] = mp ( c, t )
% Calculates instantaneous propellant mass.
% This is the exact value of mp(t) given in the problem statement.
    if t < 0
        intg = 0;
    elseif 0 <= t && t < 1
        intg = 0.5 * t^2;
    elseif 1 <= t && t < 4
        intg = t - 0.5;
    elseif 4 <= t && t < 5
        intg = 3.5 + (t-4) - 0.5 * (t-4)^2;
    else
        intg = 4;
    end
    mp = c.mp0 - intg * c.mp0 / 4;
end

function [ a ] = alpha( c, t )
% Calculates the time-dependent constant alpha.
    engines_on = 0 <= t && t < 5;
    a = -c.g + (engines_on * c.ve) * mp_dot(c,t) / (mp(c,t) + c.mc);
end

function [ b ] = beta( c, t )
% Calculates the time-dependent constant beta.
    b = -0.5 * c.rho * c.A * c.CD / (mp(c,t) + c.mc);
end

function [ dy ] = rocket( t, y, c )
% Calculates values for ODE45-based integration of the rocket equation.
    dy = zeros(2,1);
    dy(1) = alpha(c,t) + beta(c,t) * y(1) * abs(y(1));
    dy(2) = y(1);
end











