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
    
%     t = linspace(-1,6,1000);
%     x = zeros(length(t),1);
%     for i = 1:length(t)
%         x(i) = mp(c,t(i));
%     end
%     plot(t,x);
    
%     a = alpha(c,1)
%     b = beta(c,1)

    % Initialize times at which to evaluate solution.
    t0 = 0;
    tf = 85;
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
    
    % Plot velocity.
    figure();
    hold on;
    plot(t,v);
    plot([t0,tf],[0,0],'k--');
    xlim([t0,tf]);
    xlabel('Time (s)');
    ylabel('Velocity (m/s)');
    
    % Plot altitude.
    figure();
    hold on;
    plot(t,z);
    plot([t0,tf],[0,0],'k--');
    xlim([t0,tf]);
    xlabel('Time (s)');
    ylabel('Altitude (m)');
   
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
    a = -c.g + (engines_on * c.ve) * mp(c,t) / (mp(c,t) + c.mc);
end

function [ b ] = beta( c, t )
% Calculates the time-dependent constant beta.
    b = -0.5 * c.rho * c.A * c.CD / (mp(c,t) + c.mc);
end