function [] = Problem_1()

    %%%%%%
    % Solves the heat conduction equation (u_,t = a u_,xx) using the FTCS method.
    %
    % Ryan Skinner, October 2015
    %%%
    
    Set_Default_Plot_Properties();
    
    %%%
    % Define variables specific to the boundary-value problem.
    %%%
    
    % Material properties.
    alpha = 0.1515;
    
    % Solution domain: the closed interval [0,2].
    N = 21;
    x = linspace(0,2,N);
    dx = x(2) - x(1);
    
    % Boundary conditions (0 = initial, f = final).
    BC.u0 = 0;
    BC.uf = 0;
    
    % Initial conditions.
    ux0 = 100 * sin(pi * x / 2);
    
    % Condition numbers, and maximum time to reach in each case.
    d = [0.5, 1.0];
    tmax = [8,4];
    
    % Time steps.
    dt = dx^2 * d / alpha;
    
    % Maximum time steps to reach.
    tsmax = ceil(tmax ./ dt);
    
    % Store solutions indexed by time step, and spatial index.
    u = {nan(tsmax(1),length(x)), nan(tsmax(2),length(x))};
    
    % Initialize solution and set boundary conditions.
    u{1}(1,:) = ux0;
    u{2}(1,:) = ux0;
    u{1}(:,1) = BC.u0;
    u{1}(:,N) = BC.uf;
    u{2}(:,1) = BC.u0;
    u{2}(:,N) = BC.uf;
    
    %%%
    % Solve problem numerically.
    %%%
    
    % Loop over condition numbers.
    for di = 1:length(d)
        fprintf('Working with d = %3.1f, dt = %.4f\n', d(di), dt(di));
        % Iterate over time steps.
        for n = 1:tsmax(di)-1
            % Loop over spatial points.
            for j = 2:N-1
                u{di}(n+1,j) = u{di}(n,j-1) * d(di) ...
                             + u{di}(n,j  ) * (1-2*d(di)) ...
                             + u{di}(n,j+1) * d(di);
            end
        end
    end
    
    %%%
    % Process results: Part (a).
    %%%
    
    % Solution.
    
    ts_plot = 1 + round([0,1,2,4,8] / dt(1));
    t_plot = (ts_plot-1) * dt(1);
    figure();
    hold on;
    for i = 1:length(ts_plot)
        plot(x, u{1}(ts_plot(i),:), ...
            'DisplayName', sprintf('t ~ %.0f',t_plot(i)));
        plot(x, analytical(t_plot(i),x), '--k', ...
            'DisplayName', sprintf('t ~ %.0f (Analytical)',t_plot(i)));
    end
    xlabel('x');
    ylabel('u');
    hleg = legend('show');
    set(hleg, 'Location', 'eastoutside');
    
    % Relative error.
    
    figure();
    hax = axes();
    hold on;
    for i = 1:length(ts_plot)
        relerr = abs(analytical(t_plot(i),x(2:end-1)) - u{1}(ts_plot(i),2:end-1)) ./ ...
                     analytical(t_plot(i),x(2:end-1));
        plot(x(2:end-1), relerr, ...
            'DisplayName', sprintf('t ~ %.0f',t_plot(i)));
    end
    set(hax,'YScale','log');
    xlabel('x');
    ylabel('Relative Error');
    hleg = legend('show');
    set(hleg, 'Location', 'eastoutside');
    
    %%%
    % Process results: Part (b).
    %%%
    
    % Solution.
    
    ts_plot = 1 + round([0,0.1,0.5,1,2,2.1] / dt(2));
    t_plot = (ts_plot-1) * dt(2);
    figure();
    hold on;
    for i = 1:length(ts_plot)
        plot(x, u{2}(ts_plot(i),:), ...
            'DisplayName', sprintf('t ~ %.1f',t_plot(i)));
        if i == length(ts_plot)
            dn = 'Analytical Solutions';
        else
            dn = '';
        end
        plot(x, analytical(t_plot(i),x), '--k', ...
            'DisplayName', dn);
    end
    xlabel('x');
    ylabel('u');
    hleg = legend('show');
    set(hleg, 'Location', 'eastoutside');
    
    % Relative error.
    
    figure();
    hax = axes();
    hold on;
    for i = 1:length(ts_plot)
        relerr = abs(analytical(t_plot(i),x(2:end-1)) - u{2}(ts_plot(i),2:end-1)) ./ ...
                     analytical(t_plot(i),x(2:end-1));
        plot(x(2:end-1), relerr, ...
            'DisplayName', sprintf('t ~ %.1f',t_plot(i)));
    end
    set(hax,'YScale','log');
    xlabel('x');
    ylabel('Relative Error');
    hleg = legend('show');
    set(hleg, 'Location', 'eastoutside');
    
end

function [ u_exact ] = analytical( t, x )
    % Compute analytical solution.
    u_exact = 100 * exp(-0.3738*t) * sin(pi * x / 2);
end













