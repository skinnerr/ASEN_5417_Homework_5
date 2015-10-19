function [] = Problem_2()

    %%%%%%
    % Solves the heat conduction equation (u_,t = a u_,xx) using Crank-Nicolson and
    % Dirichlet boundary conditions at both ends of the domain.
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
    d = [0.5, 1.0, 10.0];
    tmax = [8,8,8];
    
    % Time steps.
    dt = dx^2 * d / alpha;
    
    % Maximum time steps to reach.
    tsmax = ceil(tmax ./ dt);
    
    % Store solutions indexed by time step, and spatial index.
    u = {nan(tsmax(1),length(x)), nan(tsmax(2),length(x)), nan(tsmax(3),length(x))};
    
    % Initialize solution and set boundary conditions.
    u{1}(1,:) = ux0;
    u{2}(1,:) = ux0;
    u{3}(1,:) = ux0;
    u{1}(:,1) = BC.u0;
    u{1}(:,N) = BC.uf;
    u{2}(:,1) = BC.u0;
    u{2}(:,N) = BC.uf;
    u{3}(:,1) = BC.u0;
    u{3}(:,N) = BC.uf;
    
    %%%
    % Solve problem numerically.
    %%%
    
    % Loop over condition numbers.
    for di = 1:length(d)
        fprintf('Working with d = %4.1f, dt = %.4f\n', d(di), dt(di));
        beta = alpha * dt(di) / (2 * dx^2);
        % Iterate over time steps.
        for n = 1:tsmax(di)-1
            [diag, sub, sup, rhs] = Assemble_u(beta, u{di}(n,:), BC);
            [sol] = Thomas(diag, sub, sup, rhs);
            u{di}(n+1,:) = [BC.u0; sol; BC.uf];
        end
    end
    
    %%%
    % Process results.
    %%%
    
    for di = 1:length(d)
    
        ts_plot = 1 + round([0,1,2,4,8] / dt(di));
        t_plot = (ts_plot-1) * dt(di);
    
        % Solution.

        hl = cell(length(ts_plot)+1,1);
        name = cell(length(ts_plot)+1,1);

        figure();
        hold on;
        for i = 1:length(ts_plot)
            hl{i} = plot(x, u{di}(ts_plot(i),:));
            name{i} = sprintf('t ~ %.1f',t_plot(i));
            if i ~= length(ts_plot)
                plot(x, analytical(t_plot(i),x), '--k');
            else
                hl{i+1} = plot(x, analytical(t_plot(i),x), '--k');
                name{i+1} = 'Analytical';
            end
        end
        xlabel('x');
        ylabel('u');
        title(sprintf('d = %.1f',d(di)));
        hleg = legend([hl{:}], name);
        set(hleg, 'Location', 'eastoutside');

        % Relative error.

        figure();
        hax = axes();
        hold on;
        for i = 1:length(ts_plot)
            relerr = abs(analytical(t_plot(i),x(2:end-1)) - u{di}(ts_plot(i),2:end-1)) ./ ...
                         analytical(t_plot(i),x(2:end-1));
            plot(x(2:end-1), relerr, ...
                'DisplayName', sprintf('t ~ %.1f',t_plot(i)));
        end
        set(hax,'YScale','log');
        xlabel('x');
        ylabel('Relative Error');
        title(sprintf('d = %.1f',d(di)));
        ylim([5*10^-4,10^-2]);
        
    end
    
end

function [ u_exact ] = analytical( t, x )
    % Compute analytical solution.
    u_exact = 100 * exp(-0.3738*t) * sin(pi * x / 2);
end













