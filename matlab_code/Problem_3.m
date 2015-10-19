function [] = Problem_3()

    %%%%%%
    % Solves the heat conduction equation (u_,t = a u_,xx) using Crank-Nicolson and an
    % adiabatic far boundary condition (Neumann).
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
    BC.upf = 0; % Where up = u prime = du/dx.
    
    % Initial conditions.
    ux0 = 100 * sin(pi * x / 2);
    
    % Condition numbers, and maximum time to reach in each case.
    d = 0.1;
    tmax = 20;
    
    % Time steps.
    dt = dx^2 * d / alpha;
    
    % Maximum time steps to reach.
    tsmax = ceil(tmax ./ dt);
    
    % Store solutions indexed by time step, and spatial index.
    u = nan(tsmax(1),length(x));
    
    % Initialize solution and set boundary conditions.
    u(1,:) = ux0;
    
    %%%
    % Solve problem numerically.
    %%%
    
    fprintf('Working with d = %4.1f, dt = %.4f\n', d, dt);
    
    beta = alpha * dt / (2 * dx^2);
    
    % Iterate over time steps, taking advantage of LU-decomposition efficiency.
    for n = 1:tsmax-1
        [diag, sub, sup, rhs] = Assemble_u_Prob3(beta, u(n,:), BC, dx);
        if n == 1
            [l_vec, u_vec] = LU_Decompose(diag, sub, sup);
        end
        [sol] = LU_Solve(sub, l_vec, u_vec, rhs);
        u(n+1,:) = [BC.u0; sol; sol(end) + dx * BC.upf];
    end
    
    %%%
    % Process results.
    %%%
    
    n_plot = 61;
    cmap = jet(n_plot);
    
    ts_plot = round(linspace(0,tmax,n_plot) / dt);
    ts_plot(1) = ts_plot(1) + 1;
    t_plot = (ts_plot) * dt;
    t_plot(1) = t_plot(1) - dt;

    % Solution.

    hl = cell(length(ts_plot),1);
    name = cell(length(ts_plot),1);

    figure();
    hold on;
    for i = 1:length(ts_plot)
        hl{i} = plot(x, u(ts_plot(i),:),'Color',cmap(i,:));
        name{i} = sprintf('t ~ %.1f',t_plot(i));
    end
    xlabel('x');
    ylabel('u');
    title(sprintf('d = %.1f',d));
    legend_indices = ceil(linspace(1,n_plot,11));
    hleg = legend([hl{legend_indices}], name{legend_indices});
    set(hleg, 'Location', 'eastoutside');
        
end