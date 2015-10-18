function [diag, sub, sup, rhs] = Assemble_u( beta, u_prev, BC )

    %%%%%%
    % Assembles the LHS matrix and the RHS vector for the g-system.
    %   diag -- diagonal
    %    sub -- sub-diagonal
    %    sup -- super-diagonal
    %    rhs -- right-hand side vector
    %
    % Ryan Skinner, October 2015
    %%%
    
    N = length(u_prev);
    
    diag_range = 2:N-1;
     sub_range = 3:N-1;
     sup_range = 2:N-2;
    
    diag = (1 + 2 * beta) * ones(length(diag_range),1);
     sub = (       -beta) * ones(length(sub_range),1);
     sup = (       -beta) * ones(length(sup_range),1);
     rhs = u_prev(diag_range) + ...
           beta * (u_prev(diag_range-1) - 2 * u_prev(diag_range) + u_prev(diag_range+1));
    
    % Account for boundary conditions.
    rhs(1)   = rhs(1)   - (-beta) * BC.u0;
    rhs(end) = rhs(end) - (-beta) * BC.uf;

end