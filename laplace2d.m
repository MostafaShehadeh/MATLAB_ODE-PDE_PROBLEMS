% m3shehad
% laplace2d
% We will use finite difference equations with insulated boundaries to
% approximate the % solution to Laplace's equation in 2D and 3D. 
% This will result in a system of equations that we can solve.
%
% Parameters
% ==========
%    U : 2D Matrix nx * ny
%        And 3D Matrix nx * ny * nz
%
% Return Values
% =============
% U_out: The matrix Uout is the array U, but with all boundary 
% conditions satisfied


function [U_out] = laplace2d( U )

 
% Error and Warning Checking
% ==========================
% The array U will be checked for errors and arguments to ensure that it
% is a 2D matrix. The second error checking would be to ensure that the
% matrix U's boundary is either Dirichlet or an insulated boundary 
% condition. The programme will run as long as these 
% conditions are met.


    if ~ismatrix( U )
        throw( MException( 'MATLAB:invalid_argument', ...
        'U is not a 2D matrix' ) );
    end
    
    if U(1,:) == -Inf + U(end,:) == -Inf + U(:,1) == -Inf + U(:,end) == -Inf >= 1
        throw( MException( 'MATLAB:invalid_argument', ...
        'Boundary conditions must be insulated' ) );
    end
 
    
% Initialization
% ==============
%
% Using the zeroes' function, initialise the initial 2D nx * ny matrix
% based on its initial parameters. We know that there is a mapping
% of m unresolved points.
 
    [n_x, n_y] = size( U );
    U_out = U;
 
% Mapping the unknown points to a unique number from 1 to m
% =========================================================
%
% Examine each point one by one. If a point is equal to -Inf, 
% it is an unknown point: assign it a unique number between 1 and m.
 
% Associate every -Inf from 1 to m to an integer
    u_to_w = zeros( n_x, n_y );
    w_to_u = zeros( 2, n_x * n_y );
    m = 0;
    for ix = 1:n_x
        for iy = 1:n_y
            if U(ix, iy) == -Inf
                m = m + 1;
                u_to_w(ix, iy) = m;
                w_to_u(:, m) = [ix, iy]';
            end
        end
    end
 
% Creating and solving a system of linear equations
% =================================================
%
% Use the spalloc and zeroes functions from class to create a system of 
% sparse linear equations. Then, find the coordinates and modify the 
% matrices for each adjacent point based on whether it's an unknown value,
% an insulated boundary point, etc.

 
    M = spalloc( m, m, 5*m );
    b = zeros( m, 1 );
 
    for iy = 1:m
        c = w_to_u(:,iy);
        p = [c + [-1 0]', c + [1 0]', c + [0 -1]', c + [0 1]'];
        
        for i = p
            % Dirichlet boundary condition
            if U(i(1), i(2)) ~= -Inf && ~isnan(U(i(1), i(2)))
                % In the ith diagonal entry of M substract 1
                M(iy, iy) = M(iy, iy) - 1;
                % In the ith entry of vector b substract the value
                b(iy) = b(iy) - U(i(1), i(2));
                
            % New unknown (-Inf) point
            elseif U(i(1), i(2)) == -Inf
                % In the ith diagonal entry of M substract 1
                M(iy, iy) = M(iy, iy) - 1;
                % In the (i,j)th entry of M add 1.
                M(iy, u_to_w(i(1), i(2))) = M(iy, u_to_w(i(1), i(2))) + 1;
            end
        end
    end
 
% Substituting the values back into the matrix U_out
% ===================================================
%
% Solve and obtain the values by reintroducing the values into the
% matrix Uout via the M=bx equation.

    w = M \ b;
    for k = 1:m
        c = w_to_u(:,k);
        U_out(c(1), c(2)) = w(k, 1);
    end
end
    
    


