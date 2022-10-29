% Mostafa Shehadeh
% wave1d
% The finite-difference method for the 1D wave equation involves locating
% a matrix U that approximates the PDE.
% Parameters
% ==========
%   c           wave speed
%   x_int       boundary spatial interval [a, b]
%   t_int       interval time [t0, tfinal]
%
%   u_init      The intitial state of the function handle (x)
%   u_bndry     The boundary condition of the function handle (t)
%   du_init     The initial speed of the function handle (x)
%
%   n_x         increment of spatial interval
%   n_t         increment of time interval
%
% Return Values
% =============
%   x_out      values of space in x_int
%   t_out      values of time in t_int
%   U_out      Approximation of the u(x,t) solution
 
function [x_out, t_out, U_out] = wave1d( c, x_int,  n_x, t_int, n_t, u_init, du_init, u_bndry )
 
% Error, arguments, and type Checking
 
    if ~isscalar( c ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'kappa is not a scalar' ) );
    end
    if ~all( size( x_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'x_int is not a 2D row vector' ) );
    end
    if ~all( size( t_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        't_int is not a 2D row vector' ) );
    end
    if ~isscalar( n_x ) || ( n_x ~= round( n_x ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'n_x is not an integer' ) );
    end
    if ~isscalar( n_t ) || ( n_t ~= round( n_t ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'n_t is not an integer' ) );
    end
    if ~isa(u_init, 'function_handle')
        throw( MException( 'MATLAB:invalid_argument', ...
        'u_init is not a funcrion handle' ) );
    end
    if ~isa(u_bndry, 'function_handle')
        throw( MException( 'MATLAB:invalid_argument', ...
        'u_bndry is not a funcrion handle' ) );
    end
    if ~isa(du_init, 'function_handle')
        throw( MException( 'MATLAB:invalid_argument', ...
        'du_init is not a funcrion handle' ) );
    end    
% =============
% Conditions checking:
 
 
    dt = range(t_int) / (n_t - 1);
    h = range(x_int) / (n_x - 1);
    r = (c * dt / h)^2;
    
    if ( r >= 1 )
        throw( MException( 'MATLAB:invalid_argument', ...
        'The ratio r = %d >= 1, use n_t = %d', ...
        r, ceil( range(t_int) * c / h + 1 ) + 1 ) );
    end
    
    % Using the vectors x and t, evaluate the solution. 
    % Using a zeros matrix, initial and boundary values.
 
    x_out = linspace(x_int(1), x_int(2), n_x).';
    t_out = linspace(t_int(1), t_int(2), n_t);
    U_out = zeros(n_x, n_t);
    U_out(:, 1) = u_init(x_out);
    U_out([1, end], :) = u_bndry(t_out);
    
    % Considering the insulated boundary, resolving the issue for the
    % general case and each special case at time t2:
 
    % For t2
    U_out(2:end-1, 2) = U_out(2:end-1, 1) + dt * du_init(U_out(2:end-1, 1));
    % For insulated boundaries:
    if isnan(U_out(1, 2))
        U_out(1, 2) = 4/3 * U_out(2, 2) - 1/3 * U_out(3, 2);
    end
    if isnan(U_out(n_x, 2))
        U_out(end, 2) = 4/3 * U_out(end-1, 2) - 1/3 * U_out(end-2, 2);
    end
    % For rest of columns:
    for k = 3:n_t
        U_out(2:end-1, k) = 2 * U_out(2:end-1, k-1) - U_out(2:end-1, k-2) + r * diff(U_out(:, k-1), 2);  
    % For insulated boundary conditions:
        if isnan(U_out(1, k))
            U_out(1, k) = 4/3 * U_out(2, k) - 1/3 * U_out(3, k);
        end
        if isnan(U_out(end, k))
            U_out(end, k) = 4/3 * U_out(end-1, k) - 1/3 * U_out(end-2, k);
        end
    end
end
