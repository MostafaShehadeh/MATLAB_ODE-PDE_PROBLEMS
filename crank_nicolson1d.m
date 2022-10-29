% crank_nicolson1d
%
% The Crank-Nicolson method finds a matrix U that approximates 
% the PDE using finite-difference to approximate solutions to the
% heat-conduction/diffusion equation.
% 
% Parameters:
% ==========
%    kappa: diffusion coefficient
%    x_rng: interval of space [a, b]
%    t_rng: interval of time [t0, tfinal]
%
%    u_init:  intitial state function
%    u_bndry: Boundary condition function
%
%    nx: increment for spacial interval
%    nt: increment for time interval
%
% Return Values:
% ==========
%    x_out: the vector x is values  of in nx space from a to b
%    t_out: the vector t is values  of in nt time from t initial to t final
%    U_out: the matrix U is an nx*nt solution that approximates u(x,t) 
 
 
function [x_out, t_out, U_out] = crank_nicolson1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )
 
% Error Checking:
% Argument Checking:
% =================
% Validate types and values of function arguments:
    if ~isscalar( kappa )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument x is not a scalar' ) );
    end
    if ~all( size( x_rng ) == [1, 2] )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument x_rng is not a 2-dimensional row vector' ) );
    end
    if ~isscalar( nx ) || ( nx ~= round( nx ) )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument nx is not an integer' ) );
    end
    if ~all( size( t_rng ) == [1, 2] )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument t_rng is not a 2-dimensional row vector' ) );
    end
    if ~isscalar( nt ) || ( nt ~= round( nt ) )
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument nt is not an integer' ) );
    end
    if ~isa(u_init, "function_handle")
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument u_init is not a funcrion handle' ) );
    end
    if ~isa(u_bndry, "function_handle")
        throw( MException( 'MATLAB:invalid_argument', ...
            'the argument u_bndry is not a funcrion handle' ));
    end
    % ==============
    % validate the types and the values of the arguments,
    % ensure that coefficient is not too large, (K?t/h^2) < 0.5.11
    dt = range(t_rng) / (nt -1);
    h = range(x_rng) / (nx - 1);
    error_check = kappa * dt / (h^2);
    
    if ( error_check >= 0.5 )
        warning( 'MATLAB:questionable_argument', ...
        'The ratio = %d >= 0.5, nt should be %e', ...
        error_check, ceil( 2*kappa*range(t_rng)/ h^2 + 1));

    end

    
    % Initialization:
    % ==============
    % Using the x and t vectors, evaluate the solution by making a
    % zeros
    % matrix and filling in the first column with initial values and
    % the first and last rows with boundary values.
 
    x_out = linspace(x_rng(1), x_rng(2), nx).';
    t_out = linspace(t_rng(1), t_rng(2), nt);
    U_out = zeros(nx, nt);
    U_out(:,1) = u_init(x_out);    
    U_out([1,nx],2:end) = u_bndry(t_out(2:end));
    
    % Tri-diagonal Matrix setup:
    
    d_low = diag( (-error_check) * ones(nx - 3, 1), -1   );
    d_middle =   diag( (2*(1+error_check)) * ones(nx - 2, 1)  );
    d_up = diag( (-error_check) * ones(nx - 3, 1), 1    );
    
    diag_tri = d_low + d_middle + d_up;
 
    for k = 2:nt
    
        b = 2 .* U_out(2:end-1, k-1) + error_check * diff(U_out(:, k-1), 2);
        % Getting the b vector 
        if isnan(U_out(1,k)) && isnan(U_out(end,k))
        % If both boundaries insulated change to tri-diagonal
            diag_tri(1,1) = 2 + error_check * (2/3);
            diag_tri(1,2) = error_check * (-2/3);
            diag_tri(end,end) = 2 + error_check * (2/3);
            diag_tri(end,end-1) = error_check * (-2/3);
        elseif isnan(U_out(1,k))
        % If first is insulated then apply:    
            diag_tri(1,1) = 2 + error_check * (2/3);
            diag_tri(1,2) = error_check * (-2/3);
            b(end) = b(end) + error_check .* U_out(end,k);   
        % If second is insulated then apply:
        elseif isnan(U_out(end,k))
            diag_tri(end,end) = 2 + error_check * (2/3);
            diag_tri(end,end-1) = error_check * (-2/3);
            b(1) = b(1) + error_check .* U_out(1,k);    
        else
            b(1) = b(1) + error_check .* U_out(1,k);
            b(end) = b(end) + error_check .* U_out(end,k);
        end
        % Matrix solving:
        U_out(2:end-1,k) = diag_tri \ b;
        % Setup matrix first element:
        if isnan(U_out(1,k))
            U_out(1,k) = U_out(2,k) * (4/3) - U_out(3,k) / 3;
        end
        % Setup matrix last element:
        if isnan(U_out(end,k))
            U_out(end,k) = U_out(end-1,k) * (4/3) - U_out(end-2,k) / 3;
        end
              
    end
    
end
 




