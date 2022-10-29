% richardson22
% The Richardson extrapolation helps us find an approximation of the derivative of a function at a specific point using multiple approximations of the derivative. This method gives us closer values to the real values with digits of error.

%
% Parameters
% ==========
%     D    The function D takes a function u and two real-values, x and h, and gives out the centered, forward, or 
%             Backwards first or second order derivative.
%             It approximates the derivative.
%   u	The function u is a real-valued function of a real variable x.
%
%   x	The function u is a real-valued function of a real variable x.
%    h       The step size, which is a very small value close to 0.
%
%    N_max   Maximum iterations to see if the Richardson extrapolation converges.

%    eps_abs  Small number that help show the convergence.
%
% Return Value
% ============
%    du   Derivative approximation 


function [du] = richardson22( D, u, x, h, N_max, eps_abs )

%function exceptions

    if ~isscalar( eps_abs ) || (eps_abs <= 0) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument eps_abs is not a positive scalar' ) );
    end
 
    if ~isscalar( N_max ) || (N_max ~= round( N_max )) || (N_max <= 0)
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument N_max is not a positive integer' ) );
    end
 
    if ~isa( D, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument D is not a function handle' ) );
    end
 
    if ~isa( u, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u is not a function handle' ) );
    end

%content of function
    M = zeros([N_max+1 N_max+1]) % zero matrix created
    R(1, 1) = D( u, x, h)        % using D, x, and h, calculate first entry of						% R sub 1,1
                                 
   
    for i=1:N_max                % loop to iterate N_max                              
        R(i + 1, 1) = D( u, x, (h/(2^i)))   
 
        for j=1:i
            R(i + 1, j + 1) = ( ( 4^j ) * R(i + 1, j) - R(i, j ) ) / ( ( 4 ^j ) - 1 );                         % Calculate all entries along the diagonal
        end
 
        if ( abs( R(i + 1, j + 1) - R(i, j) ) < eps_abs )      % If the newest %               entry - previous entry is less that the error criteria eps_abs
            du = R(i + 1, j + 1)  
            return;                              % Then return the newest entry
        end
    end
%    Otherwise, throw an exception saying the function doesn't converge
    throw(MException( 'MATLAB:numeric_exception', ...
                'The Richardson extrapolation has failed to converge' ) );
end

 
