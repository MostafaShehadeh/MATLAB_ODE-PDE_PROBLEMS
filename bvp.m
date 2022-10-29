% Boundary Value Problem (BVP)
% Approximating numerical solutions by creating a matrix in a BVP using
% linear ordinary differential equations (ODEs)
%
% Parameters
% ==========
%    c: Constant Coefficients (these are known)
%    x_int: End points (a<b)
%    u_int: Constraints for each point a and b
%
%    g: Forcing function (Known)
%
%    n: Equally spaced interval points
%
% Return Values
% =============
%    x: The vector xout is a column vector of n values
%    u: The vector uout is a column vector of u values

function [x, u] = bvp( c, x_int, u_int, g, n )

if ~all( size( c ) == [1, 3] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument c is not a 3-dimensional row vector' ) );
end
if ~all( size( x_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument x_int is not a 2-dimensional row vector' ) );
end
if ~all( size( u_int ) == [1, 2] ) 
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument u_int is not a 2-dimensional row vector' ) );
end
if ~isa( g, 'function_handle' )
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument g is not a function handle' ) );
end
if ~isscalar( n ) || ( n ~= round( n ) )  
        throw( MException( 'MATLAB:invalid_argument', ...
        'the argument n is not an integer' ) );
end

%Initialization: h takes the value of the (b-a) / (n-1).
h = (x_int(2) - x_int(1))/(n-1);
%Define d-, d, d+ (from slides).
d_minus = 2 * c(1) - h * c(2);
d_ = 2 * h^2 * c(3) - 4 *c(1);
d_plus = 2 * c(1) + h * c(2); 
%With a n increment, create a linspace of x values between a and b.
x_linspace = linspace(x_int(1), x_int(2), n)';
%Remove a and b from linspace.
x = x_linspace(2:end -1);
%Create a vector of u values with n columns.
g_ = 2 * (h^2) * g(x);
%Subtract first value with d- and last value with d+.
g_(1)   = g_(1)   - d_minus * u_int(1);
g_(end) = g_(end) - d_plus * u_int(2);
%Create a tri-diagonal matrix (n-2)*(n-2).
dminus_diag = diag(d_minus * ones(n - 3, 1), -1); 
d_diag = diag(d_ * ones(n - 2, 1));    
dplus_diag = diag(d_plus * ones(n - 3, 1), 1);  
tri_diag = dminus_diag + d_diag + dplus_diag;
%Solve system using Ax= u.
u_int_ = tri_diag \ g_;
x = x_linspace';
u = [u_int(1) u_int_' u_int(2)];
end



	




