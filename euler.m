% Euler’s method
%
% Euler's method is a numerical method for approximating the solution to a 1st order IVP. Then we will add the trapezoidal rule to improve our code.
%
% Parameters
% ==========
%    f first derivative of y with respect to t
%    t_rng row vector with t0 and tfinal as first and second values
%    y0 intianal condition where y is evaluated at t0
%    n number of points  from t0 to tfinal
% Return Values
% =============
%   t_out	row vector of n equally spaced from t0 o tfinal
%   y_out approximation of y from t to tout, where first value is y0 and last %isyfinal

function [ t_out, y_out ] = euler( f, t_rng, y0, n )
    % Argument Checking
    if ~isa(f, 'function_handle')
        error('MATLAB:invalid_argument',...
            'The argument f is not a function handle');
    end
    if ~all(size(t_rng) == [1, 2])
        error('MATLAB:invalid_argument',...
            'The argument t_rng is not 2-dimensional row vector');
    end
    if ~isscalar(y0)
        error('MATLAB:invalid_argument',...
            'The argument y0 is not a scalar');
    end
    if ~isscalar(n) || (n ~= round(n))
        error('MATLAB:invalid_argument',...
            'The argument n is not a positive integer');
    end
    
   
    h = ( t_rng( 2 ) - t_rng( 1 ) ) / ( n - 1 ); % stepsize calculation
    
    t_out = zeros(1, n); % Setting the 1 by n matrix
    t_out( 1 ) = t_rng( 1 ); % initializing tout
    
    for i = 2 : n
       t_out( i ) = t_out( i - 1 ) + h;
    end
    
    y_out = zeros(1, n); % Setting 1 by n zeros matrix
    y_out( 1 ) = y0; % initializing yout
    
    for k = 1 : ( n - 1 )
        m = f( t_out( k ), y_out( k ) ); % Finding the slope
        y_out( k + 1 ) = y_out( k ) + h * m;
    end
end

