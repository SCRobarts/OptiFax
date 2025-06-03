function y = smoothedTopHat(x, a, b, w)
% SMOOTHTOPHAT Create a smoothed top hat function.
%
%   y = smoothedTopHat(x, a, b, w) returns a vector y of the same size as x.
%   The function smoothly rises from 0 to 1 around x = a, stays near 1 for x
%   between a and b, and then smoothly falls back to 0 around x = b.
%
%   Inputs:
%       x - Vector of x values.
%       a - Lower limit of the central region (transition from 0 to 1).
%       b - Upper limit of the central region (transition from 1 to 0).
%       w - Smoothing parameter that sets the steepness of the transitions.
%           Smaller values produce steeper edges.
%
%   Example:
%       x = linspace(-10, 10, 1000);
%       y = smoothedTopHat(x, -2, 2, 0.5);
%       figure;
%       plot(x, y, 'LineWidth',2);
%       xlabel('x');
%       ylabel('f(x)');
%       title('Smoothed Top Hat Function');
%
%   The implementation is based on subtracting two logistic functions:
%       lower = 1 ./ (1 + exp(-(x - a)/w));
%       upper = 1 ./ (1 + exp(-(x - b)/w));
%       y = lower - upper;
%
%   Alternatively, you could use hyperbolic tangents:
%       y = 0.5*(tanh((x - a)/w) - tanh((x - b)/w));
%

    % Logistic approach:
    % lower = 1 ./ (1 + exp(-(x - a)/w));
    % upper = 1 ./ (1 + exp(-(x - b)/w));
    % y = lower - upper;

    % If you prefer the hyperbolic tangent version, you can comment out the
    % above code and uncomment the following line:
    y = 0.5*(tanh((x - a)/w) - tanh((x - b)/w));
end