%This script contains the function to solve the Colebrook-White equation using the Newton-Raphson method. An initial guess for the friction factor is made using the Swamee-Jain equation, and the function iteratively refines this guess until it converges to a solution within the specified tolerance or reaches the maximum number of iterations.

function [f, f_values] = colebrookFunction(Re, relativeRoughness, tolerance, maxIterations)
    if nargin < 3
        tolerance = 1e-6; % Default tolerance for convergence
    end
    if nargin < 4
        maxIterations = 100; % Default maximum number of iterations
    end
    
    %functions for easier use
    w = (relativeRoughness / 3.7) + (5.74 / (Re^0.9)); %For use in Swamee-Jain

    %Initial guess for the friction factor using the Swamee-Jain equation
    f0 = 0.25 / (log10(w)^2);

    %set boolean flag for convergence   
    converged = false;

    %create an array to store the friction factor values for each iteration 
    f_values = zeros(1, maxIterations);

    for i = 1:maxIterations
        u = (relativeRoughness / 3.7) + (2.51 / (Re * sqrt(f0))); %For use in Colebrook-White   
        %Calculate f using g(x) and g'(x)
        g = 1/sqrt(f0) + 2*log10(u);  
        g_prime = -0.5 * f0^(-1.5) + (5.02/(u * log(10) * Re));
        f_new = f0 - g / g_prime; 
        f_values(i) = f_new; % Store the current friction factor value

        if abs(f_new - f0) < tolerance
            converged = true;   
            break; % Convergence achieved
        end 
        f0 = f_new; 
    end

    f_values = f_values(1:i); % Return only the computed values up to convergence or max iterations 
     
    if ~converged
        warning('Colebrook-White function did not converge within the maximum number of iterations');   
    end   
    f = f_new; % Return the final friction factor 
end





