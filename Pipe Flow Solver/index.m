%This code aims to solve the Colebrrok-White equation using the Newton-Raphson method. We will use that result in a Darcy-Weisbach pressure drop calculation. We will compare the results to the Moody chart. 

%Define the parameters  
D = 0.05; % Diameter of the pipe in meters   
L = 100; % Length of the pipe in meters
epsilon = 0.046e-3; % Roughness of the pipe in meters
nu = 1e-6; % Kinematic viscosity of the fluid in m^2/s

%Velocity range for the calculations    
V = linspace(0.1, 10, 100); % Velocity range from 0.1 m/s to 10 m/s with 100 points

%finding the Reynolds number for each velocity  
Re = (V * D) / nu; 

%finding the relative roughness 
relativeRoughness = epsilon / D;    

f_values = zeros(size(V)); % Preallocate an array to store friction factor values       

%Solve for each Reynolds number
for k = 1:length(V)
    if Re(k) < 2300
        f_values(k) = 64 / Re(k); % Laminar flow
    else
        f_values(k) = colebrookFunction(Re(k), relativeRoughness, 1e-8, 100); % Turbulent flow using Colebrook-White function  
    end
end

rho = 1000; % Density of water in kg/m^3

% Calculate the pressure drop using the Darcy-Weisbach equation 
dP = f_values .* (L/D) .* (0.5 * rho * V.^2);   

%solve headloss
g = 9.81;
hL = dP / (rho * g); % Head loss in meters

%asymptote
fRoughLimit = 1 / (-2*log10(relativeRoughness/3.7))^2;

% %display values for debugging
% disp('Velocity (m/s):');
% disp(V);
% disp('Reynolds Number:');   
% disp(Re);
% disp('Friction Factor:');
% disp(f_values);
% disp('Pressure Drop (Pa):');    
% disp(dP);   

% Plotting the results
subplot(1,2,1)
loglog(Re, f_values, 'b-', 'LineWidth', 2);
hold on

legendEntries = {'Newton-Raphson Solution'};

if any(Re < 2300)
    loglog(Re(Re < 2300), 64 ./ Re(Re < 2300), 'g--', 'LineWidth', 1.5);
    legendEntries{end+1} = 'Laminar Flow (f = 64/Re)';
end

yline(fRoughLimit, 'r--', 'LineWidth', 1.5);
legendEntries{end+1} = 'Roughness Limit Asymptote';

xlabel('Reynolds Number (Re)');
ylabel('Friction Factor (f)');
title('Friction Factor vs Reynolds Number - Colebrook-White Solution');
legend(legendEntries);
grid on

% Plot 2 - Newton-Raphson Convergence
[f_converged, f_iter] = colebrookFunction(1e5, relativeRoughness, 1e-8, 50);
errors = abs(f_iter - f_converged);

subplot(1,2,2)
semilogy(1:length(errors), errors, 'bo-', 'LineWidth', 1.5)
xlabel('Iteration')
ylabel('Absolute Error')
title('Newton-Raphson Convergence - Colebrook-White (Re = 10^5)')
grid on
