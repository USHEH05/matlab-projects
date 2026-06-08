%This code aims to solve the Colebrrok-White equation using the Newton-Raphson method. We will use that result in a Darcy-Weisbach pressure drop calculation. Input Pipe Geometry, Fluid Properties, and flow velocity. Outputs the friction factor, pressure drop, and head loss. 

fprintf('Pipe Flow Analysis using Colebrook-White Equation\n');

%user inputs
D = input('Enter the diameter of the pipe in meters: ');
L = input('Enter the length of the pipe in meters: ');
epsilon = input('Enter the roughness of the pipe in meters: ');
nu = input('Enter the kinematic viscosity of the fluid in m^2/s [1e-6 for water at 20°C]: ');
V = input('Enter the velocity of the fluid in m/s: ');
rho = input('Enter the density of the fluid in kg/m^3 [1000 for water]: ');

% Calculations
Re = (V * D) / nu; % Reynolds number    
relativeRoughness = epsilon / D; % Relative roughness   

if Re < 2300
    f = 64 / Re; % Laminar flow
    fprintf('Flow Regime: Laminar\n');
else
    f = colebrookFunction(Re, relativeRoughness, 1e-8, 100); % Turbulent flow using Colebrook-White function  
    fprintf('Flow Regime: Turbulent\n');
end

dP = f * (L/D) * (0.5 * rho * V^2); 
g = 9.81; 
hL = dP / (rho * g); 

% Display results   
fprintf('Reynolds Number: %.2e\n', Re);
fprintf('Friction Factor: %.6f\n', f);
fprintf('Pressure Drop (Pa): %.2f\n', dP);
fprintf('Head Loss (m): %.2f\n', hL);       


