clear; close all;

%% Setup parameters for mouse lung carcinoma
k = 0.192;                  % Base tumor growth rate constant
alpha = 5.85;               % Stimulatory factor (vascular growth driven by tumor cells)
beta = 0.00873;             % Inhibitory factor (vascular limitation based on carrying capacity and surface area)

% Initial conditions: start with 1 tumor cell and a carrying capacity of 1
y0 = [1; 1]; 

% Run the simulation across an 80-day time span
tspan = [0 80];

%% Solve the coupled differential equations using ode45
[t, y] = ode45(@(t, y) tumor_ode_sys(t, y, k, alpha, beta), tspan, y0);

% Extract tumor cell population (N) and dynamic carrying capacity (M) from the results
N = y(:, 1);
M = y(:, 2);

%% Calculate the per capita growth rate over time
% Derived from the Gompertz equation: (1/N)*(dN/dt) = k * ln(M / N)
per_capita_growth = k * log(M ./ N);

%% Plot tumor population (N) and carrying capacity (M) together
figure('Name', 'Project 2a: N and M Dynamics', 'Position', [100, 100, 800, 500]);
plot(t, N, 'r-', 'LineWidth', 2);
hold on;
plot(t, M, 'b--', 'LineWidth', 2);
xlabel('Time (Days)');
ylabel('Population / Carrying Capacity');
legend('Tumor Cells ($N$)', 'Carrying Capacity ($M$)', 'Location', 'NorthWest', 'Interpreter', 'latex');
title('Tumor Growth and Dynamic Carrying Capacity ($M$) over Time');
grid on;

%% Plot per capita growth rate versus time
figure('Name', 'Project 2a: Per Capita Growth Rate', 'Position', [150, 150, 800, 500]);
plot(t, per_capita_growth, 'k-', 'LineWidth', 2);
xlabel('Time (Days)');
ylabel('Per Capita Growth Rate $\frac{1}{N}\frac{dN}{dt}$', 'Interpreter', 'latex');
title('Per Capita Growth Rate vs. Time');
grid on;

%% Helper Function: Defines the system of ODEs
function dydt = tumor_ode_sys(~, y, k, alpha, beta)
    N = y(1);
    M = y(2);
    
    % Small safety floor to prevent division by zero or log errors if N drops
    if N <= 0
        N = 1e-5;
    end
    
    % Assume a spherical tumor, meaning surface area scales with N^(2/3)
    surface_area = N^(2/3);
    
    % Rate of change for tumor cells based on dynamic carrying capacity M
    dNdt = k * N * log(M / N);
    
    % Rate of change for carrying capacity M (vascular stimulation minus inhibition)
    dMdt = (alpha * N) - (beta * M * surface_area);
    
    % Return the derivatives as a column vector for ode45
    dydt = [dNdt; dMdt];
end
