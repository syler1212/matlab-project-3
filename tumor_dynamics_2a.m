clear; close all;

%% Setup parameters for mouse lung carcinoma
k = 0.192;                  % Base tumor growth rate constant
alpha = 5.85;               % Stimulatory factor (vascular growth driven by tumor cells)
beta = 0.00873;             % Inhibitory factor (vascular limitation based on carrying capacity and surface area)

% Set up time grid for Euler's method (80 days with a 0.01-day step size)
dt = 0.01;
t = 0:dt:80;
num_steps = length(t);

% Preallocate arrays for tumor cells (N) and dynamic carrying capacity (M)
N = zeros(1, num_steps);
M = zeros(1, num_steps);

% Initial conditions: start with 1 tumor cell and a carrying capacity of 1
N(1) = 1;
M(1) = 1;

%% Solve the coupled differential equations using Euler's Method
for i = 1:(num_steps - 1)
    % Small safety floor to prevent division by zero or log errors if N drops
    if N(i) <= 0
        N(i) = 1e-5;
    end
    
    % Assume a spherical tumor, meaning surface area scales with N^(2/3)
    surface_area = N(i)^(2/3);
    
    % Calculate instantaneous rates of change at the current step
    dNdt = k * N(i) * log(M(i) / N(i));
    dMdt = (alpha * N(i)) - (beta * M(i) * surface_area);
    
    % Step forward in time using Euler's update rule
    N(i + 1) = N(i) + dNdt * dt;
    M(i + 1) = M(i) + dMdt * dt;
end

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
