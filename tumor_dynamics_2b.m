clear; close all;

%% Set up parameters for the treatment model
k = 0.192;                  % Tumor growth rate constant
alpha = 5.85;               % Vascular stimulation parameter
beta = 0.00873;             % Vascular inhibition parameter
y0 = [1; 1];                % Initial conditions: start with 1 cell and M = 1
tspan = [0 80];             % Run the simulation across an 80-day span

%% Run the simulation with a moderate treatment intensity (c = 0.1) starting at day 40
[t1, y1] = ode45(@(t, y) tumor_ode_sys_treatment(t, y, k, alpha, beta, 0.1), tspan, y0);
N1 = y1(:, 1);
M1 = y1(:, 2);

%% Run the simulation again with a stronger treatment intensity (c = 0.2) starting at day 40
[t2, y2] = ode45(@(t, y) tumor_ode_sys_treatment(t, y, k, alpha, beta, 0.2), tspan, y0);
N2 = y2(:, 1);
M2 = y2(:, 2);

%% Plot and compare the tumor cell populations under both treatments
figure('Name', 'Project 2b: Anti-Tumor Treatment Comparison', 'Position', [100, 100, 900, 500]);
plot(t1, N1, 'b-', 'LineWidth', 2);
hold on;
plot(t2, N2, 'r--', 'LineWidth', 2);
xline(40, 'k:', 'LineWidth', 1.5, 'Label', 'Treatment Starts (Day 40)');
xlabel('Time (Days)');
ylabel('Tumor Cell Count ($N$)', 'Interpreter', 'latex');
legend('Treatment ($c = 0.1$)', 'Treatment ($c = 0.2$)', 'Location', 'NorthWest');
title('Comparison of Tumor Growth Under Anti-Tumor Treatment');
grid on;

%% Helper Function: ODE system featuring the conditional treatment term
function dydt = tumor_ode_sys_treatment(t, y, k, alpha, beta, c)
    N = y(1);
    M = y(2);
    
    % Floor N slightly above zero to avoid numerical log errors
    if N <= 0
        N = 1e-5;
    end
    
    % Surface area scales proportionally to N^(2/3) for a spherical tumor
    surface_area = N^(2/3);
    
    % Check if treatment has kicked in yet (starts after day 40)
    if t >= 40
        % Apply standard growth minus the unconstrained treatment kill rate (-cN)
        dNdt = k * N * log(M / N) - (c * N);
    else
        % Standard uninhibited dynamic carrying capacity growth before day 40
        dNdt = k * N * log(M / N);
    end
    
    % Rate of change for carrying capacity (vascular stimulation minus suppression)
    dMdt = (alpha * N) - (beta * M * surface_area);
    
    % Return the column vector of derivatives
    dydt = [dNdt; dMdt];
end
