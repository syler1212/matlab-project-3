clear; close all;

%% Set up parameters for the treatment model
k = 0.192;                  % Tumor growth rate constant
alpha = 5.85;               % Vascular stimulation parameter
beta = 0.00873;             % Vascular inhibition parameter

% Set up time grid for Euler's method (80 days with a 0.01-day step size)
dt = 0.01;
t = 0:dt:80;
num_steps = length(t);

%% Run Simulation 1: Moderate treatment intensity (c = 0.1) starting at day 40
c1 = 0.1;
N1 = zeros(1, num_steps);
M1 = zeros(1, num_steps);
N1(1) = 1;
M1(1) = 1;

for i = 1:(num_steps - 1)
    if N1(i) <= 0, N1(i) = 1e-5; end
    
    surface_area = N1(i)^(2/3);
    
    % Apply standard growth before day 40, then subtract treatment kill rate (-c*N)
    if t(i) >= 40
        dNdt = k * N1(i) * log(M1(i) / N1(i)) - (c1 * N1(i));
    else
        dNdt = k * N1(i) * log(M1(i) / N1(i));
    end
    
    dMdt = (alpha * N1(i)) - (beta * M1(i) * surface_area);
    
    % Step forward in time
    N1(i + 1) = N1(i) + dNdt * dt;
    M1(i + 1) = M1(i) + dMdt * dt;
end

%% Run Simulation 2: Stronger treatment intensity (c = 0.2) starting at day 40
c2 = 0.2;
N2 = zeros(1, num_steps);
M2 = zeros(1, num_steps);
N2(1) = 1;
M2(1) = 1;

for i = 1:(num_steps - 1)
    if N2(i) <= 0, N2(i) = 1e-5; end
    
    surface_area = N2(i)^(2/3);
    
    % Apply standard growth before day 40, then subtract treatment kill rate (-c*N)
    if t(i) >= 40
        dNdt = k * N2(i) * log(M2(i) / N2(i)) - (c2 * N2(i));
    else
        dNdt = k * N2(i) * log(M2(i) / N2(i));
    end
    
    dMdt = (alpha * N2(i)) - (beta * M2(i) * surface_area);
    
    % Step forward in time
    N2(i + 1) = N2(i) + dNdt * dt;
    M2(i + 1) = M2(i) + dMdt * dt;
end

%% Plot and compare the tumor cell populations under both treatments
figure('Name', 'Project 2b: Anti-Tumor Treatment Comparison', 'Position', [100, 100, 900, 500]);
plot(t, N1, 'b-', 'LineWidth', 2);
hold on;
plot(t, N2, 'r--', 'LineWidth', 2);
xline(40, 'k:', 'LineWidth', 1.5, 'Label', 'Treatment Starts (Day 40)');
xlabel('Time (Days)');
ylabel('Tumor Cell Count ($N$)', 'Interpreter', 'latex');
legend('Treatment ($c = 0.1$)', 'Treatment ($c = 0.2$)', 'Location', 'NorthWest', 'Interpreter', 'latex');
title('Comparison of Tumor Growth Under Anti-Tumor Treatment');
grid on;
