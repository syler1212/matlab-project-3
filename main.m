clear; close all;

%% 1. Model Parameters
% Set up the core parameters for the cellular automaton simulation from Project 1.a.
n = 101;      % Grid dimension (101x101 lattice with an odd center point)[cite: 5]
k1 = 0.74;    % Base proliferation rate for cancer cells when nutrients are unlimited[cite: 4, 5]
k2 = 0.2;     % Probability that an immune effector cell binds to a resting cancer cell[cite: 5]
k3 = 0.4;     % Probability that an immune-cancer complex breaks down into a dead cell[cite: 5]
k4 = 0.4;     % Probability that a healthy normal cell takes over a dead cell site[cite: 5]
dc = 2.85;    % Critical density threshold controlling whether pressure allows 4-way or restricted growth[cite: 4, 5]
phi = 1000;   % Carrying capacity (avascular nutrient limit where division stops)[cite: 4, 5]
steps = 300;  % Total simulation duration in time steps[cite: 5]

%% 2. Grid Initialization
% Create a blank 101x101 grid representing normal healthy tissue (state 0).
grid = zeros(n, n);
center = ceil(n / 2); % Find the exact center coordinate (row/col 51)[cite: 5]

% Plant the initial seed: a small cross of 5 abnormal (cancerous) cells in the middle.
grid(center, center) = 1;
grid(center + 1, center) = 1;
grid(center - 1, center) = 1;
grid(center, center + 1) = 1;
grid(center, center - 1) = 1;

%% 3. Data Tracking Setup
% Preallocate arrays to record population sizes and tumor size over time.
time_vec = 1:steps;
num_abnormal_vec = zeros(1, steps);
num_tumor_vec = zeros(1, steps);
R_vec = zeros(1, steps);

%% 4. Video Recording Setup
% Initialize a video writer object to capture the simulation movie for Canvas submission.
vidObj = VideoWriter('TumorGrowthSimulation.mp4', 'MPEG-4');
vidObj.FrameRate = 15; % Set playback speed
open(vidObj);

%% 5. Live Animation Figure Setup
% Create a side-by-side figure showing the 2D grid map and the live population curves.
fig = figure('Name', 'Avascular Tumor Growth Simulation', 'Position', [100, 100, 1000, 500]);

subplot(1, 2, 1);
h_im = imagesc(grid);
colormap(cool);
colorbar;
title('Tumor Growth Animation');
axis square;

subplot(1, 2, 2);
hold on;
h_abn = plot(NaN, NaN, 'r-', 'LineWidth', 1.5);
h_tum = plot(NaN, NaN, 'b-', 'LineWidth', 1.5);
xlabel('Time Steps');
ylabel('Cell Count');
legend('Abnormal Cells (N)', 'Total Tumor Cells', 'Location', 'NorthWest');
title('Tumor Population Dynamics');
grid on;

%% 6. Main Simulation Loop
% Run the cellular automaton step-by-step for the specified number of time steps.
for t = 1:steps
    % Call external function to update cell states, check pressure, and handle immune interactions
    [grid, num_abnormal, num_tumor, R_val] = updateGrid(grid, k1, k2, k3, k4, dc, phi);
    
    % Log metrics for post-simulation plotting
    num_abnormal_vec(t) = num_abnormal;
    num_tumor_vec(t) = num_tumor;
    R_vec(t) = R_val;
    
    % Refresh the visual plots in real-time
    set(h_im, 'CData', grid);
    set(h_abn, 'XData', 1:t, 'YData', num_abnormal_vec(1:t));
    set(h_tum, 'XData', 1:t, 'YData', num_tumor_vec(1:t));
    drawnow;
    
    % Capture the current frame and write it into the MP4 file
    currFrame = getframe(fig);
    writeVideo(vidObj, currFrame);
    
    % Safety check: terminate early if the tumor gets completely wiped out
    if num_abnormal == 0
        disp(['Tumor eradicated at step ', num2str(t)]);
        break;
    end
end

% Finalize and close the video file
close(vidObj);

%% 7. Radius vs. Time Plot
% Plot how the average radius of the tumor expands over the 300 steps.
figure('Name', 'Average Radius R');
plot(time_vec, R_vec, 'k-', 'LineWidth', 1.5);
xlabel('Time Steps');
ylabel('Average Radius R');
title('Average Tumor Radius vs. Time');
grid on;

%% 8. Part b: Gompertz Curve Fitting & Comparison
% Compare the simulation results against the empirical Gompertz growth model[cite: 5].
N0 = 5;          % Starting number of abnormal cells[cite: 5]
k_gompertz = 0.05; % Tuned growth constant to match early simulation trajectory (adjust as needed)[cite: 5]

% Compute the theoretical Gompertz population curve over time[cite: 5]
t_vec = 1:steps;
N_gompertz = phi * exp(log(N0 / phi) * exp(-k_gompertz * t_vec));

% Plot the side-by-side comparison graph
figure('Name', 'Part b: Simulation vs. Gompertz Function', 'Position', [150, 150, 800, 500]);
plot(t_vec, num_abnormal_vec, 'r-', 'LineWidth', 2);
hold on;
plot(t_vec, N_gompertz, 'b--', 'LineWidth', 2);
xlabel('Time Steps ($t$)', 'Interpreter', 'latex');
ylabel('Number of Abnormal Cells ($N$)', 'Interpreter', 'latex');
legend('Simulation (Part a)', sprintf('Gompertz Model (k = %.3f)', k_gompertz), 'Location', 'NorthWest');
title('Comparison of Cellular Automaton Simulation and Gompertz Function');
grid on;
