function [grid, num_abnormal, num_tumor, R_val] = updateGrid(grid, k1, k2, k3, k4, dc, phi)
    % Find grid dimensions and locate the center coordinates
    [n, ~] = size(grid);
    center = ceil(n / 2);
    
    %% 1. Calculate Tumor Metrics (N, R, and Density d)
    % Find linear indices of cells based on their current state
    abnormal_indices = find(grid == 1); % Cancerous cells
    complex_indices = find(grid == 2);  % Abnormal + effector cell complexes
    dead_indices = find(grid == 3);     % Dissolved dead cells
    
    % Count abnormal cells and total tumor mass cells (N)
    num_abnormal = length(abnormal_indices);
    tumor_indices = [abnormal_indices; complex_indices; dead_indices];
    num_tumor = length(tumor_indices);
    
    % If the tumor has been completely eradicated, return zeros
    if num_tumor == 0
        R_val = 0;
        return;
    end
    
    % Convert linear indices to 2D subscripts and calculate average radius R
    [rows, cols] = ind2sub([n, n], tumor_indices);
    radii = sqrt((rows - center).^2 + (cols - center).^2);
    R_val = sum(radii) / num_tumor;
    if R_val == 0, R_val = 1; end % Prevent division by zero errors
    
    % Calculate malignant mass density d (abnormal cells / R^2)
    d = num_abnormal / (R_val^2);
    
    % Adjust proliferation probability based on avascular carrying capacity phi
    k1_prime = k1 * (1 - (num_abnormal / phi));
    if k1_prime < 0, k1_prime = 0; end
    
    % Create a copy of the grid to handle synchronous updates across the step
    new_grid = grid;
    
    %% 2. Process Cellular Automaton Rules Cell-by-Cell
    for r = 1:n
        for c = 1:n
            state = grid(r, c);
            
            if state == 1 % --- Abnormal (Cancerous) Cell ---
                if rand() < k1_prime
                    % Cell attempts to divide; determine allowed expansion directions
                    if d > dc
                        % High mechanical pressure: allow outward growth in all 4 von Neumann directions
                        dirs = [-1, 0; 1, 0; 0, -1; 0, 1];
                    else
                        % Low mechanical pressure: restrict growth toward the center of mass (origin)
                        dirs = [];
                        if r > center, dirs = [dirs; -1, 0]; end
                        if r < center, dirs = [dirs;  1, 0]; end
                        if c > center, dirs = [dirs;  0, -1]; end
                        if c < center, dirs = [dirs;  0,  1]; end
                    end
                    
                    % Filter neighbors to find valid, in-bounds normal cells (state 0)
                    valid_targets = [];
                    for k = 1:size(dirs, 1)
                        nr = r + dirs(k, 1);
                        nc = c + dirs(k, 2);
                        if nr >= 1 && nr <= n && nc >= 1 && nc <= n
                            if grid(nr, nc) == 0
                                valid_targets = [valid_targets; nr, nc];
                            end
                        end
                    end
                    
                    % If valid space exists, pick one random normal neighbor and turn it abnormal
                    if ~isempty(valid_targets)
                        choice = randi(size(valid_targets, 1));
                        target_r = valid_targets(choice, 1);
                        target_c = valid_targets(choice, 2);
                        new_grid(target_r, target_c) = 1;
                    end
                else
                    % If the cell does not divide, immune system binds with probability k2
                    if rand() < k2
                        new_grid(r, c) = 2; % Converts into an immune-cancer complex
                    end
                end
                
            elseif state == 2 % --- Complex Cell ---
                % Complex dissolves into a dead cell with probability k3
                if rand() < k3
                    new_grid(r, c) = 3;
                end
                
            elseif state == 3 % --- Dead Cell ---
                % Healthy normal cell replaces the dead cell space with probability k4
                if rand() < k4
                    new_grid(r, c) = 0;
                end
            end
        end
    end
    
    %% 3. Post-Update Refresh
    % Apply updates to grid and recalculate final metrics for this time step
    grid = new_grid;
    num_abnormal = sum(grid(:) == 1);
    num_tumor = sum(grid(:) > 0);
    if num_tumor > 0
        [rows, cols] = find(grid > 0);
        R_val = mean(sqrt((rows - center).^2 + (cols - center).^2));
    else
        R_val = 0;
    end
end
