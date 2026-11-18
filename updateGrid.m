function [grid, num_abnormal, num_tumor, R_val] = updateGrid(grid, k1, k2, k3, k4, dc, phi)
    [n, ~] = size(grid);
    
    % Grab locations for all non-healthy cells (1 = cancer, 2 = immune complex, 3 = dead)
    abnormal_indices = find(grid == 1);
    complex_indices = find(grid == 2);
    dead_indices = find(grid == 3);
    
    num_abnormal = length(abnormal_indices);
    tumor_indices = [abnormal_indices; complex_indices; dead_indices];
    num_tumor = length(tumor_indices);
    
    % Bail out early if the immune system wiped everything out
    if num_tumor == 0
        R_val = 0;
        return;
    end
    
    % Find the tumor's current center of mass and average radius
    [rows, cols] = ind2sub([n, n], tumor_indices);
    center_r = mean(rows);
    center_c = mean(cols);
    
    radii = sqrt((rows - center_r).^2 + (cols - center_c).^2);
    R_val = sum(radii) / num_tumor;
    if R_val == 0, R_val = 1; end % avoid dividing by zero on step 1
    
    % Calculate cell density and scale down division rate as we near carrying capacity
    d = num_abnormal / (R_val^2);
    k1_prime = k1 * (1 - (num_abnormal / phi));
    if k1_prime < 0, k1_prime = 0; end
    
    % Copy the grid so updates happen all at once instead of chaining within the same step
    new_grid = grid;
    dirs = [-1, 0; 1, 0; 0, -1; 0, 1]; % up, down, left, right
    
    for r = 1:n
        for c = 1:n
            state = grid(r, c);
            
            if state == 1 % cancer cell
                if rand() < k1_prime
                    curr_dist = sqrt((r - center_r)^2 + (c - center_c)^2);
                    valid_targets = [];
                    
                    % Check the 4 adjacent spots for open healthy tissue
                    for k = 1:size(dirs, 1)
                        nr = r + dirs(k, 1);
                        nc = c + dirs(k, 2);
                        
                        if nr >= 1 && nr <= n && nc >= 1 && nc <= n && grid(nr, nc) == 0
                            if d > dc
                                % Dense enough to push outward in any direction
                                valid_targets = [valid_targets; nr, nc];
                            else
                                % Too sparse: only allow growth inward toward the tumor's center
                                new_dist = sqrt((nr - center_r)^2 + (nc - center_c)^2);
                                if new_dist < curr_dist
                                    valid_targets = [valid_targets; nr, nc];
                                end
                            end
                        end
                    end
                    
                    % Pick one open neighbor at random to divide into
                    if ~isempty(valid_targets)
                        choice = randi(size(valid_targets, 1));
                        new_grid(valid_targets(choice, 1), valid_targets(choice, 2)) = 1;
                    end
                else
                    % Failed to divide, so an immune cell can latch onto it
                    if rand() < k2
                        new_grid(r, c) = 2;
                    end
                end
                
            elseif state == 2 % immune complex breaks down into a dead cell
                if rand() < k3
                    new_grid(r, c) = 3;
                end
                
            elseif state == 3 % dead cell gets cleared and replaced by healthy tissue
                if rand() < k4
                    new_grid(r, c) = 0;
                end
            end
        end
    end
    
    % Commit the new grid and update our counts/radius for plotting
    grid = new_grid;
    num_abnormal = sum(grid(:) == 1);
    num_tumor = sum(grid(:) > 0);
    if num_tumor > 0
        [rows, cols] = find(grid > 0);
        center_r = mean(rows);
        center_c = mean(cols);
        R_val = mean(sqrt((rows - center_r).^2 + (cols - center_c).^2));
    else
        R_val = 0;
    end
end
