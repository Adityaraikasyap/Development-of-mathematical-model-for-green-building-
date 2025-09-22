% Clear workspace and close all figures for a clean start
clc;
clear;
close all;

%% --- Load the Excel file ---
% IMPORTANT: The script will now prompt you to select the Excel file.
% This script assumes the file has Solar Radiation in column 1,
% Ambient Temperature in column 2, and Measured Indoor Temperature (MS2) in column 7.
[file, path] = uigetfile('*.xlsx', 'Select Excel file with solar radiation, ambient temp, and measured data');
if isequal(file, 0)
    error('No file selected. Please select a valid Excel file.');
end
fullpath = fullfile(path, file);

data = readtable(fullpath);
Solar_rad_all = data{:,1}; % W/m²
Ambient_temp_all = data{:,2}; % °C

% Load measured data for initial condition and comparison (MS2, Column 7)
try
    Measured_T_indoor_MS2_all = data{:,7}; % °C
    disp('Measured temperature data loaded for initial condition and comparison (MS2).');
    initial_temp = Measured_T_indoor_MS2_all(1);
catch
    warning('Could not load measured data from column 7 (MS2). Using initial ambient temperature for initial condition.');
    initial_temp = Ambient_temp_all(1);
end

% Time vector from the original data file
time_hours_raw = (0:length(Solar_rad_all)-1) * 0.8; % Assuming 0.8 hr interval
if size(time_hours_raw,1) > 1; time_hours_raw = time_hours_raw'; end

%% GLOBAL PARAMETERS (2D Adaptations)
L_x = 0.10;     % Wall thickness in x-direction [m]
L_y = 1.0;      % Wall height/width in y-direction [m] (e.g., considering a 1m x 1m section of wall)
Nx = 30;        % Number of nodes in x-direction (thickness)
Ny = 100;       % Number of nodes in y-direction (height/width)
dx = L_x / Nx;  % Node spatial step in x [m]
dy = L_y / Ny;  % Node spatial step in y [m]
dt = 1;         % Time step [s] - Kept small for explicit stability, potentially increase if stable
t_end = 72*3600; % 72 hours simulation

%% Material properties
% Brick Properties
k_brick = 0.5;      % [W/m.K]
rho_brick = 2100;   % [kg/m³]
cp_brick = 900;     % [J/kg.K]

% Mortar Properties
k_mortar = 0.7;     % [W/m.K]
rho_mortar = 2000;  % [kg/m³]
cp_mortar = 880;    % [J/kg.K]

% PCM Properties
k_pcm = 0.27;       % [W/m.K]
rho_pcm = (1086 + 899) / 2; % [kg/m³] (average of solid and liquid densities)
cp_pcm = (2120 + 2680) / 2; % [J/kg.K] (average of solid and liquid specific heats)
L_latent = 242e3;   % [J/kg] (latent heat of fusion)
T_melt_start = 28;  % [°C] (start of melting temperature range)
T_melt_end = 39;    % [°C] (end of melting temperature range)

%% Boundary & Environment Parameters
h_ext = 15;         % External convective coeff [W/m².K]
h_int = 5;         % Internal convective coeff [W/m².K]
sigma = 5.67e-8;    % Stefan-Boltzmann constant
eps_wall = 0.9;     % Emissivity of wall surface
alpha_surf = 0.6;   % Solar absorptivity

% Indoor air properties
V_air_room = 10; % Example room volume [m^3] - Adjust as per actual room size
rho_air = 1.2; cp_air = 1005; % Density and specific heat of air
m_air = V_air_room * rho_air; % Mass of air in the room
ACH = 2;            % Air Changes per Hour
m_dot_air = ACH * V_air_room * rho_air / 3600; % [kg/s] Mass flow rate of air due to ventilation
Q_internal = 0;    % Internal heat gain [W] (can be positive for heat sources, negative for heat loss/cooling)

% Calculate effective specific heat capacity for latent heat
delta_T_melt = T_melt_end - T_melt_start;
if delta_T_melt <= 0
    error('PCM melting range cannot be zero or negative. T_melt_end must be greater than T_melt_start.');
end
Cp_eff_latent = L_latent / delta_T_melt;

%% Time Discretization and Data Interpolation
t_fine = 0:dt:t_end; % Fine time vector for simulation
t_fine_hours = t_fine/3600; % Time vector in hours
Nt = length(t_fine); % Number of time steps

% Interpolate input data to match the simulation's fine time steps
Solar_rad = interp1(time_hours_raw, Solar_rad_all, t_fine_hours, 'linear', 'extrap');
T_ambient = interp1(time_hours_raw, Ambient_temp_all, t_fine_hours, 'linear', 'extrap');
T_measured_interp = interp1(time_hours_raw*3600, Measured_T_indoor_MS2_all, t_fine, 'linear', 'extrap');

% --- Stability Check for Explicit Finite Difference Scheme ---
% This is a critical check to ensure the simulation does not diverge.
alpha_min = min([k_brick/(rho_brick*cp_brick), k_mortar/(rho_mortar*cp_mortar), k_pcm/(rho_pcm*cp_pcm)]);
stability_condition = alpha_min * dt * (1/dx^2 + 1/dy^2);
if stability_condition > 0.5
    warning('Explicit finite difference scheme may be unstable. Stability condition alpha*dt*(1/dx^2 + 1/dy^2) = %.4f. Consider reducing dt.', stability_condition);
end

%% Store Results for all PCM Locations
pcm_locations = {'Exterior', 'Middle', 'Interior'};
num_locations = length(pcm_locations);

% Pre-allocate storage for results
all_sim_T_indoor_air = zeros(Nt, num_locations);
all_RME = zeros(1, num_locations);
all_RAE = zeros(1, num_locations);
all_RMSE = zeros(1, num_locations);
all_D_val_MS2 = zeros(1, num_locations);
all_SD_val_MS2 = zeros(1, num_locations);
all_LoA_lower_MS2 = zeros(1, num_locations);
all_LoA_upper_MS2 = zeros(1, num_locations);

%% Loop through each PCM location type
for loc_idx = 1:num_locations
    pcm_location_type = pcm_locations{loc_idx};
    
    % Initialize 2D temperature array for current simulation
    T_wall = ones(Ny, Nx) * initial_temp; % Temperature at each (y,x) node
    T_wall_old = T_wall; % Store old temperatures for explicit scheme

    % Define 2D Material Layers (k, rho, cp) for the wall
    % Reset maps for each run to avoid carrying over previous PCM placement
    k_map = zeros(Ny, Nx);
    rho_map = zeros(Ny, Nx);
    cp_map = zeros(Ny, Nx);

    % Base wall: Brick in the middle, Mortar layers at ends (in x-direction, thickness)
    mortar_Nx = max(1, round(Nx * 0.1)); % At least 1 node for mortar on each side
    brick_start_Nx = mortar_Nx + 1;
    brick_end_Nx = Nx - mortar_Nx;

    for j = 1:Ny % Fill base material properties across the height of the wall
        % Outer Mortar Layer
        k_map(j, 1:mortar_Nx) = k_mortar;
        rho_map(j, 1:mortar_Nx) = rho_mortar;
        cp_map(j, 1:mortar_Nx) = cp_mortar;
        % Brick Layer
        if brick_start_Nx <= brick_end_Nx
            k_map(j, brick_start_Nx:brick_end_Nx) = k_brick;
            rho_map(j, brick_start_Nx:brick_end_Nx) = rho_brick;
            cp_map(j, brick_start_Nx:brick_end_Nx) = cp_brick;
        end
        % Inner Mortar Layer
        k_map(j, brick_end_Nx+1:Nx) = k_mortar;
        rho_map(j, brick_end_Nx+1:Nx) = rho_mortar;
        cp_map(j, brick_end_Nx+1:Nx) = cp_mortar;
    end

    % Place the PCM layer
    pcm_thickness_nodes = 4; % Example: 4 nodes thick in x-direction (adjust as needed)
    if pcm_thickness_nodes > Nx
        error('PCM thickness nodes cannot exceed total wall thickness nodes (Nx).');
    end

    switch pcm_location_type
        case 'Exterior'
            pcm_x_start = 1;
            pcm_x_end = pcm_x_start + pcm_thickness_nodes - 1;
        case 'Middle'
            pcm_x_start = max(1, floor(Nx/2) - floor(pcm_thickness_nodes/2) + 1);
            pcm_x_end = min(Nx, pcm_x_start + pcm_thickness_nodes - 1);
        case 'Interior'
            pcm_x_end = Nx;
            pcm_x_start = pcm_x_end - pcm_thickness_nodes + 1;
        otherwise
            error('Invalid PCM location type specified. Choose ''Exterior'', ''Middle'', or ''Interior''.');
    end

    % Apply PCM properties to the designated region, overwriting existing material
    k_map(:, pcm_x_start:pcm_x_end) = k_pcm;
    rho_map(:, pcm_x_start:pcm_x_end) = rho_pcm;
    cp_map(:, pcm_x_start:pcm_x_end) = cp_pcm;

    % Store average indoor air temperature for this specific simulation run
    T_indoor_air_sim_current = zeros(Nt,1);
    T_indoor_air_sim_current(1) = initial_temp; % Initialize indoor air temp

    % --- Main time-stepping loop for the 2D configuration ---
    fprintf('\nStarting 2D simulation for PCM at: %s location...\n', pcm_location_type);
    for n = 1:Nt
        T_wall_old = T_wall; % Store temperatures from previous time step
        current_Solar_rad = Solar_rad(n);
        current_T_ambient = T_ambient(n);
        T_sky = current_T_ambient - 6; % Sky temperature for radiation (approximation)

        % Define current effective specific heat, including PCM latent effect (2D matrix)
        cp_eff_current = cp_map; % Start with base cp_map
        
        % Loop only through PCM nodes to apply effective specific heat
        for j_pcm = 1:Ny
            for i_pcm = pcm_x_start:pcm_x_end
                T_pcm_node_old = T_wall_old(j_pcm,i_pcm);
                % If node temperature is within the melting range, add latent heat effect
                if T_pcm_node_old > T_melt_start && T_pcm_node_old < T_melt_end
                    cp_eff_current(j_pcm,i_pcm) = cp_map(j_pcm,i_pcm) + Cp_eff_latent;
                end
            end
        end

        % --- Finite Difference Equations for 2D Wall Nodes ---
        for j = 1:Ny % Loop over y-direction (rows)
            for i = 1:Nx % Loop over x-direction (columns)
                % Material properties at current node
                k_node = k_map(j,i);
                rho_node = rho_map(j,i);
                cp_eff_node = cp_eff_current(j,i);
                
                % Control Volume Area for the current node (varies for boundary nodes)
                CV_area = dx * dy; % Default for internal node
                Q_sum = 0; % Sum of all heat transfer rates (W)

                % --- Internal Nodes (i.e., not on any boundary) ---
                if i > 1 && i < Nx && j > 1 && j < Ny
                    % Conduction from neighbors
                    Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - 2*T_wall_old(j,i) + T_wall_old(j,i-1)) + ...
                            k_node * (dx/dy) * (T_wall_old(j+1,i) - 2*T_wall_old(j,i) + T_wall_old(j-1,i));
                
                % --- Boundary Nodes ---
                % External Surface (Left Edge: i=1)
                elseif i == 1
                    % External surface heat fluxes (W/m^2)
                    q_ext_conv = h_ext * (current_T_ambient - T_wall_old(j,i));
                    q_ext_solar = alpha_surf * current_Solar_rad;
                    q_ext_rad = eps_wall * sigma * ((T_sky+273.15)^4 - (T_wall_old(j,i)+273.15)^4);
                    q_ext_total = q_ext_conv + q_ext_solar + q_ext_rad; % Total external surface flux
                    
                    if j == 1 % Top-Left Corner (i=1, j=1) - Quarter Volume
                        CV_area = (dx/2) * (dy/2);
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - T_wall_old(j,i)) + ... % From right
                                k_node * (dx/dy) * (T_wall_old(j+1,i) - T_wall_old(j,i)) + ... % From bottom
                                q_ext_total * dy; % External flux
                    elseif j == Ny % Bottom-Left Corner (i=1, j=Ny) - Quarter Volume
                        CV_area = (dx/2) * (dy/2);
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - T_wall_old(j,i)) + ... % From right
                                k_node * (dx/dy) * (T_wall_old(j-1,i) - T_wall_old(j,i)) + ... % From top
                                q_ext_total * dy; % External flux
                    else % Left Edge (not corner) (i=1, 1 < j < Ny) - Half Volume
                        CV_area = (dx/2) * dy;
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - T_wall_old(j,i)) + ... % From right
                                k_node * (dx/dy) * (T_wall_old(j+1,i) - 2*T_wall_old(j,i) + T_wall_old(j-1,i)) + ... % Y-conduction
                                q_ext_total * dy; % External flux
                    end
                
                % Internal Surface (Right Edge: i=Nx)
                elseif i == Nx
                    current_T_indoor_air = T_indoor_air_sim_current(max(1,n-1)); % Previous indoor air temp
                    q_int_conv = h_int * (current_T_indoor_air - T_wall_old(j,i));
                    q_int_total = q_int_conv; % Only convection for internal surface
                    
                    if j == 1 % Top-Right Corner (i=Nx, j=1) - Quarter Volume
                        CV_area = (dx/2) * (dy/2);
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i-1) - T_wall_old(j,i)) + ... % From left
                                k_node * (dx/dy) * (T_wall_old(j+1,i) - T_wall_old(j,i)) + ... % From bottom
                                q_int_total * dy; % Internal flux
                    elseif j == Ny % Bottom-Right Corner (i=Nx, j=Ny) - Quarter Volume
                        CV_area = (dx/2) * (dy/2);
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i-1) - T_wall_old(j,i)) + ... % From left
                                k_node * (dx/dy) * (T_wall_old(j-1,i) - T_wall_old(j,i)) + ... % From top
                                q_int_total * dy; % Internal flux
                    else % Right Edge (not corner) (i=Nx, 1 < j < Ny) - Half Volume
                        CV_area = (dx/2) * dy;
                        Q_sum = k_node * (dy/dx) * (T_wall_old(j,i-1) - T_wall_old(j,i)) + ... % From left
                                k_node * (dx/dy) * (T_wall_old(j+1,i) - 2*T_wall_old(j,i) + T_wall_old(j-1,i)) + ... % Y-conduction
                                q_int_total * dy; % Internal flux
                    end

                % Top Edge (j=1, 1 < i < Nx) - Half Volume (Adiabatic Boundary)
                elseif j == 1 && i > 1 && i < Nx
                    CV_area = dx * (dy/2);
                    Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - 2*T_wall_old(j,i) + T_wall_old(j,i-1)) + ... % X-conduction
                            k_node * (dx/dy) * (T_wall_old(j+1,i) - T_wall_old(j,i)); % From bottom
                
                % Bottom Edge (j=Ny, 1 < i < Nx) - Half Volume (Adiabatic Boundary)
                elseif j == Ny && i > 1 && i < Nx
                    CV_area = dx * (dy/2);
                    Q_sum = k_node * (dy/dx) * (T_wall_old(j,i+1) - 2*T_wall_old(j,i) + T_wall_old(j,i-1)) + ... % X-conduction
                            k_node * (dx/dy) * (T_wall_old(j-1,i) - T_wall_old(j,i)); % From top
                end
                
                % Update temperature for all nodes
                T_wall(j,i) = T_wall_old(j,i) + (dt / (rho_node * cp_eff_node * (CV_area * 1))) * Q_sum;
            end
        end

        % --- Indoor Air Temperature Update (Lumped Capacitance) ---
        % Average temperature of the inner wall surface (column Nx)
        avg_T_inner_wall_surface = mean(T_wall(:,Nx));
        T_air_old_step = T_indoor_air_sim_current(max(1,n-1)); % Previous indoor air temp
        
        % Total convective heat transfer from the entire inner wall surface to indoor air
        Q_conv_in_to_air_total = h_int * L_y * (avg_T_inner_wall_surface - T_air_old_step);
        % Heat transfer due to ventilation
        Q_ventilation = m_dot_air * cp_air * (current_T_ambient - T_air_old_step);
        % Total heat transfer to indoor air
        Q_total_to_air = Q_conv_in_to_air_total + Q_ventilation + Q_internal;
        
        % Update indoor air temperature
        T_indoor_air_sim_current(n) = T_air_old_step + (dt / (m_air * cp_air)) * Q_total_to_air;

        % Optional: Visualize 2D temperature distribution at certain intervals
        if mod(n, 3600*6) == 0 || n == Nt % Every 6 hours or at the very end
            figure('Name', sprintf('2D Wall Temp (%.1f hrs), PCM at %s', t_fine_hours(n), pcm_location_type), 'NumberTitle', 'off');
            imagesc([0 L_x], [0 L_y], T_wall); % x-range, y-range, data
            colorbar;
            axis xy; % Corrects y-axis direction (makes origin bottom-left)
            colormap jet;
            caxis([min(T_ambient)*0.9, max(T_ambient)*1.1]); % Set color scale based on ambient range
            xlabel('Wall Thickness (m)');
            ylabel('Wall Height (m)');
            title(sprintf('2D Wall Temperature (%.1f hrs), PCM at %s', t_fine_hours(n), pcm_location_type));
            drawnow; % Update plot immediately
        end
    end
    disp('2D simulation complete.');

    % Store results for current PCM location
    all_sim_T_indoor_air(:, loc_idx) = T_indoor_air_sim_current;

    % %% ---------------- Perform Error Analysis for current run ----------------
    % T_sim_clean = T_indoor_air_sim_current(:);
    % T_measured_clean = T_measured_interp(:);
    % 
    % valid_idx = ~isnan(T_sim_clean) & ~isnan(T_measured_clean);
    % T_sim_clean = T_sim_clean(valid_idx);
    % T_measured_clean = T_measured_clean(valid_idx);
    % 
    % T_measured_clean(T_measured_clean < 0.1 & T_measured_clean >= 0) = 0.1;
    % 
    % % Compute basic error metrics
    % [RME, RAE] = calculate_RME_RAE(T_sim_clean, T_measured_clean);
    % RMSE = sqrt(mean((T_sim_clean - T_measured_clean).^2));
    % 
    % % --- Bland-Altman Analysis ---
    % [D_val_MS2, SD_val_MS2, LoA_lower_MS2, LoA_upper_MS2, ~, ~, ~, ~] = ...
    %     calculate_bland_altman_metrics(T_sim_clean, T_measured_clean);
    % 
    % % Store current error metrics
    % all_RME(loc_idx) = RME;
    % all_RAE(loc_idx) = RAE;
    % all_RMSE(loc_idx) = RMSE;
    % all_D_val_MS2(loc_idx) = D_val_MS2;
    % all_SD_val_MS2(loc_idx) = SD_val_MS2;
    % all_LoA_lower_MS2(loc_idx) = LoA_lower_MS2;
    % all_LoA_upper_MS2(loc_idx) = LoA_upper_MS2;
    % 
    % % Display error metrics for the current PCM location
    % fprintf('\n--- Error Analysis: Simulated (2D Wall, PCM %s) vs Measured (MS2) ---\n', pcm_location_type);
    % fprintf('RME  = %.2f %%\n', RME);
    % fprintf('RAE  = %.2f %%\n', RAE);
    % fprintf('RMSE = %.2f °C\n', RMSE);
    % fprintf('Mean Difference (Bias, D) = %.2f °C\n', D_val_MS2);
    % fprintf('Standard Deviation of Differences (SD) = %.2f °C\n', SD_val_MS2);
    % fprintf('Limits of Agreement (LoA) = [%.2f, %.2f] °C\n', LoA_lower_MS2, LoA_upper_MS2);
    % 
    % % Plot with error metrics for the current run
    % figure('Name', sprintf('Simulated 2D vs Measured (MS2) Indoor Temperature with Error Metrics (PCM %s)', pcm_location_type), 'NumberTitle', 'off');
    % plot(t_fine_hours(valid_idx), T_sim_clean, 'b-', 'LineWidth', 1.5);
    % hold on;
    % plot(t_fine_hours(valid_idx), T_measured_clean, 'k--', 'LineWidth', 1.5);
    % xlabel('Time (hours)');
    % ylabel('Temperature (°C)');
    % title(sprintf('Simulated (2D Wall, PCM %s) vs Measured (MS2)', pcm_location_type));
    % legend(sprintf('Simulated (2D, PCM %s)', pcm_location_type), 'Measured (MS2)', 'Location', 'best');
    % grid on;
    % xlim([0 t_end/3600]);
    % 
    % % Add Error Metrics Text
    % xlim_vals = xlim;
    % ylim_vals = ylim;
    % error_string_metrics = sprintf('RME = %.2f%%\nRAE = %.2f%%\nRMSE = %.2f°C\nBias (D) = %.2f°C\nSD = %.2f°C', ...
    %                                RME, RAE, RMSE, D_val_MS2, SD_val_MS2);
    % text_x_pos = xlim_vals(1) + 0.02 * diff(xlim_vals);
    % text_y_pos = ylim_vals(2) - 0.02 * diff(ylim_vals);
    % text(text_x_pos, text_y_pos, error_string_metrics, ...
    %      'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
    %      'FontSize', 9, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    % hold off;
end

%% ########################################################################
% ## PLOT ALL SIMULATED INDOOR AIR TEMPERATURES VS MEASURED AND AMBIENT
% ########################################################################
figure('Name', 'Indoor Temperature Comparison (All 2D PCM Simulations)', 'NumberTitle', 'off');
hold on;
plot(t_fine_hours, T_measured_interp, 'k-', 'LineWidth', 2, 'DisplayName', 'Measured (MS2)'); % Measured MS2
plot(t_fine_hours, T_ambient, 'r:', 'LineWidth', 1.5, 'DisplayName', 'Ambient Temperature'); % Ambient for reference

colors = {'b', 'g', 'm'}; % Colors for the three simulated lines
for loc_idx = 1:num_locations
    plot(t_fine_hours, all_sim_T_indoor_air(:, loc_idx), ...
         'Color', colors{loc_idx}, 'LineWidth', 1.5, 'DisplayName', sprintf('Simulated (2D, PCM %s)', pcm_locations{loc_idx}));
end

hold off;
grid on;
title('Simulated 2D Indoor Air Temperature vs Measured for All PCM Locations');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
xlim([0 t_end/3600]);
set(gca, 'FontSize', 10); % Adjust font size for readability

% %% ########################################################################
% % ## DISPLAY FINAL ERROR METRICS SUMMARY TABLE
% % ########################################################################
% fprintf('\n\n--- Summary of Error Analysis for All PCM Locations ---\n');
% fprintf('%-10s %-8s %-8s %-8s %-8s %-8s %-10s %-10s\n', ...
%     'PCM Loc.', 'RME (%)', 'RAE (%)', 'RMSE (°C)', 'Bias (°C)', 'SD (°C)', 'LoA Lower', 'LoA Upper');
% fprintf('------------------------------------------------------------------------------------------------\n');
% for loc_idx = 1:num_locations
%     fprintf('%-10s %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-10.2f %-10.2f\n', ...
%         pcm_locations{loc_idx}, all_RME(loc_idx), all_RAE(loc_idx), all_RMSE(loc_idx), ...
%         all_D_val_MS2(loc_idx), all_SD_val_MS2(loc_idx), all_LoA_lower_MS2(loc_idx), all_LoA_upper_MS2(loc_idx));
% end
% fprintf('------------------------------------------------------------------------------------------------\n');
% 
% % --- Helper Function for Time Lag and Decrement Factor ---
% This function operates on 1D time series and remains unchanged.
function [time_lag_hours, decrement_factor] = deal_time_lag_decrement_factor(T_outer, T_indoor, time_data_hours, ~)
    % Find max/min and their times for outer surface
    [T_outer_max, idx_outer_max] = max(T_outer);
    [T_outer_min, ~] = min(T_outer); % Get min for amplitude calc
    time_outer_max_hours = time_data_hours(idx_outer_max);
    
    % Find max/min and their times for indoor air
    [T_indoor_max, idx_indoor_max] = max(T_indoor);
    [T_indoor_min, ~] = min(T_indoor); % Get min for amplitude calc
    time_indoor_max_hours = time_data_hours(idx_indoor_max);
    
    % Calculate time lag
    time_lag_hours = time_indoor_max_hours - time_outer_max_hours;
    cycle_duration = 24; % Assuming 24-hour cycle for adjustment
    
    % Adjust for cycle wraparound (if simulation spans multiple cycles)
    if ~isempty(time_data_hours) && (time_data_hours(end) - time_data_hours(1) > cycle_duration * 0.9)
        if time_lag_hours < -cycle_duration/2
            time_lag_hours = time_lag_hours + cycle_duration;
        elseif time_lag_hours > cycle_duration/2
            time_lag_hours = time_lag_hours - cycle_duration;
        end
    elseif time_lag_hours < 0 % If simulation is less than a full cycle but lag is negative
            time_lag_hours = time_lag_hours + cycle_duration;
    end
    if time_lag_hours < 0 % Ensure positive lag within a cycle
        time_lag_hours = mod(time_lag_hours, cycle_duration);
    end
    
    % Calculate amplitude for outer and indoor
    amplitude_outer = (T_outer_max - T_outer_min) / 2;
    amplitude_indoor = (T_indoor_max - T_indoor_min) / 2;
    
    % Calculate decrement factor
    if amplitude_outer > 1e-3 % Avoid division by near-zero amplitude
        decrement_factor = amplitude_indoor / amplitude_outer;
    else
        decrement_factor = NaN; % Amplitude too small to calculate meaningful DF
    end
end

% % --- Helper Functions for Error Analysis (No changes needed, already correct) ---
% function [RME, RAE] = calculate_RME_RAE(sim_data, meas_data)
%     sim_data = sim_data(:);
%     meas_data = meas_data(:);
%     if length(sim_data) ~= length(meas_data)
%         error('Simulated and measured data must have the same length.');
%     end
%     valid_indices = ~isnan(sim_data) & ~isnan(meas_data) & (meas_data ~= 0);
%     if sum(valid_indices) == 0
%         warning('No valid data points (non-NaN, non-zero measured) for RME/RAE calculation. RME and RAE will be NaN.');
%         RME = NaN;
%         RAE = NaN;
%         return;
%     end
%     sim_data_valid = sim_data(valid_indices);
%     meas_data_valid = meas_data(valid_indices);
%     pointwise_relative_errors = abs((sim_data_valid - meas_data_valid) ./ meas_data_valid) * 100;
%     pointwise_relative_errors = pointwise_relative_errors(isfinite(pointwise_relative_errors));
%     if isempty(pointwise_relative_errors)
%         warning('No valid finite pointwise relative errors to calculate RME/RAE after isfinite check.');
%         RME = NaN;
%         RAE = NaN;
%         return;
%     end
%     RME = max(pointwise_relative_errors);
%     RAE = mean(pointwise_relative_errors);
% end
% 
% function [D_val, SD_val, LoA_lower, LoA_upper, PCT1_val, D_max_abs_val, MNE_val, PCT2_val] = ...
%     calculate_bland_altman_metrics(sim_data, meas_data)
%     sim_data = sim_data(:);
%     meas_data = meas_data(:);
%     if length(sim_data) ~= length(meas_data)
%         error('Simulated and measured data must have the same length.');
%     end
%     valid_pair_indices = ~isnan(sim_data) & ~isnan(meas_data);
%     sim_data_valid = sim_data(valid_pair_indices);
%     meas_data_valid = meas_data(valid_pair_indices);
%     if isempty(sim_data_valid)
%         warning('No valid non-NaN data pairs for Bland-Altman calculation. All metrics will be NaN.');
%         D_val=NaN; SD_val=NaN; LoA_lower=NaN; LoA_upper=NaN; PCT1_val=NaN; D_max_abs_val=NaN; MNE_val=NaN; PCT2_val=NaN;
%         return;
%     end
%     differences = sim_data_valid - meas_data_valid;
%     D_val = mean(differences);
%     SD_val = std(differences);
%     if isnan(D_val) || isnan(SD_val) || length(differences) < 2 % std needs at least 2 elements
%         LoA_lower = NaN; LoA_upper = NaN; PCT1_val = NaN;
%     else
%         LoA_lower = D_val - 1.96 * SD_val;
%         LoA_upper = D_val + 1.96 * SD_val;
%         PCT1_val = sum(differences < LoA_lower | differences > LoA_upper) / length(differences) * 100;
%     end
%     if isempty(differences)
%         D_max_abs_val = NaN;
%     else
%         D_max_abs_val = max(abs(differences));
%     end
%     MNE_val = mean((sim_data_valid + meas_data_valid) / 2);
%     if ~isnan(MNE_val) && MNE_val ~= 0 && ~isnan(D_max_abs_val)
%         PCT2_val = (D_max_abs_val / MNE_val) * 100;
%     else
%         if ~isnan(MNE_val) && MNE_val == 0
%             warning('MNE is zero. PCT2 cannot be calculated and is set to NaN.');
%         end
%         PCT2_val = NaN;
%     end
% end