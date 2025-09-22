clc;
clear;
%% --- Load the Excel file ---
[file, path] = uigetfile('*.xlsx', 'Select Excel file with solar radiation, ambient temp, and measured data');
if isequal(file,0)
    error('No file selected');
end
fullpath = fullfile(path, file);
data = readtable(fullpath);
Solar_rad_all = data{:,1}; % W/m²
Ambient_temp_all = data{:,2}; % °C
%%% NEW/MODIFIED START (Loading Measured Data) %%%
% Assuming new columns for measured data:
% Col 3: Measured T_indoor_MS0
% Col 4: Measured T_outer_wall_MS0
% etc.
try
    Measured_T_indoor_MS0_all = data{:,3}; % °C
    Measured_T_outer_MS0_all = data{:,4}; % °C
    Measured_T_indoor_MS1_all = data{:,5}; % °C
    Measured_T_outer_MS1_all = data{:,6}; % °C
    Measured_T_indoor_MS2_all = data{:,7}; % °C
    Measured_T_outer_MS2_all = data{:,8}; % °C
    disp('Measured temperature data loaded.');
catch ME
    warning('Could not load all measured temperature data columns. Check Excel file structure and column numbers. Plots involving measured data might be incomplete or error out.');
    disp(ME.message);
end
%%% NEW/MODIFIED END %%%
% Time vector
time_hours = (0:length(Solar_rad_all)-1) * 0.8; 
if size(time_hours,1) > 1; time_hours = time_hours'; end 

%% GLOBAL PARAMETERS (Common to all simulations)
% --- Physical Dimensions ---
L = 0.1016; % Wall thickness [m]
A_wall = 1.045;   % Wall area [m²]
H = 1.0;          % **NEW**: Assumed height of the wall section [m] for 2D model
Depth = A_wall / H; % **NEW**: Calculated depth of the wall section [m]

% --- Discretization ---
Nx = 20;    % Number of nodes in x-direction (thickness)
Ny = 20;    % **NEW**: Number of nodes in y-direction (height)
dx = L / Nx;
dy = H / Ny; % **NEW**
dt = 10;    % Time step [s]
t_end = 72*3600; % 72 hours simulation time

%% Material properties (Provided by user)
k_base = 1.4;     % Thermal conductivity of base wall material (Concrete) [W/m.K]
rho_base = 2300;  % Density of base wall material (Concrete) [kg/m³]
cp_base = 950;    % Specific heat of base wall material (Concrete) [J/kg.K]
k_pcm = 0.318;    % Thermal conductivity of PCM [W/m.K]
rho_pcm = 880;    % Density of PCM [kg/m³]
cp_pcm = 2200;    % Specific heat of PCM [J/kg.K]
L_latent = 198e3; % J/kg - Latent Heat of Fusion of PCM (198 kJ/kg)
T_melt_start = 28; % °C - Start of PCM melting range (T_solidus)
T_melt_end = 39;   % °C - End of PCM melting range (T_liquidus)

% Heat Transfer Coefficients 
h_ext = 15;       % External convective coeff [W/m².K]
h_int = 5;        % Internal convective coeff [W/m².K]
sigma = 5.67e-8;  % Stefan-Boltzmann constant
eps_wall = 0.9;   % Emissivity of wall surface
alpha_surf = 0.5; % Solar absorptivity

% Indoor air properties
V_air = 0.09557;  % m³
rho_air = 1.2; cp_air = 1005;
m_air = V_air * rho_air;
ACH = 2;        % Air Changes per Hour
m_dot_air = ACH * V_air * rho_air / 3600;

% Interpolation of environmental and measured data
t_fine = 0:dt:t_end;
t_fine_hours = t_fine/3600;
Nt = length(t_fine);
Solar_rad = interp1(time_hours, Solar_rad_all, t_fine_hours, 'linear', 'extrap');
T_ambient = interp1(time_hours, Ambient_temp_all, t_fine_hours, 'linear', 'extrap');
%%% NEW/MODIFIED START (Interpolate Measured Data) %%%
Measured_T_indoor_MS0 = interp1(time_hours, Measured_T_indoor_MS0_all, t_fine_hours, 'linear', 'extrap');
Measured_T_outer_MS0 = interp1(time_hours, Measured_T_outer_MS0_all, t_fine_hours, 'linear', 'extrap');
Measured_T_indoor_MS1 = interp1(time_hours, Measured_T_indoor_MS1_all, t_fine_hours, 'linear', 'extrap');
Measured_T_outer_MS1 = interp1(time_hours, Measured_T_outer_MS1_all, t_fine_hours, 'linear', 'extrap');
Measured_T_indoor_MS2 = interp1(time_hours, Measured_T_indoor_MS2_all, t_fine_hours, 'linear', 'extrap');
Measured_T_outer_MS2 = interp1(time_hours, Measured_T_outer_MS2_all, t_fine_hours, 'linear', 'extrap');
%%% NEW/MODIFIED END %%%

% Calculate effective specific heat capacity for latent heat
delta_T_melt = T_melt_end - T_melt_start;
if delta_T_melt <= 0
    error('PCM melting range cannot be zero or negative (T_liquidus must be > T_solidus).');
end
Cp_eff_latent = L_latent / delta_T_melt;

% #########################################################################
% ## Simulation for Original Wall (MS0) - 2D
% #########################################################################
disp('Starting 2D simulation for Original Wall (MS0)...');
Q_internal_MS0 = 0;
rho_eff_MS0 = rho_base;
cp_eff_MS0 = cp_base;
k_eff_MS0 = k_base;
alpha_MS0 = k_eff_MS0 / (rho_eff_MS0 * cp_eff_MS0);
Fo_2D_MS0 = alpha_MS0 * dt * (1/dx^2 + 1/dy^2); % **NEW** 2D Stability Check
if Fo_2D_MS0 > 0.5
    warning('MS0: 2D Fourier number > 0.5. Fo = %f. Reduce time step dt for stability.', Fo_2D_MS0);
end
Fo_x_MS0 = alpha_MS0 * dt / dx^2;
Fo_y_MS0 = alpha_MS0 * dt / dy^2;

% **MODIFIED**: Initialize T as a 2D matrix
T_MS0 = ones(Nx, Ny) * Measured_T_indoor_MS0(1); 
T_air_MS0 = Measured_T_indoor_MS0(1); 
T_outer_sim_MS0 = zeros(Nt,1); 
T_indoor_sim_MS0 = zeros(Nt,1);

for n = 1:Nt
    T_sky = T_ambient(n) - 1; 
    T_MS0_old = T_MS0; 
    
    % **MODIFIED**: Loop over the 2D grid
    for i = 1:Nx
        for j = 1:Ny
            % --- Internal Nodes ---
            if i > 1 && i < Nx && j > 1 && j < Ny
                T_MS0(i,j) = T_MS0_old(i,j) ...
                    + Fo_x_MS0 * (T_MS0_old(i+1,j) - 2*T_MS0_old(i,j) + T_MS0_old(i-1,j)) ...
                    + Fo_y_MS0 * (T_MS0_old(i,j+1) - 2*T_MS0_old(i,j) + T_MS0_old(i,j-1));
            
            % --- Boundary Conditions ---
            elseif i == 1 && j > 1 && j < Ny % Outer face (left)
                q_ext_flux = alpha_surf*Solar_rad(n) + h_ext*(T_ambient(n)-T_MS0_old(1,j)) ...
                             + eps_wall*sigma*((T_sky+273.15)^4 - (T_MS0_old(1,j)+273.15)^4);
                T_MS0(1,j) = T_MS0_old(1,j) ...
                    + 2*Fo_x_MS0 * (T_MS0_old(2,j) - T_MS0_old(1,j) + q_ext_flux*dx/k_eff_MS0) ...
                    + Fo_y_MS0 * (T_MS0_old(1,j+1) - 2*T_MS0_old(1,j) + T_MS0_old(1,j-1));
            
            elseif i == Nx && j > 1 && j < Ny % Inner face (right)
                q_in_flux = h_int*(T_air_MS0 - T_MS0_old(Nx,j));
                T_MS0(Nx,j) = T_MS0_old(Nx,j) ...
                    + 2*Fo_x_MS0 * (T_MS0_old(Nx-1,j) - T_MS0_old(Nx,j) + q_in_flux*dx/k_eff_MS0) ...
                    + Fo_y_MS0 * (T_MS0_old(Nx,j+1) - 2*T_MS0_old(Nx,j) + T_MS0_old(Nx,j-1));

            elseif j == 1 && i > 1 && i < Nx % Bottom face (adiabatic)
                T_MS0(i,1) = T_MS0_old(i,1) ...
                    + Fo_x_MS0 * (T_MS0_old(i+1,1) - 2*T_MS0_old(i,1) + T_MS0_old(i-1,1)) ...
                    + 2*Fo_y_MS0 * (T_MS0_old(i,2) - T_MS0_old(i,1));

            elseif j == Ny && i > 1 && i < Nx % Top face (adiabatic)
                T_MS0(i,Ny) = T_MS0_old(i,Ny) ...
                    + Fo_x_MS0 * (T_MS0_old(i+1,Ny) - 2*T_MS0_old(i,Ny) + T_MS0_old(i-1,Ny)) ...
                    + 2*Fo_y_MS0 * (T_MS0_old(i,Ny-1) - T_MS0_old(i,Ny));
            end
        end
    end
    % Note: Corners are not explicitly calculated and will use the old temperature.
    % This is a simplification; for higher accuracy, corner-specific equations are needed.
    
    % **MODIFIED**: Update Indoor Air Temperature (sum heat from all inner nodes)
    Q_conv_in_MS0_to_air = 0;
    area_node_inner = dy * Depth;
    for j = 1:Ny
        Q_conv_in_MS0_to_air = Q_conv_in_MS0_to_air + h_int * area_node_inner * (T_MS0(Nx,j) - T_air_MS0);
    end
    Q_vent_MS0 = m_dot_air * cp_air * (T_ambient(n) - T_air_MS0); 
    Q_total_MS0 = Q_conv_in_MS0_to_air + Q_vent_MS0 + Q_internal_MS0; 
    T_air_MS0 = T_air_MS0 + dt * Q_total_MS0 / (m_air * cp_air);
    
    T_outer_sim_MS0(n) = mean(T_MS0(1,:)); % Store average outer temperature
    T_indoor_sim_MS0(n) = T_air_MS0;
end

% #########################################################################
% ## Simulation for Wall with 2.5% PCM (MS1) - 2D
% #########################################################################
disp('Starting 2D simulation for Wall with 2.5% PCM (MS1)...');
Q_internal_MS1 = 0; 
pcm_concentration_MS1 = 0.025; 
rho_eff_MS1 = (1 - pcm_concentration_MS1) * rho_base + pcm_concentration_MS1 * rho_pcm;
cp_eff_MS1_noPCM = (1 - pcm_concentration_MS1) * cp_base + pcm_concentration_MS1 * cp_pcm;
k_eff_MS1 = (1 - pcm_concentration_MS1) * k_base + pcm_concentration_MS1 * k_pcm;
alpha_MS1_check = k_eff_MS1 / (rho_eff_MS1 * cp_eff_MS1_noPCM);
Fo_2D_MS1 = alpha_MS1_check * dt * (1/dx^2 + 1/dy^2);
if Fo_2D_MS1 > 0.5
    warning('MS1: 2D Fourier number > 0.5. Fo = %f. Reduce dt for stability.', Fo_2D_MS1);
end
T_MS1 = ones(Nx,Ny) * Measured_T_indoor_MS1(1); 
T_air_MS1 = Measured_T_indoor_MS1(1); 
T_outer_sim_MS1 = zeros(Nt,1); 
T_indoor_sim_MS1 = zeros(Nt,1);

for n = 1:Nt
    T_sky = T_ambient(n) - 10;
    T_MS1_old = T_MS1;

    % Create matrices for node-specific properties
    cp_eff_matrix = ones(Nx,Ny) * cp_eff_MS1_noPCM;
    pcm_mask = (T_MS1_old > T_melt_start) & (T_MS1_old < T_melt_end);
    cp_eff_matrix(pcm_mask) = cp_eff_matrix(pcm_mask) + pcm_concentration_MS1 * Cp_eff_latent;
    alpha_matrix = k_eff_MS1 ./ (rho_eff_MS1 * cp_eff_matrix);
    Fo_x_matrix = alpha_matrix * dt / dx^2;
    Fo_y_matrix = alpha_matrix * dt / dy^2;

    for i = 1:Nx
        for j = 1:Ny
            Fo_x = Fo_x_matrix(i,j);
            Fo_y = Fo_y_matrix(i,j);
            % --- Internal Nodes ---
            if i > 1 && i < Nx && j > 1 && j < Ny
                T_MS1(i,j) = T_MS1_old(i,j) ...
                    + Fo_x * (T_MS1_old(i+1,j) - 2*T_MS1_old(i,j) + T_MS1_old(i-1,j)) ...
                    + Fo_y * (T_MS1_old(i,j+1) - 2*T_MS1_old(i,j) + T_MS1_old(i,j-1));
            
            % --- Boundary Conditions ---
            elseif i == 1 && j > 1 && j < Ny % Outer face
                q_ext_flux = alpha_surf*Solar_rad(n) + h_ext*(T_ambient(n)-T_MS1_old(1,j)) ...
                             + eps_wall*sigma*((T_sky+273.15)^4 - (T_MS1_old(1,j)+273.15)^4);
                T_MS1(1,j) = T_MS1_old(1,j) ...
                    + 2*Fo_x * (T_MS1_old(2,j) - T_MS1_old(1,j) + q_ext_flux*dx/k_eff_MS1) ...
                    + Fo_y * (T_MS1_old(1,j+1) - 2*T_MS1_old(1,j) + T_MS1_old(1,j-1));

            elseif i == Nx && j > 1 && j < Ny % Inner face
                q_in_flux = h_int*(T_air_MS1 - T_MS1_old(Nx,j));
                T_MS1(Nx,j) = T_MS1_old(Nx,j) ...
                    + 2*Fo_x * (T_MS1_old(Nx-1,j) - T_MS1_old(Nx,j) + q_in_flux*dx/k_eff_MS1) ...
                    + Fo_y * (T_MS1_old(Nx,j+1) - 2*T_MS1_old(Nx,j) + T_MS1_old(Nx,j-1));

            elseif j == 1 && i > 1 && i < Nx % Bottom face
                T_MS1(i,1) = T_MS1_old(i,1) ...
                    + Fo_x * (T_MS1_old(i+1,1) - 2*T_MS1_old(i,1) + T_MS1_old(i-1,1)) ...
                    + 2*Fo_y * (T_MS1_old(i,2) - T_MS1_old(i,1));

            elseif j == Ny && i > 1 && i < Nx % Top face
                T_MS1(i,Ny) = T_MS1_old(i,Ny) ...
                    + Fo_x * (T_MS1_old(i+1,Ny) - 2*T_MS1_old(i,Ny) + T_MS1_old(i-1,Ny)) ...
                    + 2*Fo_y * (T_MS1_old(i,Ny-1) - T_MS1_old(i,Ny));
            end
        end
    end
    
    Q_conv_in_MS1_to_air = 0;
    area_node_inner = dy * Depth;
    for j = 1:Ny
        Q_conv_in_MS1_to_air = Q_conv_in_MS1_to_air + h_int * area_node_inner * (T_MS1(Nx,j) - T_air_MS1);
    end
    Q_vent_MS1 = m_dot_air * cp_air * (T_ambient(n) - T_air_MS1);
    Q_total_MS1 = Q_conv_in_MS1_to_air + Q_vent_MS1 + Q_internal_MS1;
    T_air_MS1 = T_air_MS1 + dt * Q_total_MS1 / (m_air * cp_air);
    
    T_outer_sim_MS1(n) = mean(T_MS1(1,:));
    T_indoor_sim_MS1(n) = T_air_MS1;
end

% #########################################################################
% ## Simulation for Wall with 4% PCM (MS2) - 2D
% #########################################################################
disp('Starting 2D simulation for Wall with 4% PCM (MS2)...');
Q_internal_MS2 = 0; 
pcm_concentration_MS2 = 0.04; 
rho_eff_MS2 = (1 - pcm_concentration_MS2) * rho_base + pcm_concentration_MS2 * rho_pcm;
cp_eff_MS2_noPCM = (1 - pcm_concentration_MS2) * cp_base + pcm_concentration_MS2 * cp_pcm;
k_eff_MS2 = (1 - pcm_concentration_MS2) * k_base + pcm_concentration_MS2 * k_pcm;
alpha_MS2_check = k_eff_MS2 / (rho_eff_MS2 * cp_eff_MS2_noPCM); 
Fo_2D_MS2 = alpha_MS2_check * dt * (1/dx^2 + 1/dy^2);
if Fo_2D_MS2 > 0.5
    warning('MS2: 2D Fourier number > 0.5. Fo = %f. Reduce dt for stability.', Fo_2D_MS2);
end
T_MS2 = ones(Nx,Ny) * Measured_T_indoor_MS2(1); 
T_air_MS2 = Measured_T_indoor_MS2(1); 
T_outer_sim_MS2 = zeros(Nt,1); 
T_indoor_sim_MS2 = zeros(Nt,1);

for n = 1:Nt
    T_sky = T_ambient(n) - 15;
    T_MS2_old = T_MS2;
    
    cp_eff_matrix = ones(Nx,Ny) * cp_eff_MS2_noPCM;
    pcm_mask = (T_MS2_old > T_melt_start) & (T_MS2_old < T_melt_end);
    cp_eff_matrix(pcm_mask) = cp_eff_matrix(pcm_mask) + pcm_concentration_MS2 * Cp_eff_latent;
    alpha_matrix = k_eff_MS2 ./ (rho_eff_MS2 * cp_eff_matrix);
    Fo_x_matrix = alpha_matrix * dt / dx^2;
    Fo_y_matrix = alpha_matrix * dt / dy^2;

    for i = 1:Nx
        for j = 1:Ny
            Fo_x = Fo_x_matrix(i,j);
            Fo_y = Fo_y_matrix(i,j);
            % --- Internal Nodes ---
            if i > 1 && i < Nx && j > 1 && j < Ny
                T_MS2(i,j) = T_MS2_old(i,j) ...
                    + Fo_x * (T_MS2_old(i+1,j) - 2*T_MS2_old(i,j) + T_MS2_old(i-1,j)) ...
                    + Fo_y * (T_MS2_old(i,j+1) - 2*T_MS2_old(i,j) + T_MS2_old(i,j-1));
            
            % --- Boundary Conditions ---
            elseif i == 1 && j > 1 && j < Ny % Outer face
                q_ext_flux = alpha_surf*Solar_rad(n) + h_ext*(T_ambient(n)-T_MS2_old(1,j)) ...
                             + eps_wall*sigma*((T_sky+273.15)^4 - (T_MS2_old(1,j)+273.15)^4);
                T_MS2(1,j) = T_MS2_old(1,j) ...
                    + 2*Fo_x * (T_MS2_old(2,j) - T_MS2_old(1,j) + q_ext_flux*dx/k_eff_MS2) ...
                    + Fo_y * (T_MS2_old(1,j+1) - 2*T_MS2_old(1,j) + T_MS2_old(1,j-1));

            elseif i == Nx && j > 1 && j < Ny % Inner face
                q_in_flux = h_int*(T_air_MS2 - T_MS2_old(Nx,j));
                T_MS2(Nx,j) = T_MS2_old(Nx,j) ...
                    + 2*Fo_x * (T_MS2_old(Nx-1,j) - T_MS2_old(Nx,j) + q_in_flux*dx/k_eff_MS2) ...
                    + Fo_y * (T_MS2_old(Nx,j+1) - 2*T_MS2_old(Nx,j) + T_MS2_old(Nx,j-1));

            elseif j == 1 && i > 1 && i < Nx % Bottom face
                T_MS2(i,1) = T_MS2_old(i,1) ...
                    + Fo_x * (T_MS2_old(i+1,1) - 2*T_MS2_old(i,1) + T_MS2_old(i-1,1)) ...
                    + 2*Fo_y * (T_MS2_old(i,2) - T_MS2_old(i,1));

            elseif j == Ny && i > 1 && i < Nx % Top face
                T_MS2(i,Ny) = T_MS2_old(i,Ny) ...
                    + Fo_x * (T_MS2_old(i+1,Ny) - 2*T_MS2_old(i,Ny) + T_MS2_old(i-1,Ny)) ...
                    + 2*Fo_y * (T_MS2_old(i,Ny-1) - T_MS2_old(i,Ny));
            end
        end
    end
    
    Q_conv_in_MS2_to_air = 0;
    area_node_inner = dy * Depth;
    for j = 1:Ny
        Q_conv_in_MS2_to_air = Q_conv_in_MS2_to_air + h_int * area_node_inner * (T_MS2(Nx,j) - T_air_MS2);
    end
    Q_vent_MS2 = m_dot_air * cp_air * (T_ambient(n) - T_air_MS2);
    Q_total_MS2 = Q_conv_in_MS2_to_air + Q_vent_MS2 + Q_internal_MS2;
    T_air_MS2 = T_air_MS2 + dt * Q_total_MS2 / (m_air * cp_air);
    
    T_outer_sim_MS2(n) = mean(T_MS2(1,:));
    T_indoor_sim_MS2(n) = T_air_MS2;
end
disp('All 2D simulations complete.');

% You can add your plotting code here to visualize the results.
% The output vectors (T_outer_sim_MS*, T_indoor_sim_MS*) can be used
% directly with your existing plotting scripts.
% ## Error Analysis (Simulated vs. Measured Indoor Temperatures)
fprintf('\n--- Error Analysis: Simulated vs. Measured Indoor Temperatures ---\n');

% MS0 Error Analysis
[RME_MS0, RAE_MS0] = calculate_RME_RAE(T_indoor_sim_MS0, Measured_T_indoor_MS0);
[D_MS0, SD_MS0, LoA_lower_MS0, LoA_upper_MS0, PCT1_MS0, D_max_MS0, MNE_MS0, PCT2_MS0] = ...
    calculate_bland_altman_metrics(T_indoor_sim_MS0, Measured_T_indoor_MS0);
fprintf('MS0 (Original Wall):\n');
fprintf('  RME: %.2f %%, RAE: %.2f %%\n', RME_MS0, RAE_MS0);
fprintf('  Bland-Altman: D=%.2f°C, SD=%.2f°C, LoA=[%.2f, %.2f]°C\n', D_MS0, SD_MS0, LoA_lower_MS0, LoA_upper_MS0);
fprintf('                PCT1=%.2f%%, D_max=%.2f°C, MNE=%.2f°C, PCT2=%.2f%%\n', PCT1_MS0, D_max_MS0, MNE_MS0, PCT2_MS0);

% MS1 Error Analysis
[RME_MS1, RAE_MS1] = calculate_RME_RAE(T_indoor_sim_MS1, Measured_T_indoor_MS1);
[D_MS1, SD_MS1, LoA_lower_MS1, LoA_upper_MS1, PCT1_MS1, D_max_MS1, MNE_MS1, PCT2_MS1] = ...
    calculate_bland_altman_metrics(T_indoor_sim_MS1, Measured_T_indoor_MS1);
fprintf('MS1 (2.5%% PCM Wall):\n');
fprintf('  RME: %.2f %%, RAE: %.2f %%\n', RME_MS1, RAE_MS1);
fprintf('  Bland-Altman: D=%.2f°C, SD=%.2f°C, LoA=[%.2f, %.2f]°C\n', D_MS1, SD_MS1, LoA_lower_MS1, LoA_upper_MS1);
fprintf('                PCT1=%.2f%%, D_max=%.2f°C, MNE=%.2f°C, PCT2=%.2f%%\n', PCT1_MS1, D_max_MS1, MNE_MS1, PCT2_MS1);

% MS2 Error Analysis
[RME_MS2, RAE_MS2] = calculate_RME_RAE(T_indoor_sim_MS2, Measured_T_indoor_MS2);
[D_MS2, SD_MS2, LoA_lower_MS2, LoA_upper_MS2, PCT1_MS2, D_max_MS2, MNE_MS2, PCT2_MS2] = ...
    calculate_bland_altman_metrics(T_indoor_sim_MS2, Measured_T_indoor_MS2);
fprintf('MS2 (4%% PCM Wall):\n');
fprintf('  RME: %.2f %%, RAE: %.2f %%\n', RME_MS2, RAE_MS2);
fprintf('  Bland-Altman: D=%.2f°C, SD=%.2f°C, LoA=[%.2f, %.2f]°C\n', D_MS2, SD_MS2, LoA_lower_MS2, LoA_upper_MS2);
fprintf('                PCT1=%.2f%%, D_max=%.2f°C, MNE=%.2f°C, PCT2=%.2f%%\n', PCT1_MS2, D_max_MS2, MNE_MS2, PCT2_MS2);
fprintf('---------------------------------------------------------------------\n');


% ## Calculation of Time Lag and Decrement Factor
sim_period_seconds = t_fine(end) - t_fine(end-round(24*3600/dt) + 1);
if sim_period_seconds < 23.9*3600
    warning('Simulation period might be too short for accurate time lag/decrement factor calculation.');
end
start_idx_24hr = Nt - round(24*3600/dt) + 1;
if start_idx_24hr < 1; start_idx_24hr = 1; end 
end_idx_24hr = Nt;

t_24hr_hours = t_fine_hours(start_idx_24hr:end_idx_24hr);
T_outer_MS0_24hr = T_outer_sim_MS0(start_idx_24hr:end_idx_24hr);
T_indoor_MS0_24hr = T_indoor_sim_MS0(start_idx_24hr:end_idx_24hr);
T_outer_MS1_24hr = T_outer_sim_MS1(start_idx_24hr:end_idx_24hr);
T_indoor_MS1_24hr = T_indoor_sim_MS1(start_idx_24hr:end_idx_24hr);
T_outer_MS2_24hr = T_outer_sim_MS2(start_idx_24hr:end_idx_24hr);
T_indoor_MS2_24hr = T_indoor_sim_MS2(start_idx_24hr:end_idx_24hr);

calculate_thermal_factors = @(T_outer_data, T_indoor_data, time_data) deal_time_lag_decrement_factor(T_outer_data, T_indoor_data, time_data, dt/3600);

[TL_MS0, DF_MS0] = calculate_thermal_factors(T_outer_MS0_24hr, T_indoor_MS0_24hr, t_24hr_hours);
[TL_MS1, DF_MS1] = calculate_thermal_factors(T_outer_MS1_24hr, T_indoor_MS1_24hr, t_24hr_hours);
[TL_MS2, DF_MS2] = calculate_thermal_factors(T_outer_MS2_24hr, T_indoor_MS2_24hr, t_24hr_hours);

fprintf('\n--- Thermal Performance Metrics (Simulated - Last 24h) ---\n');
fprintf('Scenario       | Time Lag (hours) | Decrement Factor\n');
fprintf('--------------------------------------------------\n');
fprintf('MS0 (Original) | %16.2f | %16.3f\n', TL_MS0, DF_MS0);
fprintf('MS1 (2.5%% PCM)   | %16.2f | %16.3f\n', TL_MS1, DF_MS1); 
fprintf('MS2 (4%% PCM) | %16.2f | %16.3f\n', TL_MS2, DF_MS2); 
fprintf('--------------------------------------------------\n');

% ## Plots 
%%% NEW/MODIFIED START (Individual Scenario Plots with Measured Data and Full Comparison) %%%

% Plot 1 (Focused): Indoor Temperature for MS0 (Original Wall) - Simulated vs Measured with Error Metrics
figure;
plot(t_fine_hours, T_indoor_sim_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'MS0 Indoor Temp (Simulated)'); hold on;
% if ~all(isnan(Measured_T_indoor_MS0))
%     plot(t_fine_hours, Measured_T_indoor_MS0, 'b--', 'LineWidth', 1.5, 'DisplayName', 'MS0 Indoor Temp (Measured)');
% else
%     warning('Measured Indoor Temp (MS0) is all NaN, skipping plot for focused indoor temp.');
% end
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS0');
grid on;
set(gca, 'FontSize', 12);
% Add error metrics text
error_string_MS0 = sprintf(['Error Metrics (MS0 Indoor):\n', ...
    'RME: %.2f %%, RAE: %.2f %%\n', ...
    'D: %.2f °C, SD: %.2f °C\n', ...
    'LoA: [%.2f, %.2f] °C\n', ...
    'PCT1: %.2f %%\n', ...
    'D_max: %.2f °C, MNE: %.2f °C\n', ...
    'PCT2: %.2f %%'], ...
    RME_MS0, RAE_MS0, D_MS0, SD_MS0, LoA_lower_MS0, LoA_upper_MS0, PCT1_MS0, D_max_MS0, MNE_MS0, PCT2_MS0);
xlim_vals = xlim; ylim_vals = ylim;
text(xlim_vals(1) + 0.02*diff(xlim_vals), ylim_vals(2) - 0.02*diff(ylim_vals), error_string_MS0, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Plot 2 (Focused): Indoor Temperature for MS1 (2.5% PCM) - Simulated vs Measured with Error Metrics
figure;
plot(t_fine_hours, T_indoor_sim_MS1, 'g', 'LineWidth', 1.5, 'DisplayName', 'MS1 Indoor Temp (Simulated)'); hold on;
% if ~all(isnan(Measured_T_indoor_MS1))
%     plot(t_fine_hours, Measured_T_indoor_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'MS1 Indoor Temp (Measured)');
% else
%     warning('Measured Indoor Temp (MS1) is all NaN, skipping plot for focused indoor temp.');
% end
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS1');
grid on;
set(gca, 'FontSize', 12);
% % Add error metrics text
error_string_MS1 = sprintf(['Error Metrics (MS1 Indoor):\n', ...
    'RME: %.2f %%, RAE: %.2f %%\n', ...
    'D: %.2f °C, SD: %.2f °C\n', ...
    'LoA: [%.2f, %.2f] °C\n', ...
    'PCT1: %.2f %%\n', ...
    'D_max: %.2f °C, MNE: %.2f °C\n', ...
    'PCT2: %.2f %%'], ...
    RME_MS1, RAE_MS1, D_MS1, SD_MS1, LoA_lower_MS1, LoA_upper_MS1, PCT1_MS1, D_max_MS1, MNE_MS1, PCT2_MS1);
xlim_vals = xlim; ylim_vals = ylim;
text(xlim_vals(1) + 0.02*diff(xlim_vals), ylim_vals(2) - 0.02*diff(ylim_vals), error_string_MS1, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Plot 3 (Focused): Indoor Temperature for MS2 (4% PCM) - Simulated vs Measured with Error Metrics
figure;
plot(t_fine_hours, T_indoor_sim_MS2, 'r', 'LineWidth', 1.5, 'DisplayName', 'MS2 Indoor Temp (Simulated)'); hold on;
% if ~all(isnan(Measured_T_indoor_MS2))
%     plot(t_fine_hours, Measured_T_indoor_MS2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'MS2 Indoor Temp (Measured)');
% else
%     warning('Measured Indoor Temp (MS2) is all NaN, skipping plot.');
% end
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS2');
grid on;
set(gca, 'FontSize', 12);
% Add error metrics text
error_string_MS2 = sprintf(['Error Metrics (MS2 Indoor):\n', ...
    'RME: %.2f %%, RAE: %.2f %%\n', ...
    'D: %.2f °C, SD: %.2f °C\n', ...
    'LoA: [%.2f, %.2f] °C\n', ...
    'PCT1: %.2f %%\n', ...
    'D_max: %.2f °C, MNE: %.2f °C\n', ...
    'PCT2: %.2f %%'], ...
    RME_MS2, RAE_MS2, D_MS2, SD_MS2, LoA_lower_MS2, LoA_upper_MS2, PCT1_MS2, D_max_MS2, MNE_MS2, PCT2_MS2);
xlim_vals = xlim; ylim_vals = ylim;
text(xlim_vals(1) + 0.02*diff(xlim_vals), ylim_vals(2) - 0.02*diff(ylim_vals), error_string_MS2, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', 'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
hold off;

% Plot 1: Original Wall (MS0 - Reference) - Modelled vs Measured (Full Detail)
figure;
plot(t_fine_hours, T_indoor_sim_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Simulated)'); hold on;
plot(t_fine_hours, Measured_T_indoor_MS0, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Measured)');
% plot(t_fine_hours, T_outer_sim_MS0, 'k', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Simulated)');
% plot(t_fine_hours, Measured_T_outer_MS0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Measured)');
% plot(t_fine_hours, T_ambient, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Ambient Temp');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Thermal Performance: Original Wall (MS0) - Simulated vs Measured');
grid on;
set(gca, 'FontSize', 12);
hold off;

% Plot 2: Wall with 4% PCM (MS1) - Modelled vs Measured (Full Detail)
figure;
plot(t_fine_hours, T_indoor_sim_MS1, 'b', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Simulated)'); hold on;
plot(t_fine_hours, Measured_T_indoor_MS1, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Measured)');
% plot(t_fine_hours, T_outer_sim_MS1, 'k', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Simulated)');
% plot(t_fine_hours, Measured_T_outer_MS1, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Measured)');
% plot(t_fine_hours, T_ambient, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Ambient Temp');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Thermal Performance: Wall with 4% PCM (MS1) - Simulated vs Measured');
grid on;
set(gca, 'FontSize', 12);
hold off;

% Plot 3: Wall with 2.5% PCM (MS2) - Modelled vs Measured (Full Detail)
figure;
plot(t_fine_hours, T_indoor_sim_MS2, 'b', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Simulated)'); hold on;
plot(t_fine_hours, Measured_T_indoor_MS2, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Indoor Air Temp (Measured)');
% plot(t_fine_hours, T_outer_sim_MS2, 'k', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Simulated)');
% plot(t_fine_hours, Measured_T_outer_MS2, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Temp (Measured)');
% plot(t_fine_hours, T_ambient, 'r-.', 'LineWidth', 1.5, 'DisplayName', 'Ambient Temp');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Thermal Performance: Wall with 2.5% PCM (MS2) - Simulated vs Measured');
grid on;
set(gca, 'FontSize', 12);
hold off;

% Create a single, wide figure to hold all subplots horizontally
figure('Position', [100, 100, 1600, 500]); % [left, bottom, width, height]

% --- Plot 1: MS0 (Original Wall) ---
subplot(1, 3, 1);
plot(t_fine_hours, T_indoor_sim_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'MS0 Indoor Temp (Simulated)');
hold on;
plot(t_fine_hours, Measured_T_indoor_MS0, 'b--', 'LineWidth', 1.5, 'DisplayName', 'MS0 Indoor Temp (Measured)');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS0');
grid on;
set(gca, 'FontSize', 12);

% --- CUSTOMIZATIONS FOR MS0 ---
% 1. Set x-axis limit to 72 hours
xlim([0, 72]);
hold off; % Briefly turn off hold to allow MATLAB to recalculate axes limits

% 2. Place error metrics just above the center of the curve
ylim_vals = ylim; % Get the final y-axis limits
hold on; % Turn hold back on for adding text

% Define the error metrics string
error_string_MS0 = sprintf(['RME: %.2f %% | RAE: %.2f %%\n' ...
    'D: %.2f °C | SD: %.2f °C\n' ...
    'LoA: [%.2f, %.2f] °C'], ...
    RME_MS0, RAE_MS0, D_MS0, SD_MS0, LoA_lower_MS0, LoA_upper_MS0);

% Calculate the position for the text box
x_text_pos = mean(xlim); % Center of the x-axis
% Find the height of the simulated curve at this central point
y_curve_val = interp1(t_fine_hours, T_indoor_sim_MS0, x_text_pos);
% Position the text slightly above this point
y_text_pos = y_curve_val + 0.05 * diff(ylim_vals);

% Add the text to the plot
text(x_text_pos, y_text_pos, error_string_MS0, ...
    'HorizontalAlignment', 'center', ...  % Center the box horizontally
    'VerticalAlignment', 'bottom', ...    % Position the box *above* the y-coordinate
    'FontSize', 9, ...
    'BackgroundColor', [1 1 0.9], ... % Light yellow background
    'EdgeColor', 'k');
hold off;

% --- Plot 2: MS1 (Wall with 2.5% PCM) ---
subplot(1, 3, 2);
plot(t_fine_hours, T_indoor_sim_MS1, 'g', 'LineWidth', 1.5, 'DisplayName', 'MS1 Indoor Temp (Simulated)');
hold on;
plot(t_fine_hours, Measured_T_indoor_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'MS1 Indoor Temp (Measured)');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS1');
grid on;
set(gca, 'FontSize', 12);

% --- CUSTOMIZATIONS FOR MS1 ---
xlim([0, 72]);
hold off;
ylim_vals = ylim;
hold on;

error_string_MS1 = sprintf(['RME: %.2f %% | RAE: %.2f %%\n' ...
    'D: %.2f °C | SD: %.2f °C\n' ...
    'LoA: [%.2f, %.2f] °C'], ...
    RME_MS1, RAE_MS1, D_MS1, SD_MS1, LoA_lower_MS1, LoA_upper_MS1);

x_text_pos = mean(xlim);
y_curve_val = interp1(t_fine_hours, T_indoor_sim_MS1, x_text_pos);
y_text_pos = y_curve_val + 0.05 * diff(ylim_vals);

text(x_text_pos, y_text_pos, error_string_MS1, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 9, ...
    'BackgroundColor', [0.9 1 0.9], ... % Light green background
    'EdgeColor', 'k');
hold off;

% --- Plot 3: MS2 (Wall with 4% PCM) ---
subplot(1, 3, 3);
plot(t_fine_hours, T_indoor_sim_MS2, 'r', 'LineWidth', 1.5, 'DisplayName', 'MS2 Indoor Temp (Simulated)');
hold on;
plot(t_fine_hours, Measured_T_indoor_MS2, 'r--', 'LineWidth', 1.5, 'DisplayName', 'MS2 Indoor Temp (Measured)');
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('MS2');
grid on;
set(gca, 'FontSize', 12);

% --- CUSTOMIZATIONS FOR MS2 ---
xlim([0, 72]);
hold off;
ylim_vals = ylim;
hold on;

error_string_MS2 = sprintf(['RME: %.2f %% | RAE: %.2f %%\n' ...
    'D: %.2f °C | SD: %.2f °C\n' ...
    'LoA: [%.2f, %.2f] °C'], ...
    RME_MS2, RAE_MS2, D_MS2, SD_MS2, LoA_lower_MS2, LoA_upper_MS2);

x_text_pos = mean(xlim);
y_curve_val = interp1(t_fine_hours, T_indoor_sim_MS2, x_text_pos);
y_text_pos = y_curve_val + 0.05 * diff(ylim_vals);

text(x_text_pos, y_text_pos, error_string_MS2, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', ...
    'FontSize', 9, ...
    'BackgroundColor', [1 0.9 0.9], ... % Light red background
    'EdgeColor', 'k');
hold off;
% %%% NEW/MODIFIED END %%%


% ## Comparison Plots
% ---
% Plot 4: Comparison of Indoor Air Temperatures for all scenarios (Simulated and Measured)
figure;
subplot(2,1,1); 
plot(t_fine_hours, T_indoor_sim_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'Indoor Sim (MS0 - Original)'); hold on;
plot(t_fine_hours, T_indoor_sim_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Indoor Sim (MS1 - 4% PCM)'); 
plot(t_fine_hours, T_indoor_sim_MS2, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Indoor Sim (MS2 - 2.5% PCM)'); 
plot(t_fine_hours, T_ambient, 'r', 'LineWidth', 1, 'DisplayName', 'Ambient Temp');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Comparison of Simulated Indoor Air Temperatures');
grid on;
set(gca, 'FontSize', 10);
hold off;

subplot(2,1,2); 
plot(t_fine_hours, Measured_T_indoor_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'Indoor Meas (MS0 - Original)'); hold on;
plot(t_fine_hours, Measured_T_indoor_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Indoor Meas (MS1 - 4% PCM)'); 
plot(t_fine_hours, Measured_T_indoor_MS2, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Indoor Meas (MS2 - 2.5% PCM)'); 
plot(t_fine_hours, T_ambient, 'r', 'LineWidth', 1, 'DisplayName', 'Ambient Temp'); 
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Comparison of Measured Indoor Air Temperatures');
grid on;
set(gca, 'FontSize', 10);
hold off;
sgtitle('Indoor Air Temperature Comparisons: Simulated vs. Measured','FontSize',14); 

% Plot 5: Comparison of Outer Wall Surface Temperatures for all scenarios (Simulated and Measured)
figure;
subplot(2,1,1); 
plot(t_fine_hours, T_outer_sim_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Sim (MS0 - Original)'); hold on;
plot(t_fine_hours, T_outer_sim_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Sim (MS1 - 4% PCM)'); 
plot(t_fine_hours, T_outer_sim_MS2, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Sim (MS2 - 2.5% PCM)'); 
plot(t_fine_hours, T_ambient, 'r', 'LineWidth', 1, 'DisplayName', 'Ambient Temp'); 
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Comparison of Simulated Outer Wall Surface Temperatures');
grid on;
set(gca, 'FontSize', 10);
hold off;

subplot(2,1,2); 
plot(t_fine_hours, Measured_T_outer_MS0, 'b', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Meas (MS0 - Original)'); hold on;
plot(t_fine_hours, Measured_T_outer_MS1, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Meas (MS1 - 4% PCM)'); 
plot(t_fine_hours, Measured_T_outer_MS2, 'm:', 'LineWidth', 1.5, 'DisplayName', 'Outer Wall Meas (MS2 - 2.5% PCM)'); 
plot(t_fine_hours, T_ambient, 'r', 'LineWidth', 1, 'DisplayName', 'Ambient Temp'); 
xlabel('Time (hours)');
ylabel('Temperature (°C)');
legend('Location', 'best');
title('Comparison of Measured Outer Wall Surface Temperatures');
grid on;
set(gca, 'FontSize', 10);
hold off;
sgtitle('Outer Wall Surface Temperature Comparisons: Simulated vs. Measured','FontSize',14); 

%% --- Helper Function for Time Lag and Decrement Factor ---
function [time_lag_hours, decrement_factor] = deal_time_lag_decrement_factor(T_outer, T_indoor, time_data_hours, ~) % dt_hours not used
    % Find max/min and their times for outer surface
    [T_outer_max, idx_outer_max] = max(T_outer);
    [T_outer_min, idx_outer_min] = min(T_outer);
    time_outer_max_hours = time_data_hours(idx_outer_max);

    % Find max/min and their times for indoor air
    [T_indoor_max, idx_indoor_max] = max(T_indoor);
    [T_indoor_min, idx_indoor_min] = min(T_indoor);
    time_indoor_max_hours = time_data_hours(idx_indoor_max);

    % Calculate time lag
    time_lag_hours = time_indoor_max_hours - time_outer_max_hours;

    cycle_duration = 24; % Assuming 24-hour cycle for adjustment
    if ~isempty(time_data_hours) && (time_data_hours(end) - time_data_hours(1) > cycle_duration * 0.9) 
        if time_lag_hours < -cycle_duration/2 
            time_lag_hours = time_lag_hours + cycle_duration;
        elseif time_lag_hours > cycle_duration/2 
            time_lag_hours = time_lag_hours - cycle_duration;
        end
    elseif time_lag_hours < 0 
         time_lag_hours = time_lag_hours + cycle_duration; 
    end
    if time_lag_hours < 0 
        time_lag_hours = mod(time_lag_hours, cycle_duration); 
    end

    % Calculate amplitude for outer and indoor
    amplitude_outer = (T_outer_max - T_outer_min) / 2;
    amplitude_indoor = (T_indoor_max - T_indoor_min) / 2;

    % Calculate decrement factor
    if amplitude_outer > 1e-3 
        decrement_factor = amplitude_indoor / amplitude_outer;
    else
        decrement_factor = NaN; 
    end
end

%% --- Helper Functions for Error Analysis ---
function [RME, RAE] = calculate_RME_RAE(sim_data, meas_data)
    % Calculates Relative Maximum Error (RME) and Relative Average Error (RAE)
    % Handles NaNs in input and division by zero if meas_data contains zeros.

    sim_data = sim_data(:); 
    meas_data = meas_data(:);

    if length(sim_data) ~= length(meas_data)
        error('Simulated and measured data must have the same length.');
    end

    % Filter out pairs where sim_data or meas_data is NaN, or where meas_data is zero
    valid_indices = ~isnan(sim_data) & ~isnan(meas_data) & (meas_data ~= 0);

    if sum(valid_indices) == 0
        warning('No valid data points (non-NaN, non-zero measured) for RME/RAE calculation. RME and RAE will be NaN.');
        RME = NaN;
        RAE = NaN;
        return;
    end

    sim_data_valid = sim_data(valid_indices);
    meas_data_valid = meas_data(valid_indices);

    pointwise_relative_errors = abs((sim_data_valid - meas_data_valid) ./ meas_data_valid) * 100;

    % Ensure no Infs/NaNs remain if any edge case slipped through (e.g. extremely small meas_data_valid)
    pointwise_relative_errors = pointwise_relative_errors(isfinite(pointwise_relative_errors));

    if isempty(pointwise_relative_errors)
        warning('No valid finite pointwise relative errors to calculate RME/RAE after isfinite check.');
        RME = NaN;
        RAE = NaN;
        return;
    end

    RME = max(pointwise_relative_errors);
    RAE = mean(pointwise_relative_errors);
end

function [D_val, SD_val, LoA_lower, LoA_upper, PCT1_val, D_max_abs_val, MNE_val, PCT2_val] = ...
    calculate_bland_altman_metrics(sim_data, meas_data)
    % Calculates Bland-Altman analysis metrics.
    % Handles NaNs in input data.

    sim_data = sim_data(:); 
    meas_data = meas_data(:); 

    if length(sim_data) ~= length(meas_data)
        error('Simulated and measured data must have the same length.');
    end

    % Remove pairs where either sim_data or meas_data is NaN
    valid_pair_indices = ~isnan(sim_data) & ~isnan(meas_data);
    sim_data_valid = sim_data(valid_pair_indices);
    meas_data_valid = meas_data(valid_pair_indices);

    if isempty(sim_data_valid) % Check if any valid pairs remain
        warning('No valid non-NaN data pairs for Bland-Altman calculation. All metrics will be NaN.');
        D_val=NaN; SD_val=NaN; LoA_lower=NaN; LoA_upper=NaN; PCT1_val=NaN; D_max_abs_val=NaN; MNE_val=NaN; PCT2_val=NaN;
        return;
    end

    differences = sim_data_valid - meas_data_valid; 

    D_val = mean(differences);        
    SD_val = std(differences);         

    % Check if D_val or SD_val became NaN (e.g., if only one valid_pair_index, std is NaN)
    if isnan(D_val) || isnan(SD_val) 
        LoA_lower = NaN;
        LoA_upper = NaN;
        PCT1_val = NaN; 
    else
        LoA_lower = D_val - 1.96 * SD_val;
        LoA_upper = D_val + 1.96 * SD_val;
        if isempty(differences) % Should be caught by isempty(sim_data_valid)
             PCT1_val = NaN;
        else
             PCT1_val = sum(differences < LoA_lower | differences > LoA_upper) / length(differences) * 100;
        end
    end

    if isempty(differences)
        D_max_abs_val = NaN;
    else
        D_max_abs_val = max(abs(differences));
    end

    MNE_val = mean((sim_data_valid + meas_data_valid) / 2); % MNE can be NaN if sim_data_valid is empty

    if ~isnan(MNE_val) && MNE_val ~= 0 && ~isnan(D_max_abs_val)
        PCT2_val = (D_max_abs_val / MNE_val) * 100; % PCT2 can be negative if MNE is negative
    else
        if ~isnan(MNE_val) && MNE_val == 0
             warning('MNE is zero. PCT2 cannot be calculated and is set to NaN.');
        end
        PCT2_val = NaN;
    end
end

%% --- Data Export to Excel for Indoor Temperatures ---

% Define the time points for extraction (every 48 minutes until 72 hours)
extraction_interval_minutes = 48; % minutes
extraction_interval_hours = extraction_interval_minutes / 60; % hours
t_extract_hours = 0:extraction_interval_hours:72;

% Ensure the extraction time points are within the simulated time range
t_extract_hours = t_extract_hours(t_extract_hours <= t_fine_hours(end));

% Interpolate data at the specified extraction time points
try
    T_indoor_sim_MS0_extracted = interp1(t_fine_hours, T_indoor_sim_MS0, t_extract_hours, 'linear', NaN);
    Measured_T_indoor_MS0_extracted = interp1(t_fine_hours, Measured_T_indoor_MS0, t_extract_hours, 'linear', NaN);

    T_indoor_sim_MS1_extracted = interp1(t_fine_hours, T_indoor_sim_MS1, t_extract_hours, 'linear', NaN);
    Measured_T_indoor_MS1_extracted = interp1(t_fine_hours, Measured_T_indoor_MS1, t_extract_hours, 'linear', NaN);

    T_indoor_sim_MS2_extracted = interp1(t_fine_hours, T_indoor_sim_MS2, t_extract_hours, 'linear', NaN);
    Measured_T_indoor_MS2_extracted = interp1(t_fine_hours, Measured_T_indoor_MS2, t_extract_hours, 'linear', NaN);

    disp('Data interpolated for extraction.');

    % Create tables for each scenario
    tbl_MS0 = table(t_extract_hours', T_indoor_sim_MS0_extracted', Measured_T_indoor_MS0_extracted', ...
                    'VariableNames', {'Time_Hours', 'Simulated_Indoor_Temp_C', 'Measured_Indoor_Temp_C'});

    tbl_MS1 = table(t_extract_hours', T_indoor_sim_MS1_extracted', Measured_T_indoor_MS1_extracted', ...
                    'VariableNames', {'Time_Hours', 'Simulated_Indoor_Temp_C', 'Measured_Indoor_Temp_C'});

    tbl_MS2 = table(t_extract_hours', T_indoor_sim_MS2_extracted', Measured_T_indoor_MS2_extracted', ...
                    'VariableNames', {'Time_Hours', 'Simulated_Indoor_Temp_C', 'Measured_Indoor_Temp_C'});

    % Define file names
    filename_MS0 = 'Indoor_Temperatures_MS0 summer.xlsx';
    filename_MS1 = 'Indoor_Temperatures_MS1 summer.xlsx';
    filename_MS2 = 'Indoor_Temperatures_MS2 summer.xlsx';

    % Write tables to Excel files
    writetable(tbl_MS0, filename_MS0);
    disp(['Data for MS0 saved to: ', filename_MS0]);

    writetable(tbl_MS1, filename_MS1);
    disp(['Data for MS1 saved to: ', filename_MS1]);

    writetable(tbl_MS2, filename_MS2);
    disp(['Data for MS2 saved to: ', filename_MS2]);

catch ME
    warning('Error during data extraction and Excel export: %s', ME.message);
end
