%% 1) LOAD & PREPROCESS DATA
clear all
close all
clc

% Load time series data (positions, velocities, accelerations, forces)
load('Spring20k_1.mat','t','x_r','v_r','x_h','v_h','f_r','f_int','a_r','a_h');

% Normalize position and force offsets to start at zero
x_r = x_r - x_r(1);
x_h = x_h - x_h(1);
f_r = f_r - f_r(1);

% Define physical system constants
k  = 20000;   % Spring stiffness (N/m)
c  = 175;     % Damping coefficient (Ns/m)
mr = 16.0;    % Mass of robot (kg)
mh = 8.9;     % Mass of human (kg)

% Sampling and model size
Ts = t(2)-t(1);      % Time step
N  = numel(t);       % Number of samples
n  = 4;              % State dimension: [x_r; v_r; x_h; v_h]

%% 2) DEFINE SYSTEM MODEL (Discretized)

% Continuous-time dynamics
A_ct = [  0,       1,      0,       0;
       -k/mr, -c/mr,   k/mr,   c/mr;
          0,       0,      0,       1;
        k/mh,  c/mh,  -k/mh,  -c/mh ];

B_ct = [0; 1/mr; 0; 0];

% First-order Euler discretization
F = eye(n) + Ts*A_ct;
B = Ts * B_ct;

% Measurement model H: includes positions, velocities, accelerations
H = [ 1, 0, 0, 0;
      0, 1, 0, 0;
      0, 0, 1, 0;
      0, 0, 0, 1;
   -k/mr, -c/mr,  k/mr,  c/mr;
    k/mh,  c/mh, -k/mh, -c/mh ];

%% 3) INITIALIZE FILTER PARAMETERS

% Process noise covariance (Q) and measurement noise (R)
Q = diag([2e-7, 1e-12, 1e-12, 3e-5]);
R = diag([8e-5, 5e-6, 1e1, 1e-7, 1e2, 1e2]);

% Initial state estimate and covariance
x_est = zeros(n,N);
P_est = zeros(n,n,N);
x_est(:,1)   = [ x_r(1); v_r(1); x_h(1); v_h(1) ];
P_est(:,:,1) = diag([0.1,0.1,0.1,0.1]);

% Allocate for estimated and smoothed interaction force
F_int_est    = zeros(1,N);
F_int_smooth = zeros(1,N);

% Real-time low-pass filter constant
tau = 0.05;
alpha = exp(-Ts/tau);
F_int_smooth(1) = 0;

%% 4) KALMAN FILTER MAIN LOOP

for i = 2:N
    tic
    % Get measurements: position, velocity, acceleration
    z      = [ x_r(i); v_r(i); x_h(i); v_h(i); a_r(i); a_h(i) ];
    F_app  = f_r(i-1);        % Input: actuator force
    x_prev = x_est(:,i-1);
    P_prev = P_est(:,:,i-1);

    % --- Prediction Step ---
    x_pred = F*x_prev + B*F_app;
    P_pred = F*P_prev*F' + Q;

    % --- Update Step ---
    z_pred = H*x_pred;
    S      = H*P_pred*H' + R;
    K      = P_pred*H' / S;
    nu     = z - z_pred;

    x_upd = x_pred + K*nu;
    P_upd = (eye(n) - K*H)*P_pred;

    % Store state and covariance
    x_est(:,i) = x_upd;
    P_est(:,:,i) = P_upd;

    % Estimate interaction force from state
    F_int_est(i) = k*(x_upd(1)-x_upd(3)) + c*(x_upd(2)-x_upd(4));

    % Apply real-time low-pass filter
    F_int_smooth(i) = alpha*F_int_smooth(i-1) + (1-alpha)*F_int_est(i);
    timestamps(i) = toc;
end

%% 5) ERROR METRICS AND PLOTTING

err   = F_int_smooth.' - f_int;
rmse  = sqrt(mean(err.^2));
mae   = mean(abs(err));
bias  = mean(err);
maxE  = max(abs(err)); 
minE  = min(abs(err));
ts = mean(timestamps);
stdts = std(timestamps);
total_time = sum(timestamps);

[h, p, ci, stats] = ttest2(f_int, F_int_smooth.');

% Print error results
fprintf('\n=== Linear KF + LPF Error ===\n');
fprintf('Bias: %8.4f\nMAE:  %8.4f\nRMSE: %8.4f\n', bias, mae, rmse);
fprintf('MaxE: %8.4f\nMinE: %8.4f\n', maxE, minE);
fprintf('Sample time: %.3e\nStd sample time: %.3e\n', ts, stdts);
fprintf('Two-sample t-test results:\n');
fprintf('Hypothesis test result (h): %d\n', h);
fprintf('p-value: %.4e\n', p);
fprintf('95%% Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));
fprintf('t-statistic: %.4f, DoF: %d, Std: %.4f\n', stats.tstat, stats.df, stats.sd);

%% 6) PLOT ESTIMATED VS TRUE INTERACTION FORCE

figure; hold on;

% Plot ground-truth and estimated interaction force
plot(t, f_int - f_int(1), 'Color', [101, 67, 33]/255, 'LineWidth', 1.5); % dark brown
plot(t, F_int_smooth, 'Color', [0.3, 0.6, 1.0], 'LineWidth', 1.5);       % light blue

xlabel('Time (s)', 'FontSize', 16);
ylabel('Interaction Force (N)', 'FontSize', 16);
xlim([0 41])
ylim([-50 40])
legend('Groundtruth', 'KF', 'Location', 'northwest', 'FontSize', 16);
grid minor;
set(gca, 'FontSize', 16);
set(gcf, 'Units', 'inches', 'Position', [1 1 9 5]);

% Export figure
exportgraphics(gcf, 'KF_2k.pdf', 'ContentType', 'vector', ...
    'BackgroundColor', 'none', 'Resolution', 300);
