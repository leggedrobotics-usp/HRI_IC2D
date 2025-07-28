%% Adaptive EKF with Tighter Tuning and Real-Time Low-Pass Smoothing
clear all
close all
clc

%% 1. Load Data & Remove Offsets
% Load experimental measurements
load('Spring20k_1.mat', 't','x_r','v_r','x_h','v_h','f_r','f_int','a_r','a_h');

% Normalize positions and actuator force to start at zero
x_r = x_r - x_r(1);
x_h = x_h - x_h(1);
f_r = f_r - f_r(1);

Ts  = t(2) - t(1);  % Sampling time

%% 2. System & EKF Settings
% Masses of robot and human
mr = 16.0;    
mh = 8.9;

% State dimension: [x_r; v_r; x_h; v_h; k; c]
n = 6;  

% Measurement dimension: [x_r; v_r; x_h; v_h; a_r; a_h]
m = 6;  

% Initial process (Q) and measurement (R) noise covariances
Q = diag([1e-6, 1e-6, 1e-6, 1e-6, 1e-3, 1e-5]);  
R = diag([1e-7, 1e-7, 1e-7, 1e-7, 0.5, 0.5]);

% Adaptive copies with saturation bounds
Q_adapt = Q;
R_adapt = R;

Q_min = diag([1e-8,1e-8,1e-8,1e-8,1e-8,1e-8]);
Q_max = diag([1e2,1e2,1e2,1e2,1e4,1e3]);  

R_min = diag([1e-7,1e-7,1e-7,1e-7,1e-2,1e-2]);
R_max = diag([1e1, 1e1, 1e1, 1e1, 1e3, 1e3]);

% Buffer settings for empirical adaptation
M = 500;      % Buffer length
L = 100;      % Update interval

% Innovation and process residual buffers
innov    = zeros(m, M);
proc_res = zeros(n, M);
buf_idx  = 0;

%% 3. Preallocate Memory
N             = numel(t);
x_est         = zeros(n, N);
P_est         = zeros(n,n,N);
F_int_est     = zeros(1, N);
F_int_smooth  = zeros(1, N);
Q_hist        = zeros(n, N);   % Track Q over time
R_hist        = zeros(m, N);   % Track R over time

% Low-pass filter settings for smoothing
tau   = 0.1;
alpha = exp(-Ts/tau);               
F_int_smooth(1) = 0;

% Initial state: position/velocity from data, parameters with loose priors
x_est(:,1)   = [x_r(1); v_r(1); x_h(1); v_h(1); 2000; 175];
P_est(:,:,1) = diag([0.1,0.1,0.1,0.1,1e4,1e4]);

%% 4. Main EKF Loop
for i = 2:N
    tic

    % 4.1 Measurement
    z      = [x_r(i); v_r(i); x_h(i); v_h(i); a_r(i); a_h(i)];
    x_prev = x_est(:,i-1);
    P_prev = P_est(:,:,i-1);
    F_app  = f_r(i-1);

    % 4.2 Prediction (Runge-Kutta + Jacobian)
    x_pred = fProcess_RK4(x_prev, F_app, Ts, mr, mh);
    A_cont = processJacobian(x_prev, mr, mh);
    Fd     = eye(n) + Ts * A_cont;
    P_pred = Fd * P_prev * Fd.' + Q_adapt;

    % 4.3 Update Step
    H      = measurementJacobian(x_pred, mr, mh);
    z_pred = hMeas(x_pred, F_app, mr, mh);
    S      = H * P_pred * H.' + R_adapt;
    K      = P_pred * H.' / S;
    nu     = z - z_pred;

    x_upd  = x_pred + K * nu;
    P_upd  = (eye(n) - K*H) * P_pred;

    % 4.4 Adaptive R Update
    buf_idx = mod(buf_idx, M) + 1;
    innov(:,buf_idx) = nu;

    if i > M && mod(i, L) == 0
        S_emp   = (innov * innov.') / (M - 1);
        scaling = trace(S_emp) / trace(H * P_pred * H.' + R_adapt);
        R_adapt = min(max(R_adapt * scaling, R_min), R_max);
    end

    % 4.5 Adaptive Q Update
    proc_res(:,buf_idx) = x_upd - x_pred;

    if i > M && mod(i, L) == 0
        T_emp   = (proc_res * proc_res.') / (M - 1);
        scalingQ = trace(T_emp) / trace(Q_adapt);
        Q_adapt  = min(max(Q_adapt * scalingQ, Q_min), Q_max);
    end

    % 4.6 Store Estimates
    Q_hist(:,i) = diag(Q_adapt);
    R_hist(:,i) = diag(R_adapt);
    x_est(:,i)   = x_upd;
    P_est(:,:,i) = P_upd;

    % Compute interaction force using k and c estimates
    F_int_est(i) = x_upd(5)*(x_upd(1)-x_upd(3)) + ...
                   x_upd(6)*(x_upd(2)-x_upd(4));

    % Real-time low-pass filter
    F_int_smooth(i) = alpha * F_int_smooth(i-1) + ...
                      (1 - alpha) * F_int_est(i);

    timestamps(i) = toc;
end

%% 5. Plot Results: Interaction Force
figure; hold on;
% Define dark brown color RGB (approximate)
darkBrown = [101, 67, 33]/255;
plot(t, f_int - f_int(1), 'Color', darkBrown, 'LineWidth', 1.5);
plot(t, F_int_smooth, '-m', 'LineWidth', 1.5);
% Labels
xlabel('Time (s)', 'FontSize', 16);
xlim([0 41])
ylim([-50 40])
ylabel('Interaction Force (N)', 'FontSize', 16);
% Legend with font size
legend('Groundtruth', 'Smoothed EKF', 'Location', 'northwest', 'FontSize', 16);
grid minor;
% Set figure size (inches)
set(gcf, 'Units', 'inches', 'Position', [1 1 9 5]);
% Also update axes font size for tick labels
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'EKF_2k.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none', 'Resolution',300);

%% 6. Error Metrics (Smoothed)
err   = F_int_smooth.' - f_int;
rmse  = sqrt(mean(err.^2));
mae   = mean(abs(err));
bias  = mean(err);
maxE  = max(abs(err)); 
minE  = min(abs(err));
ts    = mean(timestamps);
stdts = std(timestamps);
total_time = sum(timestamps);

% Statistical test
[h, p, ci, stats] = ttest2(f_int, F_int_smooth.');

% Display results
fprintf('\n=== Smoothed EKF Error ===\n');
fprintf('Bias: %8.4f\nMAE:  %8.4f\nRMSE: %8.4f\n', bias, mae, rmse);
fprintf('MaxE: %8.4f\nMinE: %8.4f\n', maxE, minE);
fprintf('Sample time: %.3e\nStd sample time: %.3e\n', ts, stdts);
fprintf('Two-sample t-test results:\n');
fprintf('Hypothesis test result (h): %d\n', h);
fprintf('p-value: %.4e\n', p);
fprintf('95%% Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));
fprintf('t-statistic: %.4f\nDegrees of freedom: %d\nStd: %.4f\n', ...
        stats.tstat, stats.df, stats.sd);

%% ==== Local Functions ====

function x_next = fProcess_RK4(x, F_app, Ts, m1, m2)
    % RK4 integration of nonlinear dynamics
    function dx = f_ct(xx)
        x1 = xx(1); v1 = xx(2);
        x2 = xx(3); v2 = xx(4);
        k_ = xx(5); c_ = xx(6);
        F_int = k_*(x1 - x2) + c_*(v1 - v2);
        a1 = (F_app - F_int)/m1;
        a2 =  F_int / m2;
        dx = [v1; a1; v2; a2; 0; 0];
    end
    k1 = f_ct(x);
    k2 = f_ct(x + 0.5*Ts*k1);
    k3 = f_ct(x + 0.5*Ts*k2);
    k4 = f_ct(x + Ts*k3);
    x_next = x + Ts/6 * (k1 + 2*k2 + 2*k3 + k4);
end

function A = processJacobian(x, m1, m2)
    % Jacobian of dynamics w.r.t. state
    x1 = x(1); v1 = x(2); x2 = x(3); v2 = x(4);
    k_ = x(5); c_ = x(6);
    A = zeros(6,6);
    A(1,2) = 1;
    A(2,1) = -k_/m1;  A(2,2) = -c_/m1;
    A(2,3) =  k_/m1;  A(2,4) =  c_/m1;
    A(2,5) = -(x1 - x2)/m1; A(2,6) = -(v1 - v2)/m1;
    A(3,4) = 1;
    A(4,1) =  k_/m2;  A(4,2) =  c_/m2;
    A(4,3) = -k_/m2;  A(4,4) = -c_/m2;
    A(4,5) =  (x1 - x2)/m2; A(4,6) =  (v1 - v2)/m2;
end

function H = measurementJacobian(x, m1, m2)
    % Jacobian of measurement model w.r.t. state
    x1 = x(1); v1 = x(2); x2 = x(3); v2 = x(4);
    k_ = x(5); c_ = x(6);
    H = zeros(6,6);
    H(1,1)=1; H(2,2)=1; H(3,3)=1; H(4,4)=1;
    H(5,1) = -k_/m1;  H(5,2) = -c_/m1;  H(5,3) = k_/m1;  H(5,4) = c_/m1;
    H(5,5) = -(x1-x2)/m1; H(5,6) = -(v1-v2)/m1;
    H(6,1) =  k_/m2;  H(6,2) =  c_/m2;  H(6,3) = -k_/m2; H(6,4) = -c_/m2;
    H(6,5) =  (x1-x2)/m2; H(6,6) =  (v1-v2)/m2;
end

function z = hMeas(x, F_app, m1, m2)
    % Nonlinear measurement model
    x1 = x(1); v1 = x(2); x2 = x(3); v2 = x(4);
    k_ = x(5); c_ = x(6);
    F_int = k_*(x1 - x2) + c_*(v1 - v2);
    a1 = (F_app - F_int)/m1;
    a2 =  F_int / m2;
    z  = [x1; v1; x2; v2; a1; a2];
end
