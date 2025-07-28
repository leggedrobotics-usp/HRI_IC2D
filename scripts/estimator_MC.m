%% 1. Initialization
clc
clear all
close all

% Load experimental dataset (Spring test with 20kN/m spring)
load('Spring20k_1.mat');  % or define your vectors manually if needed

%% 2. Define Settings and Parameters
W = 20;            % Sliding window size
N = 5000;          % Monte Carlo samples
time = t;          % Time vector

% Define parameter ranges and robot mass
state = struct();
state.W = W;
state.N = N;
state.k_range = [500 30000];  % Stiffness range (N/m)
state.b_range = [10 500];     % Damping range (Ns/m)
state.mr = 16.0;              % Robot mass (kg)

% Initialize signal history buffers
state.xr_hist = zeros(W,1);
state.xh_hist = zeros(W,1);
state.vr_hist = zeros(W,1);
state.vh_hist = zeros(W,1);
state.ar_hist = zeros(W,1);
state.fr_hist = zeros(W,1);

%% 3. Preallocation
W = state.W;     
N = state.N;     
N_steps = length(x_r);  % Total number of steps
mr = state.mr;
k_range = state.k_range;
b_range = state.b_range;
timestamps = zeros(N_steps, 1);  % Timing for each iteration

% Remove position offsets
x_r = x_r - x_r(1);
x_h = x_h - x_h(1);

% Logs for estimates and confidence
k_log       = zeros(N_steps, 1);
b_log       = zeros(N_steps, 1);
f_int_log   = zeros(N_steps, 1);
k_conf_log  = zeros(N_steps, 1);
b_conf_log  = zeros(N_steps, 1);

% Buffers for adaptive filtering
filtered_err     = zeros(N, 1);
filtered_k       = zeros(N, 1);
filtered_b       = zeros(N, 1);
filtered_f_int   = zeros(N, 1);

%% 4. Main Loop for Online Estimation
for t = 1:N_steps
    tic
    % 4.1 Get current values
    xr = x_r(t); vr = v_r(t);
    xh = x_h(t); vh = v_h(t);
    ar = a_r(t); fr = f_r(t);

    % 4.2 Update history
    state.xr_hist = [state.xr_hist(2:end); xr];
    state.xh_hist = [state.xh_hist(2:end); xh];
    state.vr_hist = [state.vr_hist(2:end); vr];
    state.vh_hist = [state.vh_hist(2:end); vh];
    state.ar_hist = [state.ar_hist(2:end); ar];
    state.fr_hist = [state.fr_hist(2:end); fr];

    % 4.3 Estimate interaction force from dynamics
    f_i_meas = state.fr_hist - mr * state.ar_hist;
    dx = state.xr_hist - state.xh_hist;
    dv = state.vr_hist - state.vh_hist;

    % 4.4 Adaptive filtering of previous estimates
    M = 0;
    if isfield(state, 'err_list_old')
        mean_err = mean(state.err_list_old);
        for i = 1:N
            if state.err_list_old(i) < 1 * mean_err
                M = M + 1;
                filtered_err(M)     = state.err_list_old(i);
                filtered_k(M)       = state.k_list_old(i);
                filtered_b(M)       = state.b_list_old(i);
                filtered_f_int(M)   = state.f_int_list_old(i);
            end
        end 
    end

    % 4.5 Monte Carlo sampling
    k_list = zeros(N,1);
    b_list = zeros(N,1);
    err_list = zeros(N,1);
    f_int_list = zeros(N,1);

    k_listmc = zeros(N,1);
    b_listmc = zeros(N,1);
    err_listmc = zeros(N,1);
    f_int_listmc = zeros(N,1);

    M_new = N - M;  % Number of new MC samples

    for i = 1:M_new
        k = rand() * diff(state.k_range) + state.k_range(1);
        b = rand() * diff(state.b_range) + state.b_range(1);
        f_i_est = k * dx + b * dv;
        err = mean((f_i_est - f_i_meas).^2);

        k_listmc(i) = k;
        b_listmc(i) = b;
        err_listmc(i) = err;
        f_int_listmc(i) = mean(f_i_est);
    end

    % 4.6 Merge old and new estimates
    err_list(1:M_new)       = err_listmc(1:M_new);
    k_list(1:M_new)         = k_listmc(1:M_new);
    b_list(1:M_new)         = b_listmc(1:M_new);
    f_int_list(1:M_new)     = f_int_listmc(1:M_new);

    err_list(M_new+1:end)   = filtered_err(1:M);
    k_list(M_new+1:end)     = filtered_k(1:M);
    b_list(M_new+1:end)     = filtered_b(1:M);
    f_int_list(M_new+1:end) = filtered_f_int(1:M);

    % 4.7 Final parameter estimates and confidence
    k_est = mean(k_list);
    b_est = mean(b_list);
    f_i_out = k_est * (xr - xh) + b_est * (vr - vh);

    k_conf = std(k_list);
    b_conf = std(b_list);

    % 4.8 Save results
    k_log(t) = k_est;
    b_log(t) = b_est;
    f_int_log(t) = f_i_out;
    k_conf_log(t) = k_conf;
    b_conf_log(t) = b_conf;

    % 4.9 Store for next iteration
    state.k_list_old = k_list;
    state.b_list_old = b_list;
    state.err_list_old = err_list;
    state.f_int_list_old = f_int_list;

    timestamps(t) = toc;
end

%% 5. Performance Metrics
erro   = f_int_log - (f_int(1:length(f_int_log),1)-f_int(1));
rmse   = sqrt(mean(erro.^2));
mae    = mean(abs(erro));
bias   = mean(erro);
maxE   = max(abs(erro)); 
minE   = min(abs(erro));
ts     = mean(timestamps);
stdts  = std(timestamps);
total_time = sum(timestamps);

[h, p, ci, stats] = ttest2(f_int, f_int_log);

% Print results
fprintf('\n=== MC Error ===\n');
fprintf('Bias: %8.4f\nMAE:  %8.4f\nRMSE: %8.4f\n', bias, mae, rmse);
fprintf('MaxE: %8.4f\nMinE: %8.4f\n', maxE, minE);
fprintf('Sample time: %.3e\nStd sample time: %.3e\n', ts, stdts);
fprintf('Two-sample t-test results:\n');
fprintf('Hypothesis test result (h): %d\n', h);
fprintf('p-value: %.4e\n', p);
fprintf('95%% Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('Degrees of freedom: %d\n', stats.df);
fprintf('Estimated std deviation: %.4f\n', stats.sd);

%% 6. Plot
figure; hold on;
darkBrown = [101, 67, 33]/255;  % Ground truth color
plot(time, f_int - f_int(1), 'Color', darkBrown, 'LineWidth', 1.5);
plot(time, f_int_log, '-r', 'LineWidth', 1.2);
xlabel('Time (s)', 'FontSize', 16);
ylabel('Interaction Force (N)', 'FontSize', 16);
xlim([0 41]); ylim([-70 55]);
legend('Groundtruth', 'MC', 'FontSize', 16, 'Location', 'northwest');
grid minor;
set(gcf, 'Units', 'inches', 'Position', [1 1 9 5]);
set(gca, 'FontSize', 16);
exportgraphics(gcf, 'MC_2k.pdf', 'BackgroundColor', 'none', 'Resolution',300);
