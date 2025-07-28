%% 1. LOAD & PREPROCESS DATA
load('Spring2k_1.mat', 't','x_r','v_r','x_h','v_h','f_r','f_int','a_r','a_h');

% Remove initial offsets if necessary (optional)
% x_r = x_r - x_r(1);
% x_h = x_h - x_h(1);
% f_r = f_r - f_r(1);

% Define sampling time and number of samples
Ts = t(2) - t(1);
N  = numel(t);

% Package signals into a structured variable
data.t     = t;
data.x_r   = x_r;
data.v_r   = v_r;
data.x_h   = x_h;
data.v_h   = v_h;
data.f_r   = f_r;
data.f_int = f_int;
data.a_r   = a_r;
data.a_h   = a_h;
data.N     = N;
data.Ts    = Ts;

% Define physical system constants
mr = 16.0;    % Robot mass (kg)
mh =  8.9;    % Human mass (kg)
k  = 20000;   % Stiffness (N/m)
c  =  175;    % Damping (N*s/m)

%% 2. DEFINE OPTIMIZATION VARIABLES & BOUNDS

nQ = 4;  nR = 6;                      % Sizes of Q and R diagonals
lb = [ -12*ones(1,nQ),  -7*ones(1,nR) ];  % Lower bounds (log10 scale)
ub = [  -2*ones(1,nQ),   2*ones(1,nR) ];  % Upper bounds (log10 scale)
nVars = nQ + nR;                   % Total number of optimization variables

% Define objective function handle
objFun = @(p) kfRmseObj(p, data, mr, mh, k, c);

%% 3. RUN PARTICLE SWARM OPTIMIZATION

opts = optimoptions('particleswarm', ...
    'Display','iter',...
    'SwarmSize', 60,...
    'MaxIterations', 100,...
    'UseParallel', true);

[p_opt, fval_opt] = particleswarm(objFun, nVars, lb, ub, opts);
fprintf('\nOptimal RMSE = %.5f\n', fval_opt);

%% 4. DECODE OPTIMAL Q AND R FROM SOLUTION

q_opt = 10.^(p_opt(1:nQ));         % Convert from log scale
r_opt = 10.^(p_opt(nQ+1:end));

Q_opt = diag(q_opt);
R_opt = diag(r_opt);

fprintf('Optimal Q diag: [%s]\n', sprintf('%.2e ', q_opt));
fprintf('Optimal R diag: [%s]\n', sprintf('%.2e ', r_opt));

%% 5. RE-RUN KF WITH OPTIMAL Q, R

[F_int_est, x_est, P_est] = runKF(Q_opt, R_opt, data, mr, mh, k, c);

%% 6. PLOT INTERACTION FORCE & ERROR METRICS

figure; hold on;
plot(t, data.f_int, 'k-', 'LineWidth',1);       % Ground truth
plot(t, F_int_est, 'r--','LineWidth',1);        % KF estimated force
xlabel('Time (s)'); ylabel('Force (N)');
legend('True','KF Estimate');
title('Interaction Force');
grid on;

% Compute errors
err   = F_int_est.' - data.f_int;
rmse  = sqrt(mean(err.^2));
mae   = mean(abs(err));
bias  = mean(err);
maxE  = max(abs(err));  
minE  = min(abs(err));

% Print error metrics
fprintf('\n=== Optimized KF Error ===\n');
fprintf('Bias: %8.4f\nMAE:  %8.4f\nRMSE: %8.4f\n', bias, mae, rmse);
fprintf('MaxE: %8.4f\nMinE: %8.4f\n', maxE, minE);

%% ===== LOCAL FUNCTIONS =====

function rmse = kfRmseObj(p, data, m1, m2, k, c)
    % Unpack parameters
    nQ = 4;
    q = 10.^(p(1:nQ));
    r = 10.^(p(nQ+1:end));
    Q = diag(q);
    R = diag(r);
    % Run KF with current Q and R
    F_est = runKF(Q, R, data, m1, m2, k, c);
    % Compute RMSE
    err  = F_est.' - data.f_int;
    rmse = sqrt(mean(err.^2));
end

function [F_int_est, x_est, P_est] = runKF(Q, R, data, m1, m2, k, c)
    % Kalman Filter Implementation
    Ts = data.Ts; N = data.N; n = 4;
    x_est     = zeros(n, N);
    P_est     = zeros(n,n,N);
    F_int_est = zeros(1, N);

    % Initialization
    x_est(:,1)   = [data.x_r(1); data.v_r(1); data.x_h(1); data.v_h(1)];
    P_est(:,:,1) = diag([0.1,0.1,0.1,0.1]);

    for i = 2:N
        % Inputs
        z      = [data.x_r(i); data.v_r(i); data.x_h(i); data.v_h(i); data.a_r(i); data.a_h(i)];
        x_prev = x_est(:,i-1);
        P_prev = P_est(:,:,i-1);
        F_app  = data.f_r(i-1);

        % Prediction Step
        x_pred = fProcess_RK4_KF(x_prev, F_app, Ts, m1, m2, k, c);
        A      = processJacobian_KF(x_prev, m1, m2, k, c);
        Fd     = eye(n) + Ts * A;
        P_pred = Fd*P_prev*Fd.' + Q;

        % Update Step
        H      = measurementJacobian_KF(x_pred, m1, m2, k, c);
        z_pred = hMeas_KF(x_pred, F_app, m1, m2, k, c);
        S      = H*P_pred*H.' + R;
        K      = P_pred * H' / S;
        nu     = z - z_pred;
        x_upd  = x_pred + K*nu;
        P_upd  = (eye(n)-K*H)*P_pred;

        % Store
        x_est(:,i)   = x_upd;
        P_est(:,:,i) = P_upd;
        F_int_est(i) = k*(x_upd(1)-x_upd(3)) + c*(x_upd(2)-x_upd(4));
    end
end

function x_next = fProcess_RK4_KF(x, F_app, Ts, m1, m2, k, c)
    % Runge-Kutta 4th order integration
    function dx = f_ct(xx)
        x1 = xx(1); v1 = xx(2);
        x2 = xx(3); v2 = xx(4);
        Fint = k*(x1-x2) + c*(v1-v2);
        a1   = (F_app - Fint)/m1;
        a2   = Fint / m2;
        dx   = [v1; a1; v2; a2];
    end
    k1 = f_ct(x);
    k2 = f_ct(x + 0.5*Ts*k1);
    k3 = f_ct(x + 0.5*Ts*k2);
    k4 = f_ct(x + Ts*k3);
    x_next = x + Ts/6*(k1 + 2*k2 + 2*k3 + k4);
end

function A = processJacobian_KF(x, m1, m2, k, c)
    % Linearization of the process model
    A = zeros(4,4);
    A(1,2)=1;
    A(2,1)=-k/m1; A(2,2)=-c/m1; A(2,3)= k/m1; A(2,4)= c/m1;
    A(3,4)=1;
    A(4,1)= k/m2; A(4,2)= c/m2; A(4,3)=-k/m2; A(4,4)=-c/m2;
end

function H = measurementJacobian_KF(x, m1, m2, k, c)
    % Linearization of the measurement model
    H = zeros(6,4);
    H(1,1)=1; H(2,2)=1; H(3,3)=1; H(4,4)=1;
    H(5,1)=-k/m1; H(5,2)=-c/m1; H(5,3)= k/m1; H(5,4)= c/m1;
    H(6,1)= k/m2; H(6,2)= c/m2; H(6,3)=-k/m2; H(6,4)=-c/m2;
end

function z = hMeas_KF(x, F_app, m1, m2, k, c)
    % Nonlinear measurement model
    x1 = x(1); v1 = x(2);
    x2 = x(3); v2 = x(4);
    Fint = k*(x1-x2) + c*(v1-v2);
    a1   = (F_app - Fint)/m1;
    a2   = Fint / m2;
    z    = [x1; v1; x2; v2; a1; a2];
end
