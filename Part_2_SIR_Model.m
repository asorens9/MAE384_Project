clear all; close all; clc;

%% Initial Conditions
S0 = 990;           % number of susceptible individuals at time t=0
I0 = 10;            % number of infected individuals at time t=0
R0 = 0;             % number of recovered individuals at time t=0
N = S0 + I0 + R0;   % total population

beta = 0.3;         % parameters for seasonal influenza
gamma = 0.1;        % parameters for seasonal influenza

%% Time parameters
t0 = 0;                % Initial time
T = 100;               % Final time
h = 2;                 % time step, [days]
t_even = t0:h:T;       % Even time steps

%% Initialize values from part 1
[t, S, I, R] = sir_rk4(S0, I0, R0, beta, gamma, N, h, T);
V_model_S = S';
V_model_I = I';
V_model_R = R';

%% Solve the SIR model (using Euler's method or any method)
S = zeros(size(t_even));
I = zeros(size(t_even));
R = zeros(size(t_even));

S(1) = S0;
I(1) = I0;
R(1) = R0;

for i = 1:length(t_even)-1
    % Euler's method to update S, I, R
    dS_dt = -beta * S(i) * I(i) / N;
    dI_dt = beta * S(i) * I(i) / N - gamma * I(i);
    dR_dt = gamma * I(i);
    
    S(i+1) = S(i) + dS_dt * h;
    I(i+1) = I(i) + dI_dt * h;
    R(i+1) = R(i) + dR_dt * h;
end

%% Interpolate the values at odd time steps using linear interpolation
t_odd = t0+1:h:T;  % Odd time steps (for interpolation)

%% Interpolate S(t) using linear interpolation
V_int_S_linear = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-1
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+1)
            x0 = t_even(i);
            x1 = t_even(i+1);
            V_int_S_linear(k) = S(i) + (t_odd(k) - x0) * (S(i+1) - S(i)) / (x1 - x0);
        end
    end
end

%% Interpolate I(t) using linear interpolation
V_int_I_linear = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-1
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+1)
            x0 = t_even(i);
            x1 = t_even(i+1);
            V_int_I_linear(k) = I(i) + (t_odd(k) - x0) * (I(i+1) - I(i)) / (x1 - x0);
        end
    end
end

%% Interpolate R(t) using linear interpolation
V_int_R_linear = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-1
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+1)
            x0 = t_even(i);
            x1 = t_even(i+1);
            V_int_R_linear(k) = R(i) + (t_odd(k) - x0) * (R(i+1) - R(i)) / (x1 - x0);
        end
    end
end

%% Calculate L2 error for S(t), I(t), and R(t) using linear interpolation
L2_error_S_linear = sum(sqrt(sum((V_int_S_linear - V_model_S).^2) / length(V_int_S_linear)));
L2_error_I_linear = sum(sqrt(sum((V_int_I_linear - V_model_I).^2) / length(V_int_I_linear)));
L2_error_R_linear = sum(sqrt(sum((V_int_R_linear - V_model_R).^2) / length(V_int_R_linear)));

%% Interpolate the values at odd time steps using quadratic interpolation
V_int_S_quad = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-2
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+2)
            x0 = t_even(i);
            x1 = t_even(i+1);
            x2 = t_even(i+2);
            
            f00 = S(i);
            f01 = (S(i+1) - S(i)) / (x1 - x0);
            f02 = ((S(i+2) - S(i+1)) / (x2 - x1) - f01) / (x2 - x0);
            
            V_int_S_quad(k) = f00 + (t_odd(k) - x0) * f01 + (t_odd(k) - x0) * (t_odd(k) - x1) * f02;
        end
    end
end

%% Interpolate I(t) using quadratic interpolation
V_int_I_quad = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-2
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+2)
            x0 = t_even(i);
            x1 = t_even(i+1);
            x2 = t_even(i+2);
            
            f00 = I(i);
            f01 = (I(i+1) - I(i)) / (x1 - x0);
            f02 = ((I(i+2) - I(i+1)) / (x2 - x1) - f01) / (x2 - x0);
            
            V_int_I_quad(k) = f00 + (t_odd(k) - x0) * f01 + (t_odd(k) - x0) * (t_odd(k) - x1) * f02;
        end
    end
end

%% Interpolate R(t) using quadratic interpolation
V_int_R_quad = zeros(size(t_odd));
for k = 1:length(t_odd)
    for i = 1:length(t_even)-2
        if t_odd(k) >= t_even(i) && t_odd(k) <= t_even(i+2)
            x0 = t_even(i);
            x1 = t_even(i+1);
            x2 = t_even(i+2);
            
            f00 = R(i);
            f01 = (R(i+1) - R(i)) / (x1 - x0);
            f02 = ((R(i+2) - R(i+1)) / (x2 - x1) - f01) / (x2 - x0);
            
            V_int_R_quad(k) = f00 + (t_odd(k) - x0) * f01 + (t_odd(k) - x0) * (t_odd(k) - x1) * f02;
        end
    end
end

%% Calculate L2 error for S(t), I(t), and R(t) using quadratic interpolation
L2_error_S_quad = sum(sqrt(sum((V_int_S_quad - V_model_S).^2) / length(V_int_S_quad)));
L2_error_I_quad = sum(sqrt(sum((V_int_I_quad - V_model_I).^2) / length(V_int_I_quad)));
L2_error_R_quad = sum(sqrt(sum((V_int_R_quad - V_model_R).^2) / length(V_int_R_quad)));

%% Display L2 error matrix
Errors = [
    L2_error_S_linear, L2_error_I_linear, L2_error_R_linear;    % Linear interpolation errors
    L2_error_S_quad, L2_error_I_quad, L2_error_R_quad           % Quadratic interpolation errors
];

disp('L2 error matrix (2x3):');
disp(Errors);
disp(['The linear Newton interpolation method results in smaller error ' ...
    'values than the quadratic Newton interpolation method.'])

%% Part 1 Code - Runge-Kutta Method
function [t, S, I, R] = sir_rk4(S0, I0, R0, beta, gamma, N, h, T)
    % Time vector
    t = 0:h:T;
    n_steps = length(t);
    
    % Initialize variables
    S = zeros(1, n_steps);
    I = zeros(1, n_steps);
    R = zeros(1, n_steps);
    
    % Set initial conditions
    S(1) = S0;
    I(1) = I0;
    R(1) = R0;
    
    % Solve using Runge-Kutta
    for i = 1:n_steps-1
        k1_S = -beta * S(i) * I(i) / N;
        k1_I = beta * S(i) * I(i) / N - gamma * I(i);
        k1_R = gamma * I(i);

        k2_S = -beta * (S(i) + 0.5*h*k1_S) * (I(i) + 0.5*h*k1_I) / N;
        k2_I = beta * (S(i) + 0.5*h*k1_S) * (I(i) + 0.5*h*k1_I) / N - gamma * (I(i) + 0.5*h*k1_I);
        k2_R = gamma * (I(i) + 0.5*h*k1_I);

        k3_S = -beta * (S(i) + 0.5*h*k2_S) * (I(i) + 0.5*h*k2_I) / N;
        k3_I = beta * (S(i) + 0.5*h*k2_S) * (I(i) + 0.5*h*k2_I) / N - gamma * (I(i) + 0.5*h*k2_I);
        k3_R = gamma * (I(i) + 0.5*h*k2_I);

        k4_S = -beta * (S(i) + h*k3_S) * (I(i) + h*k3_I) / N;
        k4_I = beta * (S(i) + h*k3_S) * (I(i) + h*k3_I) / N - gamma * (I(i) + h*k3_I);
        k4_R = gamma * (I(i) + h*k3_I);

        % Update values for the next step
        S(i+1) = S(i) + (h/6) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
        I(i+1) = I(i) + (h/6) * (k1_I + 2*k2_I + 2*k3_I + k4_I);
        R(i+1) = R(i) + (h/6) * (k1_R + 2*k2_R + 2*k3_R + k4_R);
    end
end
