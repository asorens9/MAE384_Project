%% Part 4
clear all; clc;

% Givens
B0 = 0.3;
A = 5;
w = 2 * pi * (365/365);
B = @(t) (B0*(1+A*sin(w*t)));

% Inital run of SIR using function from Part I
S0 = 990;
I0 = 10;
R0 = 0;
y = 0.1;
h = 0.1;
T = 150;
N = S0 + I0 + R0;

[t, S, I, R] = sir_rk4m(S0, I0, R0, y, N, h, T, B);

% Plot Signals
figure;
plot(t, S, 'b-', t, I, 'r-', t, R, 'g-')
legend('Susceptible','Infected','Recovered')
title('SIR Model (Period = 1 day)')
ylabel('Population')
xlabel('Days')

Sfft = fft(S);
Ifft = fft(I);
Rfft = fft(R);

f = (1/T)*(0:(N/2));
figure;
plot(f,abs(Ifft(1:length(f))),'k-')
title('FFT of I(t) (Period = 1 day)')
xlabel('frequency (Hz)')
ylabel('Amplitude')

%% Run it back for step 6

% New Givens
w = 2 * pi * (100/365);
B = @(t) (B0*(1+A*sin(w*t)));

[t, S, I, R] = sir_rk4m(S0, I0, R0, y, N, h, T, B);

% Plot Signals
figure;
plot(t, S, 'b-', t, I, 'r-', t, R, 'g-')
legend('Susceptible','Infected','Recovered')
title('SIR Model (Period ~3 days)')
ylabel('Population')
xlabel('Days')

Sfft = fft(S);
Ifft = fft(I);
Rfft = fft(R);

f = (1/T)*(0:(N/2));
figure;
plot(f,abs(Ifft(1:length(f))),'k-')
title('FFT of I(t) (Period ~3 days)')
xlabel('frequency (Hz)')
ylabel('Amplitude')

%% Function from Part I with slight modifications
function [t, S, I, R] = sir_rk4m(S0, I0, R0, y, N, h, T, B)
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
    
    % RK4 Method
    for i = 1:n_steps-1

        beta = B(i * 0.1);

        % Define functions for derivatives
        dS = @(S, I) -beta * S * I / N;
        dI = @(S, I) (beta * S * I / N) - y * I;
        dR = @(I) y * I;
        
        % k1
        k1S = h * dS(S(i), I(i));
        k1I = h * dI(S(i), I(i));
        k1R = h * dR(I(i));
        
        % k2
        k2S = h * dS(S(i) + k1S/2, I(i) + k1I/2);
        k2I = h * dI(S(i) + k1S/2, I(i) + k1I/2);
        k2R = h * dR(I(i) + k1I/2);
        
        % k3
        k3S = h * dS(S(i) + k2S/2, I(i) + k2I/2);
        k3I = h * dI(S(i) + k2S/2, I(i) + k2I/2);
        k3R = h * dR(I(i) + k2I/2);
        
        % k4
        k4S = h * dS(S(i) + k3S, I(i) + k3I);
        k4I = h * dI(S(i) + k3S, I(i) + k3I);
        k4R = h * dR(I(i) + k3I);
        
        % Update
        S(i+1) = S(i) + (k1S + 2*k2S + 2*k3S + k4S) / 6;
        I(i+1) = I(i) + (k1I + 2*k2I + 2*k3I + k4I) / 6;
        R(i+1) = R(i) + (k1R + 2*k2R + 2*k3R + k4R) / 6;
    end
end

