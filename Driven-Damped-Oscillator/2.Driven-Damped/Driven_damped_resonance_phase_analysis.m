%% Driven Damped Oscillator 
% 1.Resonance frequency analysis
% 2.Phase analysis

clear; clc; close all

%% Resonance frequency analysis
w0 = 1.0;
F = 1.0;
gamma = 0.2;
w_list = 0.05:0.001:2.0;
tspan = [0 150];
A_nu = zeros(size(w_list));
A_th = zeros(size(w_list));

for i = 1:length(w_list)

    %Numerical solution
    w = w_list(i);
    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w*t)];
    x0 = [0, 0];
    [t,x] = ode45(ode_nu,tspan,x0);
    idx = t > 100;
    A_nu(i) = max(abs(x(idx,1)));
    

    %Theoretical solution
    A_th(i) = F/sqrt((w0^2-w^2)^2+(2*gamma*w)^2);

end

[resonance_apm_nu, num_nu] = max(A_nu);
resonance_freq_nu = w_list(num_nu);
[resonance_apm_th, num_th] = max(A_th);
resonance_freq_th = w_list(num_th);

fig1 = figure;

subplot(2,1,1)
hold on
plot(w_list,A_nu,'b')
plot(resonance_freq_nu,resonance_apm_nu,'o')
title('Resonance analysis-numerically')
xlabel('\omega')
ylabel('Amplitude')
grid on

subplot(2,1,2)
hold on
plot(w_list,A_th,'r')
plot(resonance_freq_th,resonance_apm_th,'o')
title('Resonance analysis-theoretically')
xlabel('\omega')
ylabel('Amplitude')
grid on

fprintf('Resonace frequency(using numerical method) = %.4f\n', resonance_freq_nu)
fprintf('Resonace frequency(using theoretical method) = %.4f\n', resonance_freq_th)

%% Phase analysis
w0 = 1.0;
F = 1.0;
gamma = 0.2;
w_list = 0.05:0.001:2.0;
tspan = [0 150];

for j = 1:length(w_list)

    %Numerical solution
    w = w_list(j);
    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w*t)];
    x0 = [0, 0];
    [t,x] = ode45(ode_nu,tspan,x0);
    idx = t > 100;

    %fitting
    t_f = t(idx);
    x_f = x(idx,1);

    M = [cos(w*t_f), sin(w*t_f)];
    coeff = M \ x_f;
    
    C = coeff(1);
    D = coeff(2);
    
    phi_fit(j) = atan2(D, C);

    %Theoretical phase calculation
    phi_th(j) = atan2(2*gamma*w, w0^2 - w^2);
end

%plot

fig2 = figure;
subplot(2,1,1)
plot(w_list,phi_fit,'b')
title('Phase vs \omega - numerically')
xlabel('\omega')
ylabel('Phase(\phi)')
grid on

subplot(2,1,2)
plot(w_list,phi_th,'r')
title('Phase vs \omega - theoretically')
xlabel('\omega')
ylabel('Phase(\phi)')
grid on

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Driven_damped_resonace_frequency.png')
saveas(fig2,'figures/Driven_damped_phase_anlaysis.png')