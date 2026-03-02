%% Driven Damped Oscillator 
% 1.Fixed γ, varying ω(Steady-state)
% 2.Fixed ω, varying γ(Steady-state)


clear; clc; close all
w0_val = 1.0;
F_val = 1.0;
gamma_fixed = 0.1;
w_fixed = 1.0;
w_list = 0.05:0.05:1.5;
gamma_list = 0.05:0.05:1.5;

%% Fixed γ, varying ω(Steady-state)
gamma = gamma_fixed;
w0 = w0_val;
F = F_val;
tspan = [0 150];

for i = 1:length(w_list)

    %Numerical solution
    w = w_list(i);
    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w*t)];
    x0 = [0, 0];
    [t,x] = ode45(ode_nu,tspan,x0);
    idx = t > 100;
    A_nu_w(i) = max(abs(x(idx,1)));

    %Theoretical solution
    A_th_w(i) = F/sqrt((w0^2-w^2)^2+(2*gamma*w)^2);
    
end

%Plot
fig1 = figure;
subplot(2,1,1)
plot(w_list, A_nu_w, 'b')
hold on
plot(w_list, A_th_w, '--r')
title('Amplitude vs \omega')
xlabel('\omega')
ylabel('Amplitude')
legend('Numerical','Theoretical')
grid on

subplot(2,1,2)
plot(w_list, A_th_w - A_nu_w, 'k')
title('Difference (Theoretical - Numerical)')
xlabel('\omega')
ylabel('Error')
grid on

%Error calculating
error = A_th_w-A_nu_w;
max_error = max(abs(error));
fprintf('MAX error for varying omega = %.6f\n', max_error)

%% Fixed ω, varying γ(Steady-state)
w = w_fixed;
w0 = w0_val;
F = F_val;
tspan = [0 150];

for j = 1:length(gamma_list)

    %Numerical solution
    gamma = gamma_list(j);
    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w*t)];
    x0 = [0, 0];
    [t,x] = ode45(ode_nu,tspan,x0);
    idx = t > 100;
    A_nu_g(j) = max(abs(x(idx,1)));

    %Theoretical solution
    A_th_g(j) = F/sqrt((w0^2-w^2)^2+(2*gamma*w)^2);

end

%Plot
fig2 = figure;
subplot(2,1,1)
plot(gamma_list, A_nu_g, 'b')
hold on
plot(gamma_list, A_th_g, '--r')
title('Amplitude vs \gamma')
xlabel('\gamma')
ylabel('Amplitude')
legend('Numerical','Theoretical')
grid on

subplot(2,1,2)
plot(gamma_list, A_th_g - A_nu_g, 'k')
title('Difference (Theoretical - Numerical)')
xlabel('\gamma')
ylabel('Error')
grid on

%Error calculating
error = A_th_g - A_nu_g;
max_error = max(abs(error));
fprintf('MAX error for varying gamma = %.6f\n', max_error)

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Driven_damped_fixed_gamma_varying_omega.png')
saveas(fig2,'figures/Driven_damped_fixed_omega_varying_gamma.png')