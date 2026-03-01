%% Driven Damped Oscillator 
% 1.Seperating transient & steady-state using FFT(Numerical solution)
% 2.Theoretical solution(steady-state)
% 3.Calulating Error

clear; clc; close all

%% Seperating transient & steady-state using FFT(Numerical solution)
color = ['b', 'r', 'g'];
w_list = [0.7, 1.0, 1.3];
gamma = 0.1;
w0 = 1.0;
F = 1.0;
tspan = [0 200];

fig1 = figure;
hold on

for i = 1:length(w_list)

    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w_list(i)*t)];
    x0 = [1, 0];

    [t,x] = ode45(ode_nu,tspan,x0);

    dt_t = 0.01;
    t_uni = t(1):dt_t:t(end);
    x_uni = interp1(t,x(:,1),t_uni);


    % Seperating Steady-state
    t_ss = t_uni(t_uni > 100);
    x_ss = x_uni(t_uni > 100);
    x_ss = x_ss - mean(x_ss);
    N = length(x_ss);

    dt = t_ss(2) - t_ss(1);
    fs = 1/dt;
    
    X_fft = fft(x_ss);
    
    freq = (0:N-1)*(fs/N)*2*pi;
    
    A = abs(X_fft)/(N);
    A = A(1:N/2);
    A(2:end) = 2*A(2:end);
    freq = freq(1:N/2);
    
    [amp_peak, idx] = max(A);
    f_peak(i) = freq(idx);
    A_nu(i) = amp_peak;
    h(i) = plot(freq, A, color(i));
    plot(f_peak(i),A_nu(i),'o')

    fprintf('Driving w = %.2f, FFT peak angular freq = %.4f, FFT peak Amplitude = %.4f\n', ...
    w_list(i), f_peak(i), A_nu(i))

end

%% Steady-state amplitude(theoretical solution)
disp('---------------------------------------------------------------------')
for j = 1:length(w_list)
    A_th(j) = F/sqrt((w0^2-w_list(j)^2)^2 + 4*gamma^2*w_list(j)^2);
    fprintf('Driving w = %.2f, Theoretical Amplitude = %.4f\n', ...
        w_list(j), A_th(j))
    yline(A_th(j),'--','Color',color(j),'Label','A_{theory}')
end

%% Calulating Error
disp('---------------------------------------------------------------------')
for k = 1:length(w_list)
    rel_error = abs(A_nu(k) - A_th(k)) / A_th(k);
    rel_error_percent = rel_error * 100;
    fprintf('Driving w =%.2f, Relative error(%%) =  %.4f\n',w_list(k),rel_error_percent)
end

title('Fourier transform-numerical solution')
xlabel('Angular Frequency(rad/s)')
ylabel('Amplitude')
axis([0 2 0 5.5])
legend([h(1), h(2), h(3)], ...
    '\omega = 0.7', ...
    '\omega = 1.0', ...
    '\omega = 1.3','Location','northwest')
grid on

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Driven_damped_FFM_amplitude_sol.png')