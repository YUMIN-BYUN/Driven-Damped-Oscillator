%% Driven Damped Oscillator 
% 1.Theoretical solution
% 2.Numerical solution(using ODE solver 'ode45')
% 3.Calulating Error

% mx'' + bx' + kx = F0cos(wt) (Newton equation)
% standard form
% x'' + 2γx' + ω0^2x = Fcos(wt)
% where γ = b/(2m), ω0 = sqrt(k/m), F = F0/m; 

clear; clc; close all
w0_val = 1.0;
gamma_val = 0.1;
w_list = [0.7, 1.0, 1.3];
t_plot = linspace(0,80,5000);
color = ['b', 'r', 'g'];
F_val = 1.0;

fig1 = figure;
subplot(2,1,1)
hold on

%% Theoretical solution
syms x(t) gamma w w0 F 
ode_general = diff(x,t,2) + 2*gamma*diff(x,t,1) + w0^2*x - F*cos(w*t);
cond1 = x(0) == 1;
cond2 = subs(diff(x,t,1),t,0) == 0;

for i = 1:length(w_list)
    ode = subs(ode_general,{gamma,w0,w,F},{gamma_val,w0_val,w_list(i),F_val});
    sol_th(i) = dsolve(ode,cond1,cond2);
    fun_th = matlabFunction(sol_th(i));
    h_th(i) = plot(t_plot,fun_th(t_plot),color(i),'LineWidth',1.5);
end

title('Driven Damped Oscillator-theoretical solution')
xlabel('Time')
ylabel('Displacement x(t)')
legend([h_th(1), h_th(2), h_th(3)], ...
    '\omega = 0.7', ...
    '\omega = 1.0', ...
    '\omega = 1.3')
grid on

disp('Theoretical solution(w=0.7, w0=1.0, gamma=0.1, F=1.0)')
disp(sol_th(1))
disp('Theoretical solution(w=1.0, w0=1.0, gamma=0.1, F=1.0)')
disp(sol_th(2))
disp('Theoretical solution(w=1.3, w0=1.0, gamma=0.1, F=1.0)')
disp(sol_th(3))

%% Numerical solution(using ODE solver 'ode45')
gamma = gamma_val;
w0 = w0_val;
F = F_val;
tspan = [0 80];

subplot(2,1,2)
hold on


for j = 1:length(w_list)
    ode_nu = @(t,x) [x(2); -2*gamma*x(2) - w0^2*x(1) + F*cos(w_list(j)*t)];
    x0 = [1, 0];
    [t,x] = ode45(ode_nu,tspan,x0);
    sol_nu{j} = x(:,1);
    T_nu{j} = t;
    h_nu(j) = plot(t, x(:,1), color(j),'LineWidth', 1.5);
end

title('Driven Damped Oscillator-numerical solution')
xlabel('Time')
ylabel('Displacement x(t)')
legend([h_nu(1), h_nu(2), h_nu(3)], ...
    '\omega = 0.7', ...
    '\omega = 1.0', ...
    '\omega = 1.3')
grid on

%% Calculating Error
for k = 1:length(w_list)
    t_nu = T_nu{k};
    fun_th = matlabFunction(sol_th(k));
    x_th = fun_th(t_nu);
    x_nu = sol_nu{k};
    error = x_nu-x_th;

    L2_error = sqrt(trapz(t_nu,error.^2) / (t_nu(end)-t_nu(1)));
    max_error = max(abs(error));

    fprintf('omega = %.2f | L2 error = %.3e | Max error = %.3e\n', ...
        w_list(k), L2_error, max_error);

end

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Driven_damped_theoretical_numerical_sol.png')
