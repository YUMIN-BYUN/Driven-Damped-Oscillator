%% Free Damped Oscillator
% 1.Under,Over,Critical Damping(theoretical solution)
% 2.Under,Over,Critical Damping(numerical solution)(ode 45)
% 3.Comparison theoretical solution and numerical solution(ode 45)

% mx'' + bx' + kx = 0
% standard form
% x'' + 2γx' + ω0^2x = 0
% γ = b/(2m), ω0 = sqrt(k/m)

clear; clc; close all
color = ['b', 'r', 'g'];
gamma_list  = [0.2, 1.0, 2.0];
w0_val = 1.0 ;

%% Under,Over,Critical Damping(theoretical solution)
syms x(t) gamma w0
t_plot = linspace(0,20,2000);

ode_general = diff(x,t,2) + 2*gamma*diff(x,t) + w0^2*x == 0;
cond1 = x(0) == 1;
cond2 = subs(diff(x,t), t, 0) == 0;

fig1 = figure;
hold on

for i = 1:length(gamma_list)
    ode = subs(ode_general,{gamma,w0},{gamma_list(i),w0_val});
    sol(i) = dsolve(ode,cond1,cond2);
    fun = matlabFunction(sol(i));
    h_th(i) = plot(t_plot,fun(t_plot),color(i),'LineWidth',1.5);

    if gamma_list(i) < w0_val %Under Damping
        wd = sqrt(w0_val^2 - gamma_list(i)^2);
        A  = sqrt(1 + (gamma_list(i)/wd)^2);
        env = A * exp(-gamma_list(i)*t_plot);
        h_env = plot(t_plot,env,'--k');
        plot(t_plot,-env,'--k')
    end

end

sol_th = sol;

xlabel('Time')
ylabel('Displacement x(t)')
title('Free Damped Oscillator-theoretical sol')
legend([h_th(1), h_env, h_th(2), h_th(3)], ...
    'Under (\gamma(=0.2) < \omega_0)', ...
    'Envelope of UnderDamping', ...
    'Critical (\gamma(=1.0) = \omega_0)', ...
    'Over (\gamma(=2.0) > \omega_0)')
grid on

disp('Underdamping solution:')
disp(sol(1))
disp('Critical damping solution:')
disp(sol(2))
disp('Overdamping solution:')
disp(sol(3))

%% Under, Over, Critical Damping(numerical solution)
fig2 = figure; 
hold on

tspan = [0 20];
w0 = w0_val;

for j = 1:length(gamma_list)
    gamma = gamma_list(j);
    ode = @(t, x) [x(2); -2*gamma*x(2) - w0^2*x(1)];
    x0 = [1; 0];
    [t,x] = ode45(ode,tspan,x0);
    sol_nu{j} = x(:,1);
    T_nu{j} = t;
    h_nu(j) = plot(t,x(:,1),color(j),'LineWidth',1.5);

    if gamma_list(j) < w0
        [pks, locs] = findpeaks(abs(x(:,1)), t);
        log_pks = log(pks);
        coef = polyfit(locs, log_pks, 1);
        gamma_est = -coef(1);    
        A_est = exp(coef(2));

        env_est = A_est * exp(-gamma_est*t);
        h_env_nu = plot(t, env_est,'--k');
        plot(t,-env_est,'--k')
    end
end

xlabel('Time')
ylabel('Displacement x(t)')
title('Free Damped Oscillator-numerical sol')
legend([h_nu(1), h_env_nu, h_nu(2), h_nu(3)], ...
       'Under (\gamma(=0.2) < \omega_0)', ...
    'Envelope of UnderDamping', ...
    'Critical (\gamma(=1.0) = \omega_0)', ...
    'Over (\gamma(=2.0) > \omega_0)')
grid on
fprintf('True gamma = %.3f, Estimated gamma = %.3f\n', ...
         gamma_list(1), gamma_est);

%% Comparison theoretical solution and numerical solution
for k = 1:length(gamma_list)
    fun_th = matlabFunction(sol_th(k));
    x_nu = sol_nu{k};
    t_nu = T_nu{k};
    x_th = fun_th(t_nu);
    error = x_nu-x_th;

    L2_error = sqrt(trapz(t_nu,error.^2));
    max_error = max(abs(error));

    fprintf('gamma = %.2f | L2 error = %.3e | Max error = %.3e\n', ...
        gamma_list(k), L2_error, max_error);

end

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Free_damped_theoretical_sol.png')
saveas(fig2,'figures/Free_damped_numerical_sol.png')