%% Free Damped Oscillator
% 1.Non-dimensionalization(theoretical solution)
% 2.Comparing system same zeta, but different ω0 & γ (Underdamping case)
% 3.Envelope analysis(Underdamping case)

% def τ = ω0*t (tau is dimensionless)
% ω0^2*x(τ)'' + 2*γ*ω0*x' + ω0^2*x = 0
% x'' + 2*(γ/ω0)*x' + x = 0
% def ζ = γ/ω0 (zeta is dimensionless)
% x'' + 2*ζ*x' + x = 0

clear; clc; close all
color = ['b', 'r', 'g'];
zeta_list = [0.2, 1.0, 2.0]; %Under, Critical, Overdamping
tau_plot = linspace(0,20,2000);

%% Solving theoretically(Non-dimensionalization)
syms x(tau) zeta
ode_general = diff(x,tau,2) + 2*zeta*diff(x,tau) + x == 0;
cond1 = x(0) == 1;
cond2 = subs(diff(x,tau), tau, 0) == 0;

fig1 = figure;
hold on

for i = 1:length(zeta_list)
    ode = subs(ode_general,{zeta},{zeta_list(i)});
    sol(i) = dsolve(ode,cond1,cond2);
    fun = matlabFunction(sol(i));
    h(i) = plot(tau_plot,fun(tau_plot),color(i),'LineWidth',1.5);

    if zeta_list(i) < 1.0
        zeta_d = sqrt(1-(zeta_list(i))^2);
        A = sqrt(1+(zeta_list(i)/zeta_d)^2);
        env = A*exp(-zeta_list(i)*tau_plot);
        h_env = plot(tau_plot,env,'--k');
        plot(tau_plot,-env,'--k');
    end

end
grid on
title('Free Damped Oscillator Solutions');
xlabel('Dimensionless Time (\tau)');
ylabel('Displacement x(\tau)');
legend([h(1) h_env h(2) h(3)], ...
    'Under (\zeta < 1.0)', ...
    'Envelope of UnderDamping', ...
    'Critical (\zeta = 1.0)', ...
    'Over (\zeta > 1.0)')
disp('Underdamping solution:')
disp(sol(1))
disp('Critical damping solution:')
disp(sol(2))
disp('Overdamping solution:')
disp(sol(3))

%% Comparing same zeta, but different ω0 & γ system  (Underdamping case)
w0_list = [5.0, 1.0];
gamma_list = [1.0, 0.2];
zeta_same = gamma_list./w0_list;
Linewidth_list = [5.0, 0.5];
color_tau_list = ['y','k'];
syms x(t) gamma w0 x_nd(tau) zeta
t_plot = linspace(0,20,2000);
tau_plot = linspace(0,20,2000);

ode_general_tau = diff(x_nd,tau,2) + 2*zeta*diff(x_nd,tau) + x_nd == 0;
cond1_tau = x_nd(0) == 1;
cond2_tau = subs(diff(x_nd,tau), tau, 0) == 0;

ode_general_t = diff(x,t,2) + 2*gamma*diff(x,t) + w0^2*x == 0;
cond1_t = x(0) == 1;
cond2_t = subs(diff(x,t), t, 0) == 0;

fig2 = figure;
subplot(2,1,1)
hold on

for j = 1:length(gamma_list)
    ode_t = subs(ode_general_t,{gamma,w0},{gamma_list(j),w0_list(j)});
    sol_t(j) = dsolve(ode_t,cond1_t,cond2_t);
    fun_t = matlabFunction(sol_t(j));
    h_t(j) = plot(t_plot, fun_t(t_plot), color(j), 'LineWidth', 1.5);
end
title('Comparing same \zeta, different \gamma, \omega0 system - t-plot')
xlabel('Time')
ylabel('Displacement x(t)')
legend([h_t(1), h_t(2)], ...
    '\omega0 = 5.0, \gamma = 1.0', ...
    '\omega0 = 1.0, \gamma = 0.2')

disp('----------------------------------------------------')
disp('\omega0 = 5.0, \gamma = 1.0 system solution - x(t)')
disp(sol_t(1))
disp('\omega0 = 1.0, \gamma = 0.2 system solution - x(t)')
disp(sol_t(2))


subplot(2,1,2)
hold on
for k = 1:length(zeta_same)
    ode_tau = subs(ode_general_tau,{zeta},{zeta_same(k)});
    sol_tau(k) = dsolve(ode_tau,cond1_tau,cond2_tau);
    fun_tau = matlabFunction(sol_tau(k));
    h_tau(k) = plot(tau_plot, fun_tau(tau_plot), color_tau_list(k), 'LineWidth', Linewidth_list(k));
end
title('Comparing same \zeta, different \gamma, \omega0 system - tau-plot')
xlabel('Dimensionless Time (\tau)')
ylabel('Displacement x(\tau)')
legend([h_tau(1), h_tau(2)], ...
    '\omega0 = 5.0, \gamma = 1.0', ...
    '\omega0 = 1.0, \gamma = 0.2')

disp('----------------------------------------------------')
disp('\omega0 = 5.0, \gamma = 1.0 system solution - x(\tau)')
disp(sol_tau(1))
disp('\omega0 = 1.0, \gamma = 0.2 system solution - x(\tau)')
disp(sol_tau(2))

%% Envelope Analysis (Underdamped Case)
tau = linspace(0,40,4000);

zeta_list = [0.05 0.2 0.5];   % 여러 감쇠비 비교
colors = ['b','r','g'];

fig3 = figure;
hold on

for m = 1:length(zeta_list)
    zeta = zeta_list(m);
    env = exp(-zeta*tau);
    h_env_a(m) = plot(tau,env,['--',colors(m)]);
    tau_decay = 1/zeta;
    xline(tau_decay,'--', ...
        ['\tau_{decay} = ' num2str(tau_decay)], ...
        'LabelVerticalAlignment','bottom');
end
yline(exp(-1),'--')
title('Envelope Analysis')
xlabel('Dimensionless Time (\tau)')
ylabel('Envelope(exp(-\zeta\tau))')
legend([h_env_a(1), h_env_a(2), h_env_a(3)], ...
    '\zeta = 0.05', ...
    '\zeta = 0.2', ...
    '\zeta = 0.5')

if ~exist('figures','dir')
    mkdir('figures')
end
saveas(fig1,'figures/Free_damped_Non_demensionalization.png')
saveas(fig2,'figures/Free_damped_Same_zeta_different_gamma&omega_system.png')
saveas(fig3,'figures/Free_damped_Envelopde_analysis.png')