
clc
clear
close all
rng(123);

set(0,'defaultAxesFontSize', 12)
set(0, 'DefaultLineLineWidth', 2);

%% Define a grid prior f(\pi)

grid_res = 0.1; % grid resolution
values_pi = [0:grid_res:1]; % values of \pi (probability of head)
N_pi = size(values_pi, 2); % numeber of \pi values to test

% prior probability on \theta 
f_theta_first_half = [0:grid_res:0.5]; 
f_theta_second_half = flip([0:grid_res:0.4]);
f_theta = [f_theta_first_half f_theta_second_half]; 
f_theta = f_theta./sum(f_theta); % make distribution sum to 1

sum(f_theta) % check if prior sums to one

%% Observed data and compute likelihood f(Y | \pi)

N = 10; % number of observed data
N_head = 7; % number of successes (should be < N)
N_cross = N - N_head; % number of insuccesses

Y = [ones(1, N_head) zeros(1, N_cross)]'; % data

f_Y_given_theta = ones(1, N_pi); % init likelihood

for pp = 1 : 1 : N_pi % for each value of \pi 
    for ii = 1 : 1 : N % for each datum
        
        % Bernoulli likelihood
        f_Y_ii = values_pi(pp)^Y(ii) * (1-values_pi(pp))^(1-Y(ii));
        f_Y_given_theta(pp) = f_Y_given_theta(pp) * f_Y_ii; % compute likelihood for pp-th parameter value
        
    end
end

% Note: the likelihood is not a probability distribution, so it does not
% sums to one
sum(f_Y_given_theta)

%% Compute posterior f(\theta | Y)

f_Y = sum(f_Y_given_theta .* f_theta); % marginal likelihood f(Y)

% Apply Bayes theorem
f_theta_given_Y = (f_Y_given_theta .* f_theta) ./ f_Y; % Posterior f(\theta | Y)

sum(f_theta_given_Y) % check if posterior sums to one

%% Plot single

figure
subplot(3,1,1)
stem(values_pi, f_theta, 'b','LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{\theta}(\pi)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); ylim([0, 0.35]);
xticks(values_pi);
title('$\textbf{Prior distribution}$', 'interpreter', 'latex');


figure
stem(values_pi, f_Y_given_theta, 'b','LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{Y\vert\theta}(Y\vert\pi)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); 
ylim([0, 0.003]);
xticks(values_pi);
title('$\textbf{Likelihood function}$', 'interpreter', 'latex');


figure
stem(values_pi, f_theta_given_Y, 'b','LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{\theta\vert Y}(\pi\vert Y)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); 
ylim([0, 0.35]);
xticks(values_pi);
title('$\textbf{Posterior distribution}$', 'interpreter', 'latex');

%% Plot all together

figure
subplot(3,1,1)
stem(values_pi, f_theta, 'b','LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{\theta}(\pi)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); ylim([0, 0.35]);
xticks(values_pi);
title('$\textbf{Prior distribution}$', 'interpreter', 'latex');

 
subplot(3,1,2)
stem(values_pi, f_Y_given_theta, 'r','LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{Y\vert\theta}(Y\vert\pi)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); 
ylim([0, 0.003]);
xticks(values_pi);
title('$\textbf{Likelihood function}$', 'interpreter', 'latex');

 
subplot(3,1,3)
stem(values_pi, f_theta_given_Y, 'Color', [128 0 128]./255, 'LineWidth', 2);
xlabel('$\pi$', 'interpreter', 'latex')
ylabel('$f_{\theta\vert Y}(\pi\vert Y)$', 'interpreter', 'latex')
grid on; xlim([-0.1, 1.1]); 
ylim([0, 0.35]);
xticks(values_pi);
title('$\textbf{Posterior distribution}$', 'interpreter', 'latex');
