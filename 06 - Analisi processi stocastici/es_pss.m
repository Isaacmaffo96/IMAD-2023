
clc
clear
close all
rng('default');

set(0, 'defaultAxesFontSize', 12)
set(0, 'DefaultLineLineWidth', 2);
set(0, 'defaultAxesFontSize', 14)
set(0, 'defaultAxesTickLabelInterpreter','latex'); 
set(0, 'defaultlegendInterpreter','latex')


%% Load data 

load('data.mat');


%% Pre-process data

%% Qual'è la varianza del processo?

Ts = 1;
fs = 1/Ts; % sampling frequency [Hz]



lags = 30; % lags up to which compute the estimated autocovariance function
[Ryy_hat, ~, bounds] = autocorr(y, lags); % sample autocorrelation function
% - bounds contains the 95% confidence intervals on the estimate

gammayy_hat_0 = var(y); % variance of y(t)
gammayy_hat = gammayy_hat_0 * Ryy_hat; % autocovariance function

figure
hold on
stem(0:1:lags, gammayy_hat, 'b', 'linewidth', 2); 
l = line([0 lags], [bounds(1) bounds(1)]); l.LineStyle = '-'; l.LineWidth = 1.5; l.Color = 'r';
l = line([0 lags], [bounds(2) bounds(2)]); l.LineStyle = '-'; l.LineWidth = 1.5; l.Color = 'r';
grid on; xlabel('Lags', 'interpreter', 'latex');  ylabel('$\gamma_{yy}(\tau)$', 'interpreter', 'latex')
title('Sample autocovariance', 'interpreter', 'latex'); 
xticks(0:1:lags);

% Partial autocorrelation function
[gammayy_PAR_hat, ~, bounds] = parcorr(y, lags);

figure
hold on
stem(0:1:lags, gammayy_PAR_hat, 'b', 'linewidth', 2); 
l = line([0 lags], [bounds(1) bounds(1)]); l.LineStyle = '-'; l.LineWidth = 1.5; l.Color = 'r';
l = line([0 lags], [bounds(2) bounds(2)]); l.LineStyle = '-'; l.LineWidth = 1.5; l.Color = 'r';
grid on; xlabel('Lags', 'interpreter', 'latex'); ylabel('$\gamma^\textrm{PAR}_{yy}(\tau)$', 'interpreter', 'latex')
title('Sample partial autocorrelation', 'interpreter', 'latex'); 
xticks(0:1:lags);

%% Il processo è modellabile come un AR, MA o ARMA?

%% Che forma ha la densità spettrale di potenza?

N = size(y, 1);

% Estimation with the fft
bin = fs/N; % frequency resolution of the fft
freqs = [0 bin:bin:(N/2-1)*bin]; % grid of discrete frequencies

Gamma_yy_hat_fft = 1/N * abs(fft(y)).^2; % alternative command
Gamma_yy_hat_fft_useful = Gamma_yy_hat_fft(1:floor(N/2)); % retain only useful data from the fft


% Smoothed periodogram (10 folds with 1000 elements each)
n_folds  = 300;               % folds number
M = round(N/n_folds);        % elements in each fold
bin_smooth = fs/M;           % frequency resolution for smoothed periodogram
freqs_smooth = [0 bin_smooth:bin_smooth:(M/2-1)*bin_smooth]; % grid of discrete frequencies for smoothed periodogram

Gamma_yy_hat_temp = zeros(M, n_folds);
for k = 1 : n_folds
    temp = fft(y(1+(k-1)*M:M+(k-1)*M));
    Gamma_yy_hat_temp(:, k) = abs(temp).^2 / M;
end
Gamma_yy_hat_smooth = mean(Gamma_yy_hat_temp, 2);  % estimate as the mean of the single spectra
Gamma_yy_hat_smooth_useful = Gamma_yy_hat_smooth(1:floor(M/2));


 

% Estimation with pwelch
w_norm_eval = 0:0.1:pi;
WINDOW = 20;
Gamma_yy_hat_pwelch = pwelch(y, WINDOW, [], w_norm_eval);
Gamma_yy_hat_pwelch = 2*pi * Gamma_yy_hat_pwelch; % the pwelch method returns a scaled psd, see documentation at https://it.mathworks.com/help/signal/ref/pwelch.html

figure;
plot(freqs, Gamma_yy_hat_fft_useful, 'r:', 'linewidth', 1);
hold on
plot(freqs_smooth, Gamma_yy_hat_smooth_useful, 'g', 'linewidth', 1.5); 
plot(w_norm_eval/2/pi, Gamma_yy_hat_pwelch , 'b', 'linewidth', 1.5); 
xlabel('Frequency [Hz]', 'interpreter', 'latex'); ylabel('Amplitude', 'interpreter', 'latex');
grid on;
xlim([0, fs/2]);
legend('Estimated PSD with $\texttt{fft}$', ...
       'Estimated PSD with $\texttt{fft}$ smoothed', ...
       'Estimated PSD with $\texttt{pwelch}$');


