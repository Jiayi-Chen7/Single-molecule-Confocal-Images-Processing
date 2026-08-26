function [mu_log, sigma_log, mu_gauss, sigma_gauss, pi_log, pi_gauss, gamma_log, threshold] = mixed_em_gmm_ori(data, varargin)
% EM_ALGORITHM_FOR_LOGNORMAL_GAUSSIAN_MIXTURE
% Estimates parameters of a mixture with one lognormal and one Gaussian component
%
% Inputs:
%   data       - Observed measurements (positive values for lognormal part)
% Optional Parameters:
%   'MaxIter'  - Maximum iterations (default: 100)
%   'Tol'      - Convergence tolerance (default: 1e-6)
%   'InitMethod' - Initialization method: 'kmeans'(default), 'percentile', 'manual'
%   'InitParams' - Manual initialization [mu_log, sigma_log, mu_gauss, sigma_gauss, pi_log]
%
% Outputs:
%   Parameters for both components and classification threshold

%% Input Validation
validateattributes(data, {'numeric'}, {'vector', 'real', 'finite'});
if any(data <= 0)
    warning('Negative values will be assigned to Gaussian component');
end

% Parse optional parameters
p = inputParser;
addParameter(p, 'MaxIter', 100, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive', 'integer'}));
addParameter(p, 'Tol', 1e-6, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive'}));
addParameter(p, 'InitMethod', 'kmeans', @(x) ismember(x, ...
    {'kmeans', 'percentile', 'manual'}));
addParameter(p, 'InitParams', [], @(x) isempty(x) || ...
    (isnumeric(x) && numel(x)==5));
parse(p, varargin{:});

max_iter = p.Results.MaxIter;
tol = p.Results.Tol;
init_method = p.Results.InitMethod;
init_params = p.Results.InitParams;

N = length(data);

%% Initialization
switch lower(init_method)
    case 'kmeans'
        % k-means initialization
        [idx, centers] = kmeans(data, 2);
        
        % Assign component types based on skewness
        skew1 = skewness(data(idx == 1));
        skew2 = skewness(data(idx == 2));
        
        if skew1 > skew2
            log_idx = 1;
            gauss_idx = 2;
        else
            log_idx = 2;
            gauss_idx = 1;
        end
        
        % Lognormal component (log-space params)
        log_data = log(data(idx == log_idx));
        mu_log = mean(log_data);
        sigma_log = std(log_data);
        
        % Gaussian component
        mu_gauss = centers(gauss_idx);
        sigma_gauss = std(data(idx == gauss_idx));
        pi_log = sum(idx == log_idx)/N;
        
    case 'percentile'
        % Percentile-based initialization
        p30 = prctile(data, 30);
        p70 = prctile(data, 70);
        
        % Assume lower values are lognormal
        log_data = log(data(data < p30));
        mu_log = mean(log_data);
        sigma_log = std(log_data);
        
        % Higher values are Gaussian
        gauss_data = data(data > p70);
        mu_gauss = mean(gauss_data);
        sigma_gauss = std(gauss_data);
        pi_log = sum(data < p30)/N;
        
    case 'manual'
        % User-provided parameters
        if isempty(init_params)
            error('InitParams required for manual initialization');
        end
        mu_log = init_params(1);
        sigma_log = init_params(2);
        mu_gauss = init_params(3);
        sigma_gauss = init_params(4);
        pi_log = init_params(5);
end

% Parameter validation
sigma_log = max(sigma_log, eps);
sigma_gauss = max(sigma_gauss, eps);
pi_log = max(min(pi_log, 1), 0);
pi_gauss = 1 - pi_log;

fprintf('Initial parameters:\n');
fprintf('Lognormal: μ_log=%.4f, σ_log=%.4f\n', mu_log, sigma_log);
fprintf('Gaussian: μ=%.4f, σ=%.4f\n', mu_gauss, sigma_gauss);
fprintf('Mixing proportions: π_log=%.4f, π_gauss=%.4f\n', pi_log, pi_gauss);

%% EM Algorithm
log_likelihood = -inf;
converged = false;

for iter = 1:max_iter
    % E-step: Compute component responsibilities
    log_prob = zeros(N,1);
    gauss_prob = zeros(N,1);
    
    % Lognormal component (only for positive data)
    pos_idx = data > 0;
    log_prob(pos_idx) = pi_log .* (1./(data(pos_idx)*sigma_log*sqrt(2*pi))) .* ...
                        exp(-(log(data(pos_idx))-mu_log).^2./(2*sigma_log^2));
    
    % Gaussian component (all data)
    gauss_prob = pi_gauss .* normpdf(data, mu_gauss, sigma_gauss);
    
    total_prob = log_prob + gauss_prob;
    total_prob(total_prob < eps) = eps; % Avoid division by zero
    
    gamma_log = log_prob ./ total_prob;
    gamma_gauss = 1 - gamma_log;
    
    % Log-likelihood
    current_log_likelihood = sum(log(total_prob));
    
    % Check convergence
    if abs(current_log_likelihood - log_likelihood) < tol
        converged = true;
        break;
    end
    log_likelihood = current_log_likelihood;
    
    % M-step: Update parameters
    N_log = sum(gamma_log);
    N_gauss = sum(gamma_gauss);
    
    % Prevent component collapse
    if N_log < eps || N_gauss < eps
        warning('Component collapse detected');
        break;
    end
    
    % Update mixing proportions
    pi_log = N_log / N;
    pi_gauss = N_gauss / N;
    
    % Update lognormal parameters (only using positive data)
    if N_log > 0
        weighted_log = sum(gamma_log(pos_idx) .* log(data(pos_idx))) / N_log;
        mu_log = weighted_log;
        sigma_log = sqrt(sum(gamma_log(pos_idx) .* (log(data(pos_idx)) - mu_log).^2) / N_log);
    end
    
    % Update Gaussian parameters
    mu_gauss = sum(gamma_gauss .* data) / N_gauss;
    sigma_gauss = sqrt(sum(gamma_gauss .* (data - mu_gauss).^2) / N_gauss);
end

%% Calculate Classification Threshold
% Find where densities are equal
[threshold, found] = find_mixed_threshold(mu_log, sigma_log, mu_gauss, sigma_gauss, ...
                                         pi_log, pi_gauss, min(data), max(data));
if ~found
    warning('Using approximate threshold');
end

%% Visualization
visualize_mixed_results(data, mu_log, sigma_log, mu_gauss, sigma_gauss, pi_log, pi_gauss, threshold);

%% Display Results
fprintf('\n=== EM Algorithm Results ===\n');
if converged
    fprintf('Converged in %d iterations\n', iter);
else
    fprintf('Stopped after %d iterations\n', iter);
end

fprintf('\nLognormal Component (log-space):\n');
fprintf('μ=%.4f, σ=%.4f\n', mu_log, sigma_log);
fprintf('Original scale median=%.4f, geometric std=%.4f\n', exp(mu_log), exp(sigma_log));

fprintf('\nGaussian Component:\n');
fprintf('μ=%.4f, σ=%.4f\n', mu_gauss, sigma_gauss);

fprintf('\nMixing proportions:\n');
fprintf('π_lognormal=%.4f, π_gaussian=%.4f\n', pi_log, pi_gauss);
fprintf('Classification threshold: %.4f\n', threshold);
end

%% Helper Functions
function [threshold, found] = find_mixed_threshold(mu_log, sigma_log, mu_gauss, sigma_gauss, ...
                                                  pi_log, pi_gauss, min_val, max_val)
% Find where lognormal and Gaussian densities are equal
syms x
eqn = pi_log*(1/(x*sigma_log*sqrt(2*pi)))*exp(-(log(x)-mu_log)^2/(2*sigma_log^2)) == ...
      pi_gauss*normpdf(x, mu_gauss, sigma_gauss);
  
try
    threshold = double(vpasolve(eqn, x, [max(min_val,eps), max_val]));
    found = ~isempty(threshold);
    if ~found
        error('No symbolic solution');
    end
catch
    % Numerical approximation
    x_vals = linspace(max(min_val,eps), max_val, 1000);
    lognorm_dens = pi_log*(1./(x_vals*sigma_log*sqrt(2*pi))).*exp(-(log(x_vals)-mu_log).^2./(2*sigma_log^2));
    gauss_dens = pi_gauss*normpdf(x_vals, mu_gauss, sigma_gauss);
    [~, idx] = min(abs(lognorm_dens - gauss_dens));
    threshold = x_vals(idx);
    found = false;
end
end

function visualize_mixed_results(data, mu_log, sigma_log, mu_gauss, sigma_gauss, pi_log, pi_gauss, threshold)
figure('Position', [100, 100, 1000, 800]);

% Panel 1: Histogram with fitted densities
subplot(2,2,1);
h = 2.*iqr(data).*length(data).^(-1/3);
BinNum = round((max(data)-min(data))./h);
histogram(data,BinNum,'Normalization', 'pdf', 'EdgeColor', 'none');
hold on;
x_vals = linspace(max(min(data),eps), max(data), 1000);

% Lognormal density
lognorm_dens = pi_log*(1./(x_vals*sigma_log*sqrt(2*pi))).*exp(-(log(x_vals)-mu_log).^2./(2*sigma_log^2));
plot(x_vals, lognorm_dens, 'r', 'LineWidth', 2);

% Gaussian density
gauss_dens = pi_gauss*normpdf(x_vals, mu_gauss, sigma_gauss);
plot(x_vals, gauss_dens, 'g', 'LineWidth', 2);

% Mixture density
plot(x_vals, lognorm_dens + gauss_dens, 'b', 'LineWidth', 2);

xline(threshold, '--k', sprintf('Threshold=%.2f', threshold), 'LineWidth', 1.5);
xlabel('Value');
ylabel('Probability Density');
title('Mixture Model Fit');
legend('Data', 'Lognormal', 'Gaussian', 'Mixture', 'Location', 'best');
grid on;

% Panel 2: QQ plots
subplot(2,2,2);
probplot('lognormal', data(data > 0));
hold on;
probplot('normal', data);
title('QQ Plots for Distribution Checking');
legend('Lognormal QQ', 'Normal QQ', 'Location', 'best');
grid on;

% Panel 3: Cumulative distribution
subplot(2,2,3);
ecdf(data);
hold on;
plot(x_vals, pi_log*logncdf(x_vals, mu_log, sigma_log) + pi_gauss*normcdf(x_vals, mu_gauss, sigma_gauss), ...
     'r', 'LineWidth', 2);
title('Empirical vs. Fitted CDF');
legend('Empirical', 'Fitted', 'Location', 'best');
grid on;

% Panel 4: Posterior probabilities
subplot(2,2,4);
[~, sort_idx] = sort(data);
posterior_log = zeros(size(data));
posterior_log(data > 0) = pi_log*(1./(data(data>0)*sigma_log*sqrt(2*pi))).*...
                          exp(-(log(data(data>0))-mu_log).^2./(2*sigma_log^2));
posterior_gauss = pi_gauss*normpdf(data, mu_gauss, sigma_gauss);
posterior = posterior_log./(posterior_log + posterior_gauss);

plot(data(sort_idx), posterior(sort_idx), 'LineWidth', 2);
hold on;
yline(0.5, '--r', 'Decision Boundary');
xline(threshold, '--k', 'Threshold', 'LineWidth', 1.5);
xlabel('Value');
ylabel('P(Lognormal | Data)');
title('Component Membership Probability');
grid on;
end