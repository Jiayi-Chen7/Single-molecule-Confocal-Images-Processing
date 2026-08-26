function [mu1, sigma1, mu2, sigma2, pi1, pi2, gamma1, threshold] = em_gmm_length_analysis(data, varargin)
% EM_ALGORITHM_FOR_TWO_COMPONENT_GAUSSIAN_MIXTURE_MODEL 
% Estimates parameters of a 2-component GMM using EM algorithm with flexible initialization
%
% Inputs:
%   data       - Observed length measurements (N×1 vector)
% Optional Parameters:
%   'MaxIter'  - Maximum iterations (default: 100)
%   'Tol'      - Convergence tolerance (default: 1e-6)
%   'InitMethod' - Initialization method: 
%                 'kmeans'(default), 'random', 'percentile', 'manual'
%   'InitParams' - Manual initialization parameters [mu1, sigma1, mu2, sigma2, pi1]
%
% Outputs:
%   mu1, sigma1 - Mean and std of component 1
%   mu2, sigma2 - Mean and std of component 2  
%   pi1, pi2    - Mixing proportions
%   gamma1      - Posterior probabilities for component 1
%   threshold   - Optimal classification threshold

%% Input Validation and Parameter Parsing
p = inputParser;
addRequired(p, 'data', @(x) validateattributes(x, {'numeric'}, ...
    {'vector', 'real', 'finite'}));
addParameter(p, 'MaxIter', 100, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive', 'integer'}));
addParameter(p, 'Tol', 1e-6, @(x) validateattributes(x, {'numeric'}, ...
    {'scalar', 'positive'}));
addParameter(p, 'InitMethod', 'kmeans', @(x) ismember(x, ...
    {'kmeans', 'random', 'percentile', 'manual'}));
addParameter(p, 'InitParams', [], @(x) isempty(x) || ...
    (isnumeric(x) && numel(x)==5));
parse(p, data, varargin{:});

% Extract parsed parameters
max_iter = p.Results.MaxIter;
tol = p.Results.Tol;
init_method = p.Results.InitMethod;
init_params = p.Results.InitParams;

data = data(:); % Ensure column vector
N = length(data);

%% Initialization Methods
switch lower(init_method)
    case 'kmeans'
        % Method 1: k-means clustering initialization (default)
        [idx, centers] = kmeans(data, 2);
        mu1 = centers(1);
        mu2 = centers(2);
        sigma1 = std(data(idx == 1));
        sigma2 = std(data(idx == 2));
        pi1 = sum(idx == 1)/N;
        
    case 'random'
        % Method 2: Random initialization within data range
        mu1 = min(data) + (max(data)-min(data))*rand();
        mu2 = min(data) + (max(data)-min(data))*rand();
        sigma1 = std(data)*rand();
        sigma2 = std(data)*rand();
        pi1 = rand();
        
    case 'percentile'
        % Method 3: Percentile-based initialization
        p25 = prctile(data, 25);
        p75 = prctile(data, 75);
        mu1 = mean(data(data < p25));
        mu2 = mean(data(data > p75));
        sigma1 = std(data(data < p25));
        sigma2 = std(data(data > p75));
        pi1 = sum(data < p25)/N;
        
    case 'manual'
        % Method 4: User-provided parameters
        if isempty(init_params)
            error('InitParams must be provided when using manual initialization');
        end
        mu1 = init_params(1);
        sigma1 = init_params(2);
        mu2 = init_params(3);
        sigma2 = init_params(4);
        pi1 = init_params(5);
        
    otherwise
        error('Unsupported initialization method');
end

% Ensure parameter validity
sigma1 = max(sigma1, eps);  % Prevent zero standard deviation
sigma2 = max(sigma2, eps);
pi1 = max(min(pi1, 1), 0);  % Keep proportion in [0,1]
pi2 = 1 - pi1;

fprintf('Initialization method: %s\n', init_method);
fprintf('Initial params: μ1=%.4f, σ1=%.4f, μ2=%.4f, σ2=%.4f, π1=%.4f\n', ...
        mu1, sigma1, mu2, sigma2, pi1);

%% EM Algorithm Main Loop
log_likelihood = -inf;
converged = false;

for iter = 1:max_iter
    % E-step: Compute posterior probabilities
    prob1 = pi1 * normpdf(data, mu1, sigma1);
    prob2 = pi2 * normpdf(data, mu2, sigma2);
    total_prob = prob1 + prob2;
    total_prob(total_prob < eps) = eps; % Avoid division by zero
    
    gamma1 = prob1 ./ total_prob;
    gamma2 = 1 - gamma1;
    
    % Compute current log-likelihood
    current_log_likelihood = sum(log(total_prob));
    
    % Check convergence
    if abs(current_log_likelihood - log_likelihood) < tol
        converged = true;
        break;
    end
    log_likelihood = current_log_likelihood;
    
    % M-step: Update parameters
    N1 = sum(gamma1);
    N2 = N - N1;
    
    % Prevent component collapse
    if N1 < eps || N2 < eps
        warning('One component is collapsing - try different initialization');
        break;
    end
    
    % Update mixing proportions
    pi1 = N1 / N;
    pi2 = N2 / N;
    
    % Update component 1 parameters
    mu1 = sum(gamma1 .* data) / N1;
    sigma1 = sqrt(sum(gamma1 .* (data - mu1).^2) / N1);
    
    % Update component 2 parameters
    mu2 = sum(gamma2 .* data) / N2;
    sigma2 = sqrt(sum(gamma2 .* (data - mu2).^2) / N2);
end

%% Calculate Classification Threshold
[threshold, found] = find_gmm_threshold(mu1, sigma1, mu2, sigma2, pi1, pi2, min(data), max(data));
if ~found
    warning('Analytical solution not found - using approximate threshold');
end

%% Visualization
visualize_results(data, mu1, sigma1, mu2, sigma2, pi1, pi2, threshold);

%% Display Results
if converged
    fprintf('EM converged in %d iterations\n', iter);
else
    fprintf('EM stopped after %d iterations (may not have fully converged)\n', iter);
end
fprintf('Final parameters: μ1=%.4f, σ1=%.4f, μ2=%.4f, σ2=%.4f, π1=%.4f\n', ...
        mu1, sigma1, mu2, sigma2, pi1);
fprintf('Classification threshold=%.4f\n', threshold);
end

%% Helper Function: Find Optimal Threshold
function [threshold, found] = find_gmm_threshold(mu1, sigma1, mu2, sigma2, pi1, pi2, min_val, max_val)
% Attempt analytical solution using symbolic math
syms x
eqn = pi1 * normpdf(x, mu1, sigma1) == pi2 * normpdf(x, mu2, sigma2);
solutions = double(vpasolve(eqn, x, [min_val, max_val]));

if ~isempty(solutions)
    threshold = solutions(1);
    found = true;
else
    % Numerical approximation if no analytical solution
    x_vals = linspace(min_val, max_val, 1000);
    diff = pi1 * normpdf(x_vals, mu1, sigma1) - pi2 * normpdf(x_vals, mu2, sigma2);
    [~, idx] = min(abs(diff));
    threshold = x_vals(idx);
    found = false;
end
end

%% Helper Function: Visualization
function visualize_results(data, mu1, sigma1, mu2, sigma2, pi1, pi2, threshold)
% Create figure window
figure;
set(gcf, 'Position', [100, 100, 800, 600]);

% Subplot 1: Data histogram and fitted distributions
subplot(2,2,1);
h = 2.*iqr(data).*length(data).^(-1/3);
BinNum = round((max(data)-min(data))./h);
histogram(data,BinNum, 'Normalization', 'pdf', 'EdgeColor', 'none');
hold on;
x_range = linspace(min(data), max(data), 1000);
plot(x_range, pi1 * normpdf(x_range, mu1, sigma1), 'r', 'LineWidth', 2);
plot(x_range, pi2 * normpdf(x_range, mu2, sigma2), 'g', 'LineWidth', 2);
plot(x_range, pi1 * normpdf(x_range, mu1, sigma1) + pi2 * normpdf(x_range, mu2, sigma2), ...
     'b', 'LineWidth', 2);
xline(threshold, '--k', sprintf('Threshold=%.2f', threshold), 'LineWidth', 1.5);
legend('Data', 'Component 1', 'Component 2', 'Mixture', 'Location', 'best');
xlabel('Length');
ylabel('Probability Density');
title('EM Algorithm for Two-Component GMM');

% Subplot 2: Classification probabilities
subplot(2,2,2);
[~, sort_idx] = sort(data);
plot(data(sort_idx), pi1 * normpdf(data(sort_idx), mu1, sigma1) ./ ...
    (pi1 * normpdf(data(sort_idx), mu1, sigma1) + pi2 * normpdf(data(sort_idx), mu2, sigma2)), ...
    'LineWidth', 2);
hold on;
yline(0.5, '--r', 'Decision Boundary');
xline(threshold, '--k', sprintf('Threshold=%.2f', threshold), 'LineWidth', 1.5);
xlabel('Length');
ylabel('Posterior Probability (Component 1)');
title('Classification Probability');
grid on;
end