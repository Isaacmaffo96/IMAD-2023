
clc
clear
close all
set(0,'defaultAxesFontSize', 12)
set(0, 'DefaultLineLineWidth', 2);

%% Load Data
%  The first two columns contains the exam scores and the third column
%  contains the label.

data = load('Data\ex2data1.txt');
X = data(:, [1, 2]); % regressors matrix
y = data(:, 3); % output

%% ==================== Part 1: Plotting ====================

fprintf(['Plotting data with + indicating (y = 1) examples and o ' ...
         'indicating (y = 0) examples.\n']);

plotData(X, y);
hold on;
xlabel('Exam 1 score')
ylabel('Exam 2 score')
grid on
legend('Admitted', 'Not admitted')

%% ============ Part 2: Compute Cost and Gradient ============

%  Setup the data matrix appropriately, and add ones for the intercept term
[N, d] = size(X);

X = [ones(N, 1) X]; % Add intercept term 

initial_theta = zeros(d + 1, 1); % Initialize fitting parameters

% Compute and display initial cost and gradient
[cost, grad] = costFunction(initial_theta, X, y);


fprintf('Cost at initial theta (zeros): %f\n', cost);
fprintf('Gradient at initial theta (zeros): \n');
fprintf(' %f \n', grad);


%% ============= Part 3: Optimizing using fminunc  =============

% fminunc is a nonlinear optimization function in Matlab that attempts to 
% find the minimum of a single or multi-variable objective function.
% It uses a gradient descent algorithm to search for the minimum of the 
% function, starting from an initial user-specified starting point. 
% fminunc attempts to find the global minimum of the function, although 
% it may only find a local minimum depending on the starting point.

%  Set options for fminunc
options = optimset('GradObj', 'on', 'MaxIter', 400);

%  Run fminunc to obtain the optimal theta
%  This function will return theta and the cost 
[theta, cost] = fminunc(@(t)(costFunction(t, X, y)), initial_theta, options);

% Print theta to screen
fprintf('Cost at theta found by fminunc: %f\n', cost);
fprintf('theta: \n');
fprintf(' %f \n', theta);

% Plot Boundary
plotDecisionBoundary(theta, X, y, []);

% Put some labels 
hold on;
% Labels and Legend
xlabel('Exam 1 score')
ylabel('Exam 2 score')

% Specified in plot order
legend('Admitted', 'Not admitted')
hold off;


%% ============== Part 4: Predict and Accuracies ==============
%  After learning the parameters, you'll like to use it to predict the outcomes
%  on unseen data. In this part, you will use the logistic regression model
%  to predict the probability that a student with score 45 on exam 1 and 
%  score 85 on exam 2 will be admitted.
%
%  Furthermore, you will compute the training and test set accuracies of 
%  our model.
%

%  Predict probability for a student with score 45 on exam 1 
%  and score 85 on exam 2 

prob = sigmoid([1 45 85] * theta);
fprintf(['For a student with scores 45 and 85, we predict an admission ' ...
         'probability of %f\n\n'], prob);

% Compute accuracy on our training set
% make the prediction using the model on the training data set
p = predict(theta, X);

fprintf('Train Accuracy: %f\n', mean(p == y) * 100);

