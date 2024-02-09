clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform CCA for brain-behavior analysis
%
%  Jin
%  11/20/2023
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% loading data
% setting output folder
output_path = fullfile(working_path,'outputs','Stanford',filesep);
mkdir(output_path);

% loading the behavioral data file
numopt_mathres_N219 = importdata(fullfile(working_path,'data','Stanford cohort','numopt_mathres_N219.mat'));

% imaging data
input_path = fullfile(working_path,'data','Stanford cohort',filesep); % 
X = importdata(fullfile(input_path,'Stanford_GMV_BN246_N219.mat'));
Y = numopt_mathres_N219;

% Number of folds for cross-validation
numFolds = 10;

% Number of times to run the cross-validation
numRuns = 100;

% Number of permutations for significance testing
numPermutations = 1000;

originalMeanCorrelation = performCrossValidation(X, Y, numFolds, numRuns);
disp(originalMeanCorrelation);

permutedCorrelations = zeros(numPermutations, 1);
for p = 1:numPermutations
    disp(p)
    Y_permuted = Y(randperm(size(Y, 1)), :);
    permutedCorrelations(p) = performCrossValidation(X, Y_permuted, numFolds, numRuns);
    disp(permutedCorrelations(p))
end

% Calculate the p-value
pValue = mean(originalMeanCorrelation <= permutedCorrelations);

% Display the results
disp('Original Mean Test Correlation:');
disp(originalMeanCorrelation);
disp('P-Value from Permutation Test:');
disp(pValue);

function meanCorrelation = performCrossValidation(X, Y, numFolds, numRuns)
    allTestCorrelations = zeros(numRuns, numFolds);

    for run = 1:numRuns
        cv = cvpartition(size(X, 1), 'KFold', numFolds);

        for i = 1:numFolds
            X_train = X(training(cv, i), :);
            Y_train = Y(training(cv, i), :);
            X_test = X(test(cv, i), :);
            Y_test = Y(test(cv, i), :);

            [coeff, score, ~, ~, explained] = pca(X_train);
            
            % Determine the number of components to retain
            cumulativeVariance = cumsum(explained);
            numComponents = find(cumulativeVariance >= 85, 1, 'first');

            X_train_pca = score(:, 1:numComponents);
            X_test_pca = (X_test - mean(X_train)) * coeff(:, 1:numComponents);

            [A, B, r, ~, ~] = canoncorr(X_train_pca, Y_train);

            U_test = X_test_pca * A;
            V_test = Y_test * B;

            for j = 1:length(r)
                corr_r(j) = corr(U_test(:, j), V_test(:, j));
            end
            allTestCorrelations(run, i) = mean(corr_r);
        end
    end

    meanCorrelation = mean(allTestCorrelations, 'all');
end

