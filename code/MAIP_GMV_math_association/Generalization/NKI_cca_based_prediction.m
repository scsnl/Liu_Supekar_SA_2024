clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform prediction of NKI math scores based on
%  CCA of brain-behavior analysis from Stanford cohort
%
%  Jin
%  7/25/2022  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% loading data
% setting output path
output_path = fullfile(working_path,'outputs','NKI');
mkdir(output_path);

% loading the numopt data from NKI
numopt_sub91=importdata(fullfile(working_path,'data','NKI-RS','numopt_N91.mat'));

%% CCA-based prediction
% loading the gray matter volume data from NKI
NKI_GMVsub91_BN246=importdata(fullfile(working_path,'data','NKI-RS','NKI_GMV_BN246_N91.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(working_path,'outputs','Stanford','cca_output_pca85.mat'));
numopt_mathres_N219 = importdata(fullfile(working_path,'data','Stanford Cohort','numopt_mathres_N219.mat'));

% number of subject
N = length(numopt_sub91);

% scale data using PCA from stanford
X = NKI_GMVsub91_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% define Y by filling the missing math reasoning data
Y = [numopt_sub91 zeros(N,1)];

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A;
V = (Y - repmat(mean(Y),N,1))*B;

% scale U to get the predict Y
pred_Y = (U(:,1).*std(numopt_mathres_N219(:,1)) + mean(numopt_mathres_N219(:,1)));

% corr between predict Y and real Y
[r p]=corr(pred_Y(:,1),Y(:,1))
% [r p]=corr(U(:,1),V(:,1))

%% figure
output_tiff=fullfile(output_path,'CCA_based_prediction.tif');
beh_corr_scatter_NKI(Y(:,1),pred_Y(:,1),'Observed scores',['r = ' num2str(r(1),2)],['p = ' num2str(p(1),3)],'Predicted scores','Prediction Numerical Operations',output_tiff)
