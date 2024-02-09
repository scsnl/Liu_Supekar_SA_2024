clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform prediction of another math measure based 
%  on CCA of brain-behavior analysis
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
output_path = fullfile(working_path,'outputs','Stanford',filesep);

% loading behavioral data 
Applied_N110 = importdata(fullfile(working_path,'data','Stanford cohort','validation_AppliedRes_subN110.mat'));

%% CCA-based prediction
% loading the gray matter volume data
Stanford_GMVsub110_BN246=importdata(fullfile(working_path,'data','Stanford cohort','Stanford_GMV_BN246_subN110_appliedres.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(output_path,'cca_output_pca85.mat'));
numopt_mathres_N219 = importdata(fullfile(working_path,'data','Stanford cohort','numopt_mathres_N219.mat'));

% number of subject
N = length(Applied_N110);

% scale data using PCA from stanford
X = Stanford_GMVsub110_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% define Y by filling the missing math reasoning data
Y = [Applied_N110 zeros(110,1)];

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A
% V = (Y - repmat(mean(Y),N,1))*B

% scale U to get the predict Y
pred_Y = (U(:,1).*std(numopt_mathres_N219(:,1)) + mean(numopt_mathres_N219(:,1)))

% corr between predict Y and real Y
[r1 p]=corr(pred_Y(:,1),Y(:,1))
% [r p]=corr(U(:,1),V(:,1))

%% figure
output_tiff=fullfile(output_path,'CCA_based_wjiii_prediction.tif');
beh_corr_scatter_NKI(Y(:,1),U(:,1),'Applied Reasoning',['r = ' num2str(r1(1),2)],['p = ' num2str(p(1),3)],'MAIP brain score','Prediction',output_tiff)

