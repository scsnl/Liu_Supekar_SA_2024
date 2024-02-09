clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to validate the specificity of CCA mode
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
WM_N94 = importdata(fullfile(working_path,'data','Stanford cohort','validation_WM_subN94.mat'));

%% CCA-based prediction for WM
% loading the gray matter volume data
Stanford_GMVsub94_BN246=importdata(fullfile(working_path,'data','Stanford cohort','Stanford_GMV_BN246_subN94_WM.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(output_path,'cca_output_pca85.mat'));

% number of subject
N = length(WM_N94);

% scale data using PCA from stanford
X = Stanford_GMVsub94_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A
Y = [WM_N94 zeros(N,1)]
V = (Y - repmat(mean(Y),N,1))*B

% corr between predict Y and real Y
[r1 p]=corr(U(:,1),Y(:,1))

%% figure
output_tiff=fullfile(output_path,'CCA_based_WMprediction.tif');
beh_corr_scatter_NKI(Y(:,1),U(:,1),'Working Memory',['r = ' num2str(r1(1),2)],['p = ' num2str(p(1),3)],'MAIP brain score','Prediction',output_tiff)

