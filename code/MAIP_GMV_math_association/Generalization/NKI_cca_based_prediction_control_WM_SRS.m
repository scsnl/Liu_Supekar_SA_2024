clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform CCA for brain-behavior analysis
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

%% CCA-based prediction
% loading the gray matter volume data from NKI
NKI_GMVsub91_BN246=importdata(fullfile(working_path,'data','NKI-RS','NKI_GMV_BN246_N91.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(working_path,'outputs','Stanford','cca_output_pca85.mat'));

% number of subject
N = 91;

% scale data using PCA from stanford
X = NKI_GMVsub91_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A;

%% loading data
% setting output path
output_path = fullfile(working_path,'outputs','NKI');

% SRS
NKI_SRS_N91 = importdata(fullfile(working_path,'data','NKI-RS','SRS_N91.mat'));

Y = [NKI_SRS_N91 zeros(N,1)];
Y(12,:)=[]; % invalid data with zero
U1=U(:,1);U1(12,:)=[]; % invalid data with zero

[r p]=corr(U1(:,1),Y(:,1))

%% WM
WM_score = importdata(fullfile(working_path,'data','NKI-RS','WM_N91.mat'))

% excluded NaN
valid_ind = find(~isnan(WM_score))
Y = [WM_score(valid_ind)];

U2 = U(valid_ind,1)
[r p]=corr(U2(:,1),Y(:,1))
