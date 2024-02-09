clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform CCA validation 
%
%  Jin
%  8/10/2022  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% loading data
output_path = fullfile(working_path,'outputs','Stanford',filesep);
numopt_mathres_N219 = importdata(fullfile(working_path,'data','Stanford cohort','numopt_mathres_N219.mat'));

input_path = fullfile(working_path,'data','Stanford cohort',filesep); % 
img_path = fullfile(input_path,'Stanford_GMV_BN246_N219.mat'); % the path of imaging features with a matrix of subject x features
component_num = 1; % only check the first CCA component
pca_ind = 1; % PCA will be used before CCA when pca_ind is set to 1
perm_time = 1000; % permuation test 
beh_path = fullfile(working_path,'data','Stanford cohort','numopt_mathres_N219.mat'); % the path of behavior features with a matrix of subject x features
perm_pfwer = 0; % set to 1 if used pfwer permutation approach
pca_percentage = 0.85; % set a threshold such as 0.9 or 0.8 if pca_ind is set to 1

%% CCA with cov age
cov_path =  fullfile(working_path,'data','Stanford cohort','age_N219.mat');  % the path of covariates with a matrix of subject x covariates
cca_output_path = fullfile(output_path,'cca_output_pca85_withcov_age.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

cov_path =  fullfile(working_path,'data','Stanford cohort','study_N219.mat');  % the path of covariates with a matrix of subject x covariates
cca_output_path = fullfile(output_path,'cca_output_pca85_withcov_study.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

cov_path =  fullfile(working_path,'data','Stanford cohort','gender_N219.mat');  % the path of covariates with a matrix of subject x covariates
cca_output_path = fullfile(output_path,'cca_output_pca85_withcov_sex.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

%% CCA with pca diff parameters
cov_path =  [];  % the path of covariates with a matrix of subject x covariates
pca_percentage = 0.90; % set a threshold such as 0.9 or 0.8 if pca_ind is set to 1
cca_output_path = fullfile(output_path,'cca_output_pca90.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

pca_percentage = 0.80; % set a threshold such as 0.9 or 0.8 if pca_ind is set to 1
cca_output_path = fullfile(output_path,'cca_output_pca80.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

%% compared to main results
cca_output_main = importdata(fullfile(output_path,'cca_output_pca85.mat'));
cca_output_age = importdata(fullfile(output_path,'cca_output_pca85_withcov_age.mat'));
cca_output_sex = importdata(fullfile(output_path,'cca_output_pca85_withcov_sex.mat'));
cca_output_study = importdata(fullfile(output_path,'cca_output_pca85_withcov_study.mat'));
cca_output_pca90 = importdata(fullfile(output_path,'cca_output_pca90.mat'));
cca_output_pca80 = importdata(fullfile(output_path,'cca_output_pca80.mat'));

[r_U ~]=corr(cca_output_main.U(:,1),[cca_output_age.U(:,1),cca_output_sex.U(:,1),cca_output_study.U(:,1), cca_output_pca90.U(:,1),cca_output_pca80.U(:,1)])
[r_V ~]=corr(cca_output_main.V(:,1),[cca_output_age.V(:,1),cca_output_sex.V(:,1),cca_output_study.V(:,1),cca_output_pca90.V(:,1),cca_output_pca80.V(:,1)])

