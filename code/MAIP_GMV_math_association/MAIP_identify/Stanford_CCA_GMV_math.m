clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform CCA for brain-behavior analysis
%
%  Jin
%  11/3/2021  
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

%% CCA
input_path = fullfile(working_path,'data','Stanford cohort',filesep); % 
img_path = fullfile(input_path,'Stanford_GMV_BN246_N219.mat'); % the path of imaging features with a matrix of subject x features
component_num = 1; % only check the first CCA component
cov_path = [];  % the path of covariates with a matrix of subject x covariates
pca_ind = 1; % PCA will be used before CCA when pca_ind is set to 1
perm_time = 1000; % permuation test
beh_path = fullfile(working_path,'data','Stanford cohort','numopt_mathres_N219.mat'); % the path of behavior features with a matrix of subject x features
pca_percentage = 0.85; % set a threshold from 0.9 to 0.8 if pca_ind is set to 1
perm_pfwer = 1; % set to 1 if used pfwer permutation approach
% cca_output_path = fullfile(output_path,'cca_output_pca85_component2.mat');
cca_output_path = fullfile(output_path,'cca_output_pca85.mat');
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

%% CCA figures
% plot brain and behavioral canonical variate correlation
output_tiff=fullfile(output_path,'CCA_beain_beh_r_pca85_1.tif');
cca_mode_scatter(cca_output.U(:,component_num),cca_output.V(:,component_num),'Brain canonical variate','Behavioral canonical variate',output_tiff)

% plot behavioral weight
output_pic = fullfile(output_path,'CCA_beh_loading_pca85_1.tif');
cca_mode_beh_weight(weight_beh_r,{'Numerical operations','Math reasoning'},output_pic)

% generate nii for brain weight 
OutName =fullfile(output_path,'CCA_brain_loading_GMVbn246_pca85_1.nii');
Mask_path = fullfile(working_path,'data','BN_Atlas_246.nii');
BN2brainmap(weight_img_r',Mask_path,OutName);

% plot brain weight
cfg_file = 'CCA_loading_Cfg_blue_red.mat';
output_figure = fullfile(output_path,'CCA_brain_loading_GMVbn246_pca85_1.jpg');
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',OutName,cfg_file,output_figure);

% check output
cd(output_path)
