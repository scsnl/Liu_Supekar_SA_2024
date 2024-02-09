clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to estimate effect of PCA threshold on CCA
%
%  Jin
%  11/3/2021  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% iMac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');
% addpath(genpath(fullfile(filesep,'Volumes','menon','toolboxes','spm12')));
% addpath(genpath(fullfile(filesep,'Volumes','menon','toolboxes','BrainNetViewer_20191031')));

% Windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');
% addpath(genpath(fullfile('Y:','toolboxes','spm12')));
% addpath(genpath(fullfile('Y:','toolboxes','BrainNetViewer_20191031')));

% path for code
% addpath(genpath(fullfile(filesep,'Volumes','menon','projects','jinliu5','toolbox','customcolormap')))
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','GMV_math_association_analysis','PermCCA-master')))
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','GMV_math_association_analysis')))

%% loading data
output_path = fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219',filesep);
mkdir(output_path)

beh_input_path=fullfile(oak_path,'data','behavior','Stanford_cohort','subjectlist_behavior_all.xlsx');
N='N219';
beh_file_generation(beh_input_path,N,output_path)
numopt_mathres_N219 = importdata(fullfile(output_path,'numopt_mathres_N219.mat'));

%% CCA
input_path = fullfile(oak_path,'results','smri','vbm','Stanford_cohort','GMV_BN246_N219',filesep);
img_path = fullfile(input_path,'ROISignals_GMV_BN246_N219.mat'); % the path of imaging features with a matrix of subject x features
component_num = 1; % only check the first CCA component
cov_path = [];  % the path of covariates with a matrix of subject x covariates
pca_ind = 1; % PCA will be used before CCA when pca_ind is set to 1
perm_time = 1000; % permuation test 
beh_path = fullfile(output_path,'numopt_mathres_N219.mat'); % the path of behavior features with a matrix of subject x features
perm_pfwer = 1; % set to 1 if used pfwer permutation approach

x = 0.8:0.01:0.9;
for i = 1:11
pca_percentage = x(i); % set a threshold from 0.9 to 0.8 if pca_ind is set to 1
cca_output_path = fullfile(output_path,['cca_output_pca' num2str(x(i)*100) '.mat']);
[component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)
all_r(i) = component_r
end
