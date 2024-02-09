clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to predict the learning outcomes of week 4
%  intervetion
%
%  Jin
%  8/17/2022  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% loading data
output_path = fullfile(working_path,'outputs','w4');
mkdir(output_path)

MAIP=importdata(fullfile(working_path,'outputs','Stanford','CCA_math_brainmap_N219.txt'));

PLS1_brainloading = importdata(fullfile(working_path,'outputs','GeneExpression','PLS1_brainloading.txt'));
PLS2_brainloading = importdata(fullfile(working_path,'outputs','GeneExpression','PLS2_brainloading.txt'));
PLS3_brainloading = importdata(fullfile(working_path,'outputs','GeneExpression','PLS3_brainloading.txt'));

% loading behavior data
pre_RT = importdata(fullfile(working_path,'data','week4_intervention','pre_nde_RT_N61.mat'));
post_RT = importdata(fullfile(working_path,'data','week4_intervention','post_nde_RT_N61.mat'));
RTdiff = importdata(fullfile(working_path,'data','week4_intervention','diffRT_N61.mat'));

% visualization of behavioral performances
[h,p,ci,stats]=ttest2(pre_RT,post_RT)

output_name = fullfile(output_path,'learninggains.tiff');
line_pre_post(pre_RT,post_RT,[233 233 233]./255,[145 145 145]./255,[-400 600],[-400:200:600],[-400:200:600], {'Pre','Post'}, ' ','NDE of RT (ms)', output_name)

output_name = fullfile(output_path,'learninggains_distribution.tiff');
histogram_tutoring(RTdiff,[200 200 200], ' ', ' ', ' ',output_name)

% loading brain data
img_path = fullfile(working_path,'data','week4_intervention','w4_GMV_BN246_N61.mat'); % the path of imaging features with a matrix of subject x features
Individual_gmv = importdata(img_path);

% calculating TSI
pls_T=[PLS1_brainloading,PLS2_brainloading,PLS3_brainloading];
valid_value_ind=find(pls_T(:,1)~=0);
genetic_index = corr(Individual_gmv(:,valid_value_ind)',pls_T(valid_value_ind,:));

%% correlation between plasticity index and learning outcomes
x = [genetic_index];
y = RTdiff;

% TSI prediction
% TSI and PLS1 only
[r p]=corr(x,y)
% plot(x,y,'o')
output_name = fullfile(output_path,'corr_genePLS1_learninggains.tiff');
beh_corr_scatter(x(:,1),y,'Transcriptome similarity index (PLS 1)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'NDE changes (ms)',[223,219,65]./255,[93 91 0]./255,' ',output_name)

% TSI and all PLS
[bb,dev,stats]=glmfit(x,y)
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)
output_name = fullfile(output_path,'corr_genePLSall_learninggains.tiff');
beh_corr_scatter(y_predict,y,'Predicted NDE changes (ms)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'NDE changes (ms)',[223,219,65]./255,[93 91 0]./255,' ',output_name)

%% validation and control analysis
% age effect
age = importdata(fullfile(working_path,'data','week4_intervention','age_N61.mat'));
[bb,dev,stats]=glmfit([x age],y);
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

% sex effect
sex = importdata(fullfile(working_path,'data','week4_intervention','sex_N61.mat'));
[bb,dev,stats]=glmfit([x sex],y);
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

mean_gene_expression = importdata(fullfile(working_path,'outputs','GeneExpression','mean_gene_expression_for_control.txt'));
% control analysis with mean gene expression
valid_value_ind=find(~isnan(mean_gene_expression(:,1)));
genetic_random_index = corr(Individual_gmv(:,valid_value_ind)',mean_gene_expression(valid_value_ind,:));
[r p]=corr(genetic_random_index,y)

% GMV effect
GMV_index = corr(Individual_gmv',MAIP');
[r p]=corr(GMV_index,y)
% structure only
% [bb,dev,stats]=glmfit([GMV_index],y)
% y_predict =  GMV_index*bb(2) + bb(1)
% [r p]=corr(y_predict,y)

% math at pre
fluency = importdata(fullfile(working_path,'data','week4_intervention','pre_fluency_N61.mat'));
[r p]=corr(fluency,y)




