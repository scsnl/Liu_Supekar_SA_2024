clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to perform CCA for brain-behavior analysis
%
%  Jin
%  7/25/2022  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
%iMac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% Windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

% path for code
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','GMV_math_association_analysis','PermCCA-master')))
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','GMV_math_association_analysis')))
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','figures_code')))
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','stats_code')))

%% loading data
% setting output path
output_path = fullfile(oak_path,'results','smri','vbm','NKI','CCA91');
mkdir(output_path)

% loading behavioral data from NKI
Beh_N91 = readtable(fullfile(oak_path,'data','behavior','NKI-RS','VBM_QC_unique_T91.dat'));

PID_sub91=Beh_N91.PID;
for i=1:length(PID_sub91)
 new_PID_sub155{i,1}=strip(PID_sub91{i},'''');
end

read_T = readtable(fullfile(oak_path,'data','behavior','NKI-RS','wiat_age_image_TD_T.dat'));

[C,IA,IB] = intersect(new_PID_sub155,read_T.AnonymizedID);
wordread_sub91 = read_T.WordReading_Sta(IB);


%% CCA-based prediction
% loading the gray matter volume data from NKI
NKI_GMVsub91_BN246=importdata(fullfile(oak_path,'results','smri','vbm','NKI','GMV_BN246_N91','ROISignals_GMV_BN246_N91.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','cca_output_pca85.mat'));
numopt_mathres_N219 = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','numopt_mathres_N219.mat'));

% number of subject
N = length(wordread_sub91);

% scale data using PCA from stanford
X = NKI_GMVsub91_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% define Y by filling the missing math reasoning data
Y = [wordread_sub91 zeros(N,1)];

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A
V = (Y - repmat(mean(Y),N,1))*B

% scale U to get the predict Y
pred_Y = (U(:,1).*std(numopt_mathres_N219(:,1)) + mean(numopt_mathres_N219(:,1)))

% corr between predict Y and real Y
[r1 p]=corr(pred_Y(:,1),Y(:,1))
% [r p]=corr(U(:,1),V(:,1))

%% compared with numopts
numopt_sub91=Beh_N91.Numopts_Sta;
[r2 p]=corr(pred_Y(:,1),numopt_sub91)
% p = compare_correlation_coefficients(r1,r2,91,91)

% y = [pred_Y(:,1);pred_Y(:,1)];
% type_score = [ones(91,1);ones(91,1)*2];
% all_score = [Y(:,1);numopt_sub91];
% x = [all_score, type_score, all_score.*type_score];
% [bb,dev,stats] = glmfit(x,y)

[r3 p]=corr(numopt_sub91,wordread_sub91)

