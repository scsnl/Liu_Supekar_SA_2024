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
output_path = fullfile(oak_path,'results','smri','vbm','Stanford_cohort','wjiii','CCA110');
mkdir(output_path)

% loading behavioral data from NKI
Beh_N127 = readtable(fullfile(oak_path,'data','subjectlist','wjiii','T1_wjiii_TDlist_N218.xlsx'),'Sheet','N110');

% stanford 219 cohort
Stanford_N219 = readtable(fullfile(oak_path,'data','behavior','Stanford_cohort','subjectlist_behavior_all.xlsx'),'Sheet','Stanford_cohort_N219');

[C,IA,IB] = intersect(Stanford_N219.PID,Beh_N127.record_id)
length(find(Stanford_N219.visit(IA)==Beh_N127.Visit(IB)))

% mathfluency_sub127 = Beh_N127.MathFluency_Standard_;
applied_sub127 = Beh_N127.AppliedReason_Standard_;
% calculate_sub127 = Beh_N127.Calculation_Standard_;


%% CCA-based prediction
% loading the gray matter volume data from NKI
Stanford_GMVsub127_BN246=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','wjiii','GMV_BN246_N1110','ROISignals_GMV_BN246_N110.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','cca_output_pca85.mat'));
numopt_mathres_N219 = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','numopt_mathres_N219.mat'));

% number of subject
N = length(applied_sub127);

% scale data using PCA from stanford
X = Stanford_GMVsub127_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% define Y by filling the missing math reasoning data
Y = [applied_sub127 zeros(110,1)];

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


