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
% loading behavioral data from NKI
Beh_N91 = readtable(fullfile(oak_path,'data','behavior','NKI-RS','VBM_QC_unique_T91.dat'));

PID_sub91=Beh_N91.PID;
for i=1:length(PID_sub91)
 new_PID_sub91{i,1}=strip(PID_sub91{i},'''');
end

read_T = readtable(fullfile(oak_path,'data','behavior','NKI-RS','wiat_age_image_TD_T.dat'));

[C,IA,IB] = intersect(new_PID_sub91,read_T.AnonymizedID);
IQ_sub91 = read_T.FSIQ(IB);
VIQ_sub91 = read_T.VIQ(IB);
PIQ_sub91 = read_T.PIQ(IB);


%% CCA-based prediction
% loading the gray matter volume data from NKI
NKI_GMVsub91_BN246=importdata(fullfile(oak_path,'results','smri','vbm','NKI','GMV_BN246_N91','ROISignals_GMV_BN246_N91.mat'));

% loading the cca mode from Stanford
CCA_output=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','cca_output_pca85.mat'));
numopt_mathres_N219 = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','numopt_mathres_N219.mat'));

% number of subject
N = length(IQ_sub91);

% scale data using PCA from stanford
X = NKI_GMVsub91_BN246/CCA_output.PCA_coeff(:,1:size(CCA_output.A,1))';

% loading efficient
A = -CCA_output.A;
B = -CCA_output.B;

% define Y by filling the missing math reasoning data
Y = [IQ_sub91 zeros(N,1)];

% multiply centered X/Y and coefficient
U = (X - repmat(mean(X),N,1))*A
V = (Y - repmat(mean(Y),N,1))*B

% scale U to get the predict Y
% pred_Y = (U(:,1).*std(numopt_mathres_N219(:,1)) + mean(numopt_mathres_N219(:,1)))

% corr between predict Y and real Y
% [r1 p]=corr(pred_Y(:,1),Y(:,1))
[r p]=corr(U(:,1),Y(:,1))
% [r p]=corr(U(:,1),V(:,1))

Y = [PIQ_sub91 zeros(N,1)];
[r p]=corr(U(:,1),Y(:,1))

Y = [VIQ_sub91 zeros(N,1)];
[r p]=corr(U(:,1),Y(:,1))

%% SRS
SRS_list = readtable(fullfile(oak_path,'data','behavior','NKI-RS','8100_SRS_-_Parent_Report_20191009.csv'));

for i=1:size(Beh_N91,1)
 Beh_N91.PID{i}=strip(Beh_N91.PID{i},'''');
end

for i=1:size(Beh_N91,1)
 Beh_N91.SubStudyLabel{i}=strip(Beh_N91.SubStudyLabel{i},'''');
end

for i=1:size(Beh_N91,1)
    if find(strcmp(SRS_list.AnonymizedID,Beh_N91.PID{i}) & strcmp(SRS_list.SubStudyLabel,Beh_N91.SubStudyLabel{i}));
        SRS_N91{i,1} = SRS_list.SRSTotalScore_TScore{find(strcmp(SRS_list.AnonymizedID,Beh_N91.PID{i}) & strcmp(SRS_list.SubStudyLabel,Beh_N91.SubStudyLabel{i}))};
    end
end

NKI_SRS_N91=[Beh_N91(:,1:4) SRS_N91];
NKI_SRS_N91.Properties.VariableNames(5)= {'SRSTotalScore_TScore'};
writetable(NKI_SRS_N91,fullfile(oak_path,'data','behavior','NKI-RS','T_SRS_N91.dat'));

% double check the 0 and found that the 'A00040815' may have error to
% generate the total T score
% then delete this participant here 
% inde = 12
SRS_Tscore = str2double(NKI_SRS_N91.SRSTotalScore_TScore);
Y = [SRS_Tscore zeros(N,1)];Y(12,:)=[]
U1=U(:,1);U1(12,:)=[];

[r p]=corr(U1(:,1),Y(:,1))

%% WM
WM_list = readtable(fullfile(oak_path,'data','behavior','NKI-RS','8100_Digit_Span_20191009.csv'));

for i=1:size(Beh_N91,1)
 Beh_N91.PID{i}=strip(Beh_N91.PID{i},'''');
end

for i=1:size(Beh_N91,1)
 Beh_N91.SubStudyLabel{i}=strip(Beh_N91.SubStudyLabel{i},'''');
end

for i=1:size(Beh_N91,1)
    if find(strcmp(WM_list.AnonymizedID,Beh_N91.PID{i}) & strcmp(WM_list.SubStudyLabel,Beh_N91.SubStudyLabel{i}));
        Backward_span_N91{i,1} = WM_list.DS_BKW_TC{find(strcmp(WM_list.AnonymizedID,Beh_N91.PID{i}) & strcmp(WM_list.SubStudyLabel,Beh_N91.SubStudyLabel{i}))};
    end
end

NKI_WM_N91=[Beh_N91(:,1:4) Backward_span_N91];
NKI_WM_N91.Properties.VariableNames(5)= {'Backward_span'};
writetable(NKI_WM_N91,fullfile(oak_path,'data','behavior','NKI-RS','T_Backward_span_N91.dat'));

% double check the 0 and found that the 'A00040815' may have error to
% generate the total T score
% then delete this participant here 
% inde = 12
WM_score = str2double(NKI_WM_N91.Backward_span);
valid_ind = find(~isnan(WM_score))
Y = [WM_score(valid_ind)];

U1 = U(valid_ind,1)
[r p]=corr(U1(:,1),Y(:,1))
