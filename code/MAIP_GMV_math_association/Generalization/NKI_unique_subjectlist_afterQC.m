clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is developed to generate the subjectlist for NKI-RS 
%
% Jin Liu
% 7-24-2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting path
% imac
% box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% windows
box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

%% QC for VBM
VBM_T = readtable(fullfile(box_path,'0. ELN','1. log of results','data, scripts, and results','data','NKI','subject267_info_QC.xlsx'));
subjectlist = readtable(fullfile(oak_path,'data','behavior','NKI-RS','wiat_age_image_TD_T.dat'));

% additional diagnosisi table
diagnosis_list = readtable(fullfile(oak_path,'data','behavior','NKI-RS','8100_Diagnostic_Summary_20191009.csv'));
diagnosis_ind = diagnosis_list.Diagnosis1_CODE(2:end);
diagnosis_ind(find(~ismember(diagnosis_ind,{'V71.09','V62.82','V61.20','307.6',''})))={'1'};
diagnosis_ind(find(ismember(diagnosis_ind,{'V71.09','V62.82','V61.20','307.6',''})))={'0'};
diagnosis_ind = cellfun(@str2num,diagnosis_ind,'UniformOutput',false);
diagnosis_ind = cell2mat(diagnosis_ind);

diagnosis_list_PID=diagnosis_list.AnonymizedID(2:end);
for i=1:length(diagnosis_list_PID)
    diagnosis_list_PID(i,1)=strip(diagnosis_list_PID(i),'''');
end
diagnosis_table = table(diagnosis_list_PID,diagnosis_list.Visit(2:end),diagnosis_ind);
diagnosis_table.Properties.VariableNames={'PID','visit','diagnosis'}

for i=1:size(VBM_T,1)
    if find(strcmp(diagnosis_table.PID,strip(VBM_T.PID{i},'''')) & (strcmp(diagnosis_table.visit,{'V1'}) | strcmp(diagnosis_table.visit,{'VA'})))
     VBM_T{i,17}= diagnosis_table.diagnosis(find(strcmp(diagnosis_table.PID,strip(VBM_T.PID{i},'''')) & (strcmp(diagnosis_table.visit,{'V1'}) | strcmp(diagnosis_table.visit,{'VA'})))); 
    else
     VBM_T{i,17}= 'NaN';
    end   
end

filter_QC=VBM_T.QCVBM;
filter_QC_T1=VBM_T.QCT1RawCheck;
% VBM_QC_T=VBM_T(find(VBM_T.FSIQ>80 & filter_QC>2 & filter_QC_T1>2),:);
VBM_QC_T=VBM_T(find(VBM_T.FSIQ>80 & filter_QC>2 & filter_QC_T1>2 & ~strcmp(VBM_T.Visit,'''V4''') & ~strcmp(VBM_T.Visit,'''V5''') & VBM_T.Var17==0),:);

unique_PID = unique(VBM_QC_T.PID);
for i=1:length(unique_PID)
    
    if length(find(strcmp(VBM_QC_T.PID,unique_PID(i))))==1
        
        VBM_QC_unique_T(i,:)=VBM_QC_T(find(strcmp(VBM_QC_T.PID,unique_PID(i))),:);
    else
        temp=VBM_QC_T(find(strcmp(VBM_QC_T.PID,unique_PID(i))),:);
        % VBM_QC_unique_T(i,:)=temp(find(min(temp.Age)),:);
        VBM_QC_unique_T(i,:)=temp(find(max(temp.QCVBM)),:);
    end
end
output = fullfile(oak_path,'data','behavior','NKI-RS','VBM_QC_unique_T91.dat');
writetable(VBM_QC_unique_T,output)
