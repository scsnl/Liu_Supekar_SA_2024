clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to predict mathFun's fluency based on CCA for 
%  brain-behavior analysis
%
%  Jin
%  12/4/2023  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% iMac
% box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','backups','2021_Longt_math_gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% Windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

%% check if structure image exist
input_file = fullfile(oak_path,'data','subjectlist','wjiii','T1_wjiii_TDlist_N218.xlsx');
T1_wjiii_all=readtable(input_file,'Sheet','N131');
j=0
missing_list=[];
subjectlist=[];
subjectlist_T1=[];
for  i=1:length(T1_wjiii_all.record_id)
    
    if exist(fullfile(filesep,'Volumes','menon','rawdata','scsnl',num2str(T1_wjiii_all.record_id(i)),['visit' num2str(T1_wjiii_all.Visit(i))],'session1','anatomical','spgr.nii.gz'))
        j=j+1;
        subjectlist(j,1)=T1_wjiii_all.record_id(i);
        subjectlist(j,2)=T1_wjiii_all.Visit(i);
        subjectlist(j,3)=1;
        subjectlist_T1{j,1}='spgr';
    elseif exist(fullfile(filesep,'Volumes','menon','rawdata','scsnl',num2str(T1_wjiii_all.record_id(i)),['visit' num2str(T1_wjiii_all.Visit(i))],'session1','anatomical','spgr_1.nii.gz'))
        j=j+1;
        subjectlist(j,1)=T1_wjiii_all.record_id(i);
        subjectlist(j,2)=T1_wjiii_all.Visit(i);
        subjectlist(j,3)=1;
        subjectlist_T1{j,1}='spgr_1';
    elseif exist(fullfile(filesep,'Volumes','menon','rawdata','scsnl',num2str(T1_wjiii_all.record_id(i)),['visit' num2str(T1_wjiii_all.Visit(i))],'session2','anatomical','spgr.nii.gz'))
        j=j+1;
        subjectlist(j,1)=T1_wjiii_all.record_id(i);
        subjectlist(j,2)=T1_wjiii_all.Visit(i);
        subjectlist(j,3)=2;
        subjectlist_T1{j,1}='spgr';
    else
        missing_list = [missing_list;num2str(T1_wjiii_all.record_id(i)) '_visit' num2str(T1_wjiii_all.Visit(i))];
    end
    
end
