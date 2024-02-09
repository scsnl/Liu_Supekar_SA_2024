clear,clc
%% QC for VBM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is used to copy pdf from VBM
%
% Jin Liu
% 5/25/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% imac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

%% copy pdf
mkdir(fullfile(oak_path,'results','smri','vbm','NKI','pdf_vbm'))
subjectlist = readtable(fullfile(oak_path,'data','behavior','NKI-RS','wiat_age_image_TD_T.dat'));

missing_list =[];
for i=1:size(subjectlist)
    
    input_path = fullfile(oak_path,'data','imaging','participants',subjectlist.AnonymizedID{i},['visit' num2str(subjectlist.Img_visit(i))],'session1','anatomical','VBM_spm12','report','catreport_spgr.pdf');

    if exist(input_path)
        output_path = fullfile(oak_path,'results','smri','vbm','NKI','pdf_vbm',[subjectlist.AnonymizedID{i},'_visit' num2str(subjectlist.Img_visit(i)), '_catreport_spgr.pdf'])
        copyfile(input_path,output_path)
    else
        temp = subjectlist.AnonymizedID{i};
        missing_list = [missing_list temp];
    end
    
end

% checked. A00074447 have poor quality of T1 spgr.nii  
