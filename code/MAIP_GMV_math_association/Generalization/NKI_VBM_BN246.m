%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to calculate the averaged grey matter volume on
%  atlas for each participant
%
%  10-10-2021
%  Jin Liu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
clear,clc
% iMac
% oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% window
oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

%% loading subjectlist
NKI_unique_subjectlist = readtable(fullfile(oak_path,'data','behavior','NKI-RS','VBM_QC_unique_T91.dat'));

%% %% extract GMV based on atlas
file_path = fullfile(oak_path,'data','imaging','participants',filesep);
missing_list = [];
for i=1:size(NKI_unique_subjectlist,1)
    file_path1 = fullfile(file_path,strrep(NKI_unique_subjectlist.PID{i},'''',''),['visit' num2str(NKI_unique_subjectlist.Img_visit(i))],'session1','anatomical','VBM_spm12','mri',['smwp1spgr.nii']);
    if exist(file_path1)
    AllVolume{i,1} = file_path1;   
    else
        temp = [NKI_unique_subjectlist.PID(i),NKI_unique_subjectlist.Img_visit(i)] 
        missing_list = [missing_list; temp];
        ['data not found for ' file_path1]
    end
end
output_path = fullfile(oak_path,'results','smri','vbm','NKI','GMV_BN246_N91',filesep);
mkdir(output_path);
OutputName = fullfile(output_path,'GMV_BN246_N91')
ROIDef{1} = fullfile(oak_path,'data','imaging','roi','Reslice_BN_Atlas_246_2mm.nii');
MaskData = fullfile(oak_path,'data','imaging','roi','Reslice_BN_Atlas_246_2mm.nii');
IsMultipleLabel = 1;
[ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, IsMultipleLabel);

cd(output_path)

