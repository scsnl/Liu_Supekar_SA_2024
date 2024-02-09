%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to calculate the averaged grey matter volume on
%  brainnectome 246 atlas for each participant
%
%  6-12-2022
%  Jin Liu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
clear,clc
% iMac
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
% addpath(genpath(fullfile(filesep,'Volumes','menon','projects','jinliu5','toolbox','DPABI_V6.1_220101','DPARSF','Subfunctions')));
% addpath(genpath(fullfile(filesep,'Volumes','menon','projects','jinliu5','toolbox','DPABI_V6.1_220101','Subfunctions')));
% addpath(genpath(fullfile(filesep,'Volumes','menon','toolboxes','spm12')));

% windows
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% addpath(genpath(fullfile('Y:','projects','jinliu5','toolbox','DPABI_V6.1_220101','DPARSF','Subfunctions')));
% addpath(genpath(fullfile('Y:','projects','jinliu5','toolbox','DPABI_V6.1_220101','Subfunctions')));
% addpath(genpath(fullfile('Y:','toolboxes','spm12')));

%% loading data
M = importdata(fullfile(oak_path,'data','behavior','week8','week8_N25.mat'));
M = M(:,1:4); % PID visit session

%% extract GMV based on atlas
file_path = fullfile(oak_path,'data','imaging','participants',filesep);
missing_list = [];
for i=1:size(M,1)
    file_path1 = fullfile(file_path,num2str(M.PID(i),'%04d'),['visit' num2str(M.visit(i))],['session' num2str(M.session(i))],'anatomical','VBM_spm12','mri',['smwp1' M.file{i} '.nii']);
    if exist(file_path1)   
    AllVolume{i,1} = file_path1;
    else
        temp = [M.PID(i),M.visit(i),M.session(i),M.file(i)] 
        missing_list = ([missing_list; temp]);
        ['data not found for ' file_path1]
    end
end

output_path = fullfile(oak_path,'results','smri','vbm','week8','GMV_BN246_N25',filesep);
mkdir(output_path);
OutputName = fullfile(output_path,'GMV_BN246_N25')
ROIDef{1} = fullfile(oak_path,'data','imaging','roi','Reslice_BN_Atlas_246_2mm.nii');
MaskData = fullfile(oak_path,'data','imaging','roi','Reslice_BN_Atlas_246_2mm.nii');
IsMultipleLabel = 1;
[ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, IsMultipleLabel);

cd(output_path)

