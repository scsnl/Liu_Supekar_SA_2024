clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to generate figure for the mean map of grey matter volume on
%  brainnectome 246 atlas across participants
%
%  7-22-2022
%  Jin Liu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% figure for meanmap of GMV
% output folder
output_path = fullfile(working_path,'outputs','NKI',filesep);
mkdir(output_path);

% loading ROI mask
Mask_path = fullfile(working_path,'data','BN_Atlas_246.nii');

% loading participants' ROI signals (219 participants x 246 ROIs)
GMV_BN246=importdata(fullfile(working_path,'data','NKI-RS','NKI_GMV_BN246_N91.mat'));

% Mean map nii output name
OutName =fullfile(output_path,'mean_GMV_BN246_N91.nii');

% generating mean map nii
BN2brainmap(mean(GMV_BN246),Mask_path,OutName);

% BNV setting
cfg_file = fullfile(working_path,'code','MAIP_GMV_math_association','MAIP_identify','meanmap_GMV_Cfg.mat');

% BNV figure name
output_figure = fullfile(output_path,'mean_GMV_BN246_N91.jpg');

% generating BNV image
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',OutName,cfg_file,output_figure);

