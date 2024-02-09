clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to evaluate the spatial correlation between 
%  CCA-weighted brain map and math-related meta maps
%
%  Jin
%  7/21/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');
% addpath('DPABI'));

%% loading CCA results
% loading CCA results
% CCA result folder
output_path = fullfile(working_path,'outputs','Stanford',filesep);

% CCA brain loading nii file
V = spm_vol(fullfile(output_path,'CCA_brain_loading_GMVbn246_pca85_1.nii'));

% reading the CCA brain loding nii 
CCA_math_brainloading = spm_read_vols(V);

% reading other CCA results
CCA_output=importdata(fullfile(output_path,'cca_output_pca85.mat'));

%% loading meta map
input_path = fullfile(working_path,'data','Stanford Cohort',filesep);

% meta map folder
list_metamap = dir(fullfile(input_path,'meta_map_childrenonly'))
list_metamap(1:2,:)=[];

%% calculating correlation between CCA brain loading and meta-mask
% calculating correlation between CCA brain loading and each meta-mask

% meta map file
V_temp = spm_vol(fullfile(list_metamap(1).folder,list_metamap(1).name));

% reading meta map
meta_mask = spm_read_vols(V_temp);

% mean weight within meta mask
meanweight_meta(1)=mean(CCA_math_brainloading(meta_mask~=0));

% calculating the correlation between CCA brain loading and meta map
[corr_coef(1) corr_coef(2)]= corr(CCA_math_brainloading(CCA_math_brainloading~=0),meta_mask(CCA_math_brainloading~=0),'type','Spearman');

save(fullfile(output_path,'corr_CCAloading_metanmap_childrenonly.mat'),'corr_coef');


%% extract mean Z-score based on 246 ROIs with 20 bins label 
% Input data setting
output_bin_path = fullfile(output_path,'meta_bins20_childrenonly');
mkdir(output_bin_path)

% Extrating setting
file_path1 = fullfile(list_metamap(1).folder,list_metamap(1).name);
AllVolume{1,1} = file_path1;

% Output name
OutputName = fullfile(output_bin_path,'zscore_meta_brainloading_bins20');

% ROI mask
ROIDef{1} = fullfile(output_path,'CCA_brain_loading_bin20_ROI.nii');
MaskData = fullfile(working_path,'data','BN_Atlas_246.nii');
IsMultipleLabel = 1; 

% Extracting ROI signal
[ROISignals] = y_ExtractROISignal(AllVolume, ROIDef, OutputName, MaskData, IsMultipleLabel);

% Check output data
cd(output_bin_path);

%% permutation test compared to random maps
% permutation times
num_perm=1000;

for i=1:length(list_metamap)
    i
    
    % read each meta mask
    V_temp = spm_vol(fullfile(list_metamap(i).folder,list_metamap(i).name));
    meta_mask = spm_read_vols(V_temp);
    
    % permutation with the same size
    for j=1:num_perm
        j
        
        % reset null
        null_mask = zeros(121,145,121);
        
        % find all valid index within mask
        loading_ture_index=find(CCA_math_brainloading~=0);
        
        % calculate the size of meta_mask
        N_size=length(find(meta_mask~=0));
        
        % ramdon index
        ind_perm(:,1)=randperm(numel(loading_ture_index));
        
        % ramdon mask
        null_mask(loading_ture_index(ind_perm(1:N_size)))=1;
        
        % calculating the correlation between CCA brain loading and meta map
        [corr_coef_random(i,j) ~]= corr(CCA_math_brainloading(CCA_math_brainloading~=0),null_mask(CCA_math_brainloading~=0),'type','Spearman');
        
    end
end

% permutation test for p of meta map
for i=1:length(list_metamap)
    p_r_metamap_MAIP(i,1)=length(find(corr_coef_random(i,:)>corr_coef(i,1)))./num_perm
end
save(fullfile(output_path,'corr_coef_permutation_p_childrenonly.mat'),'corr_coef_random','p_r_metamap_MAIP');

%% visualization
% brain maps for each terms
for i=1:length(list_metamap)
    
    % nii file
    input_nii=fullfile(list_metamap(i).folder,list_metamap(i).name);
    
    % BNV setting
    cfg_file = 'meta_term_Cfg.mat';
    
    % BNV figure name
    output_figure = fullfile(output_path,[list_metamap(i).name '.jpg']);
    
    % generating BNV image
    BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',input_nii,cfg_file,output_figure);
    
end

% distribution matrix
% term name
term_name = {'Math children'}

% importing mean Z-score based on 246 ROIs with 20 bins label 
ROISignals= importdata(fullfile(output_bin_path,'ROISignals_zscore_meta_brainloading_bins20.mat'))

% rank the data based on correlation
[B,I] = sort(corr_coef(:,1),'descend');
ROISignals_rank=ROISignals(I,:);
term_name_rank = term_name(I);
corr_coef_rank = corr_coef(I,1);
corr_coef_random_rank = corr_coef_random(I,:);

% matrix for math
output_name = fullfile(output_path,'matrix_math_children.tiff')
matrix_math(ROISignals_rank(1,:),term_name_rank(1),[480,150],output_name)

% boxplot for math
output_name = fullfile(output_path,'boxplot_math_chilren.tiff')
boxplot_r_jin(corr_coef_random_rank(1,:),corr_coef_rank(1),term_name_rank(1),[350,150],output_name)

% check figure
cd(output_path)

