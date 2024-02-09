clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to evaluate the spatial correlation between 
%  CCA-weighted brain map and math-related meta maps
%
%  Jin
%  7/21/2022  
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
list_metamap_topic = dir(fullfile(input_path,'meta_map_association_uni_topic'));
list_metamap_topic(1:2,:)=[];

% reorder the meta map (first 5 for math, and last 2 for reading as control)
list_metamap_topic = list_metamap_topic([4,5,1,2,3,6],:); % first 2 for math, last 4 for reading
list_metamap = [list_metamap_topic];

%% calculating correlation between CCA brain loading and meta-mask
% calculating correlation between CCA brain loading and each meta-mask
for i=1:length(list_metamap)
    
    % meta map file
    V_temp = spm_vol(fullfile(list_metamap(i).folder,list_metamap(i).name));
    
    % reading meta map
    meta_mask = spm_read_vols(V_temp);
    
    % mean weight within meta mask
    meanweight_meta(i)=mean(CCA_math_brainloading(meta_mask~=0));
    
    % calculating the correlation between CCA brain loading and meta map
    [corr_coef(i,1) corr_coef(i,2)]= corr(CCA_math_brainloading(CCA_math_brainloading~=0),meta_mask(CCA_math_brainloading~=0),'type','Spearman');

end
save(fullfile(output_path,'corr_CCAloading_metanmap_topic.mat'),'corr_coef');


%% extract mean Z-score based on 246 ROIs with 20 bins label 
% setting output folder 
output_bin_path = fullfile(output_path,'meta_bins20_topic');
mkdir(output_bin_path)

% Extrating setting
for i=1:length(list_metamap)
    file_path1 = fullfile(list_metamap(i).folder,list_metamap(i).name);
    AllVolume{i,1} = file_path1;
end

% Output name
OutputName = fullfile(output_bin_path,'zscore_meta_brainloading_bins20_topic');

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
        ind_perm=randperm(numel(loading_ture_index));
        
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
save(fullfile(output_path,'corr_coef_permutation_p_topic.mat'),'corr_coef_random','p_r_metamap_MAIP');

%% comparing the correlation using permutation tests
num_perm = 1000;

for i=3:6 % first two math terms
    i
    % read math map
    V_temp = spm_vol(fullfile(list_metamap(i).folder,list_metamap(i).name));
    meta_mask = spm_read_vols(V_temp);
    
    for l = 1:2
        % read reading map
        V_temp_math1 = spm_vol(fullfile(list_metamap(l).folder,list_metamap(l).name));
        meta_mask_math1 = spm_read_vols(V_temp_math1);
        
        % permutation with the same size
        for j=1:num_perm
            j
            
            % reset null
            null_mask = zeros(121,145,121);
            null_mask_math1 = zeros(121,145,121);
            
            % find all valid index within mask
            loading_ture_index=find(CCA_math_brainloading~=0);
            
            % calculate the size of meta_mask
            N_size=length(find(meta_mask~=0));
            N_size_math=length(find(meta_mask_math1~=0));
            
            % ramdon index
            ind_perm=randperm(numel(loading_ture_index));
            ind_perm_math=randperm(numel(loading_ture_index));
            
            % ramdon mask
            null_mask(loading_ture_index(ind_perm(1:N_size)))=meta_mask(find(meta_mask~=0));
            null_mask_math1(loading_ture_index(ind_perm_math(1:N_size_math)))=meta_mask_math1(find(meta_mask_math1~=0));
            
            % calculating the correlation between CCA brain loading and meta map
            [corr_coef_random_read ~]= corr(CCA_math_brainloading(CCA_math_brainloading~=0),null_mask(CCA_math_brainloading~=0),'type','Spearman');
            [corr_coef_random_math ~]= corr(CCA_math_brainloading(CCA_math_brainloading~=0),null_mask_math1(CCA_math_brainloading~=0),'type','Spearman');
            
            corr_coef_diff(j,:)=corr_coef_random_math-corr_coef_random_read;
        end
        all_corr_coef_diff{i,l}=corr_coef_diff;
    end
end

for i=3:6 % first four math terms
    for l = 1:2   
       real_diff =  corr_coef(l)-corr_coef(i);
       ramdon_diff=all_corr_coef_diff{i,l};
       p_r_metamap_bw(i,l)=length(find(ramdon_diff>abs(real_diff)))./num_perm
    end
end

save(fullfile(output_path,'corr_bw_metanmap_topic.mat'),'all_corr_coef_diff','p_r_metamap_bw');

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
term_name = {'arithmetic','number','sentences','phonological','words','language'};

% importing mean Z-score based on 246 ROIs with 20 bins label 
ROISignals= importdata(fullfile(output_bin_path,'ROISignals_zscore_meta_brainloading_bins20_topic.mat'));

% rank the data based on correlation
[B,I] = sort(corr_coef(:,1),'descend');
ROISignals_rank=ROISignals(I,:);
term_name_rank = term_name(I);
corr_coef_rank = corr_coef(I,1);
corr_coef_random_rank = corr_coef_random(I,:);

% matrix for math
output_name = fullfile(output_path,'matrix_math.tiff')
matrix_math(ROISignals_rank(1:2,:),term_name_rank(1:2),[480,280],output_name)

% matrix for read
output_name = fullfile(output_path,'matrix_read.tiff')
matrix_math(ROISignals_rank(3:6,:),term_name_rank(3:6),[480,480],output_name)

% boxplot for math
output_name = fullfile(output_path,'boxplot_math.tiff')
boxplot_r_jin(corr_coef_random_rank(1:2,:),corr_coef_rank(1:2),term_name_rank(1:2),[350,180],output_name)

output_name = fullfile(output_path,'boxplot_read.tiff')
boxplot_r_jin(corr_coef_random_rank(3:6,:),corr_coef_rank(3:6),term_name_rank(3:6),[350,400],output_name)

% check figure
cd(output_path)