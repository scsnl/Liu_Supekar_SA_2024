clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Dr Petra Vértes and is taken from Whitaker and Vértes, PNAS 2016 
% please cite that paper if you this code in your own work.
%
% Jin Liu modified to fit the data format in this study
% 2022-8-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% loading data
% loading math-related brain map
output_path = fullfile(working_path,'outputs','Stanford',filesep);

CCA_output=importdata(fullfile(output_path,'cca_output_pca85.mat'));
math_related_brainmap=CCA_output.weight_img_r;

% loading prerpocessed gene expression data
Table_gene=importdata(fullfile(working_path,'data','AllenBrainAtlas','BN246_geneexpression.csv'));
genes=Table_gene.textdata(2:end);

%  loading math gene 
math_gene=importdata(fullfile(working_path,'code','MAIP_GeneExpression','GeneEnrichment','math_gene_list.xlsx'));

[C,IA,IB] =intersect(genes,[math_gene(2:end,1);math_gene(2:end,2)]);
genes_math_overlap = C;

% the first column is index
gene_regional_expression=Table_gene.data(:,2:end);

% removed the region with NaN value of gene expression
num_region=find(~isnan(gene_regional_expression(:,1))); % only included region covered by donor samples N = 235

%% PLS
% setting the output folder
output_path = fullfile(working_path,'outputs','GeneExpression','mathgene_only')
mkdir(output_path)

X=gene_regional_expression(num_region,IA); % Predictors
Y=math_related_brainmap(num_region)'; % Response variable

% z-score:
clear CCA_output gene_regional_expression Table_gene
X=zscore(X);
Y=zscore(Y);

%perform full PLS and plot variance in Y explained by top 15 components
%typically top 2 or 3 components will explain a large part of the variance
%(hopefully!)
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);

% plot explained ratio
dim=15;
output_name = fullfile(output_path,'explainedratio.tiff');
plot_explained_ratio(PCTVAR,dim,output_name)

% plot three scatter
dim=3; % top 3 components above 10% variance
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim); % no need to do this but it keeps outputs tidy

% plot corr pls1 and MAIP
[r p]=corr(XS(:,1),math_related_brainmap(num_region)');
output_name = fullfile(output_path,'corr_pls1_brainloading.tiff');
beh_corr_scatter(XS(:,1),math_related_brainmap(num_region)','PLS1 score',['r = ' num2str(r,2)],['p < 0.001'],'Brain measure weights',[195 195 195]./255,[55 55 55]./255,[],output_name)
% [195 195 195]

% plot corr pls2 and MAIP
[r p]=corr(XS(:,2),math_related_brainmap(num_region)');
output_name = fullfile(output_path,'corr_pls2_brainloading.tiff');
beh_corr_scatter(XS(:,2),math_related_brainmap(num_region)','PLS2 score',['r = ' num2str(r,2)],['p < 0.001'],'Brain measure weights',[195 195 195]./255,[55 55 55]./255,[],output_name)

% plot corr pls3 and MAIP
[r p]=corr(XS(:,3),math_related_brainmap(num_region)');
output_name = fullfile(output_path,'corr_pls3_brainloading.tiff');
beh_corr_scatter(XS(:,3),math_related_brainmap(num_region)','PLS3 score',['r = ' num2str(r,2)],['p < 0.001'],'Brain measure weights',[195 195 195]./255,[55 55 55]./255,[],output_name)


%% plot pls brain loading map
% PLS1
%  the brain loading
brainmeasure1(num_region,:) = XS(:,1);
% making the nii for the brain loading
output_name = fullfile(output_path,'PLS1_brainloading_mathgeneonly.nii');
Mask_path = fullfile(working_path,'data','BN_Atlas_246.nii');
BN2brainmap(brainmeasure1,Mask_path,output_name)
% making brain figure for the brain loading
cfg_file = 'PLS_loading_Cfg.mat';
output_figure = fullfile(output_path,'PLS1_brainloading_mathgeneonly.jpg');
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',output_name,cfg_file,output_figure);
save(fullfile(output_path,'PLS1_brainloading_mathgeneonly.txt'),'brainmeasure1','-ascii','-double');


% PLS2
%  the brain loading
brainmeasure2(num_region,:) = XS(:,2);
% making the nii for the brain loading
output_name = fullfile(output_path,'PLS2_brainloading.nii');
Mask_path = fullfile(working_path,'data','BN_Atlas_246.nii');
BN2brainmap(brainmeasure2,Mask_path,output_name)
% making brain figure for the brain loading
cfg_file = 'PLS_loading_Cfg.mat';
output_figure = fullfile(output_path,'PLS2_brainloading.jpg');
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',output_name,cfg_file,output_figure);
save(fullfile(output_path,'PLS2_brainloading.txt'),'brainmeasure2','-ascii','-double');

% PLS3
%  the brain loading
brainmeasure3(num_region,:) = XS(:,3);
% making the nii for the brain loading
output_name = fullfile(output_path,'PLS3_brainloading.nii');
Mask_path = fullfile(working_path,'data','BN_Atlas_246.nii');
BN2brainmap(brainmeasure3,Mask_path,output_name)
% making brain figure for the brain loading
cfg_file = 'PLS_loading_Cfg.mat';
output_figure = fullfile(output_path,'PLS3_brainloading.jpg');
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv',output_name,cfg_file,output_figure);
save(fullfile(output_path,'PLS3_brainloading.txt'),'brainmeasure3','-ascii','-double');

%% permutation testing to assess significance of PLS result as a function of
% the number of components (dim) included:
perm_brainmap = importdata(fullfile(working_path,'code','MAIP_GeneExpression','AHBA_GeneExpression','BrainSMASH','surrogatebrainmap.txt'));

% setting permutation times
num_permutation=1000;

for dim=1:3
    
    % real PLS
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
    
    %  real PLS accumulate R2
    temp=cumsum(100*PCTVAR(2,1:dim));
    Rsquared = temp(dim);
    R(dim)=Rsquared;
    
    % permutation PLS
    for j=1:num_permutation
        
        % only the regions covered by gene expression data were considered
        Yp=perm_brainmap(j,num_region)';
        
        % zscore
        Yp=zscore(Yp);
        
        % permutation PLS
        [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
        
        % permutation PLS  accumulated R2
        temp=cumsum(100*PCTVAR(2,1:dim));
        perm_Rsq_all(j,dim) = temp(dim);
        
    end
    
    % p value for each PLS component
     culmulate_p(dim)=length(find(perm_Rsq_all(:,dim)>=Rsquared))/num_permutation;
end
% save data
output_name = fullfile(output_path,'BrainSMASH_cumulateP.mat');
save(output_name,'culmulate_p','R','perm_Rsq_all');

