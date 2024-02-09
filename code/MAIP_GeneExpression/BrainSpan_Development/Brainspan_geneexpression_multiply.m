clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was developed to perform the expression plot of brainspan for
% genes of interest
%
% Jin Liu
% 2022-8-16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% setting path
% imac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');
oak_toolbox_path = fullfile(filesep,'Volumes','menon','projects','jinliu5');

% windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');
% oak_toolbox_path = fullfile('Y:','projects','jinliu5');

% addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','figures_code')))
% addpath(genpath(fullfile(oak_toolbox_path,'toolbox','customcolormap')))

%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% all genes
output_path = fullfile(working_path,'outputs','GeneExpression')

[PLSweighted_genelist] = readtable(fullfile(working_path,'outputs','GeneExpression','PLS1_geneWeights_descend.csv'));
input_expression_matrix = readtable(fullfile(working_path,'data','BrainSpan','expression_matrix.csv'));
input_gene_rows = readtable(fullfile(working_path,'data','BrainSpan','rows_metadata.csv'));
input_sample_columns = readtable(fullfile(working_path,'data','BrainSpan','columns_metadata.csv'));
output_name = fullfile(output_path,'PLS1_weight_brainspan.tif');

brainspan_analysis_multiply(PLSweighted_genelist,input_expression_matrix,input_gene_rows,input_sample_columns,output_name)


%% math genes only
%  loading math gene 
math_gene=importdata(fullfile(working_path,'code','MAIP_GeneExpression','GeneEnrichment','math_gene_list.xlsx'));
[C,IA,IB] =intersect(table2array(PLSweighted_genelist(:,1)),[math_gene(2:end,1);math_gene(2:end,2)]);
mathgene_PLSweighted_genelist = PLSweighted_genelist(IA,:);
output_name = fullfile(output_path,'mathgeneonly_PLS1_weight_brainspan.tif');

brainspan_analysis_multiply(mathgene_PLSweighted_genelist,input_expression_matrix,input_gene_rows,input_sample_columns,output_name)


