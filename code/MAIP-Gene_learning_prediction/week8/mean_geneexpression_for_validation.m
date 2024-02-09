clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was develop to calculate the mean gene expression across all
% genes as control map
%
% Jin Liu 
% 2022-8-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% getting average gene expression
% loading prerpocessed gene expression data
Table_gene=importdata(fullfile(working_path,'data','AllenBrainAtlas','BN246_geneexpression.csv'));

% the first column is index
gene_regional_expression=Table_gene.data(:,2:end);

% removed the region with NaN value of gene expression
mean_gene_regional_expression=mean(gene_regional_expression,2); % only included region covered by donor samples N = 235

% output path
output_path = fullfile(working_path,'outputs','GeneExpression',filesep);

% save output
save(fullfile(output_path,'mean_gene_expression_for_control.txt'),'mean_gene_regional_expression','-ascii','-double');


