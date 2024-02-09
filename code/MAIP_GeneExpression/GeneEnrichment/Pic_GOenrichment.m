clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was generated to plot the enrichment results
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

%% pic for GO enrichment
output_path = fullfile(working_path,'outputs','GeneExpression')

data = readtable('PLS1_GOCOMPONENT_descend.xlsx')
output_name = fullfile(output_path,'PLS1minus_GOComp_nolabel.tiff');
bar_GOenrich(data,'PLS1- enriched','Cellular component',10,output_name);

data = readtable('PLS1_GOFUNCTION_descend.xlsx')
output_name = fullfile(output_path,'PLS1minus_GOFunc_nolabel.tiff');
bar_GOenrich(data,'PLS1- enriched','Molecular Function',10,output_name);

data = readtable('PLS1_GOPROCESS_descend.xlsx')
output_name = fullfile(output_path,'PLS1minus_GOProc_nolabel.tiff');
bar_GOenrich(data,'PLS1- enriched','Biological process',9,output_name);


