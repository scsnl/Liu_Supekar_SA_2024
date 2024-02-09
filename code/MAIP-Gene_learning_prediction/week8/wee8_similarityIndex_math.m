clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to predict the learning outcomes of week 8
%  intervetion
%
%  Jin
%  6/27/2023 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% iMac
% box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');

% Windows
box_path = fullfile('C:','Users','jinliu5','Box','backups','2021_Longt_math_gene');

% path for code
addpath(genpath(fullfile(box_path,'scripts','smri','vbm','final','figures_code')))

%% week8
GMV_index_w8 = importdata(fullfile(box_path,'results','smri','vbm','week8','GMV_index.mat'))
genetic_index_w8 = importdata(fullfile(box_path,'results','smri','vbm','week8','genetic_index.mat'))
behavior_w8 = importdata(fullfile(box_path,'data','behavior','week8','week8_N25.mat'))

% numops - CCA similarity
nonNAN_ind=find(~isnan(behavior_w8.numopsStd))
[r p]=corr(GMV_index_w8(nonNAN_ind),behavior_w8.numopsStd(nonNAN_ind))

% mathres - CCA similarity
nonNAN_ind=find(~isnan(behavior_w8.mathresStd))
[r p]=corr(GMV_index_w8(nonNAN_ind),behavior_w8.mathresStd(nonNAN_ind))

% numops - gene expression similarity
nonNAN_ind=find(~isnan(behavior_w8.numopsStd))
[r p]=corr(genetic_index_w8(nonNAN_ind,:),behavior_w8.numopsStd(nonNAN_ind))

% mathres - gene expression similarity
nonNAN_ind=find(~isnan(behavior_w8.mathresStd))
[r p]=corr(genetic_index_w8(nonNAN_ind,:),behavior_w8.mathresStd(nonNAN_ind))

% numops - gene expression similarity + CCA similarity
y = behavior_w8.numopsStd(nonNAN_ind)
[bb,dev,stats]=glmfit([genetic_index_w8(nonNAN_ind,:) GMV_index_w8(nonNAN_ind)],y)
y_predict = genetic_index_w8(nonNAN_ind,1)*bb(2) + genetic_index_w8(nonNAN_ind,2)*bb(3) + genetic_index_w8(nonNAN_ind,3)*bb(4) + GMV_index(nonNAN_ind)*bb(5) + bb(1)
[r p]=corr(y,y_predict)

% mathres - gene expression similarity + CCA similarity
y = behavior_w8.mathresStd(nonNAN_ind)
[bb,dev,stats]=glmfit([genetic_index_w8(nonNAN_ind,:) GMV_index_w8(nonNAN_ind)],y)
y_predict = genetic_index_w8(nonNAN_ind,1)*bb(2) + genetic_index_w8(nonNAN_ind,2)*bb(3) + genetic_index_w8(nonNAN_ind,3)*bb(4) + GMV_index(nonNAN_ind)*bb(5) + bb(1)
[r p]=corr(y,y_predict)

% mathres - gene expression similarity 
[bb,dev,stats]=glmfit([genetic_index_w8(nonNAN_ind,:)],y)
y_predict = genetic_index_w8(nonNAN_ind,1)*bb(2) + genetic_index_w8(nonNAN_ind,2)*bb(3) + genetic_index_w8(nonNAN_ind,3)*bb(4) + bb(1)
[r p]=corr(y,y_predict)

% learning gains - CCA similarity



