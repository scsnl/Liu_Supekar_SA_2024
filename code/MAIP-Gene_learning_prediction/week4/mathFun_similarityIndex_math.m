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
GMV_index_mathfun = importdata(fullfile(box_path,'results','smri','vbm','mathFUN','GMV_index.mat'))
genetic_index_mathfun= importdata(fullfile(box_path,'results','smri','vbm','mathFUN','genetic_index.mat'))
behavior_mathfun = importdata(fullfile(box_path,'data','behavior','mathFUN','mathfun_all_data.mat'))
mathFUNtutoringN62 = importdata(fullfile(box_path,'data','behavior','mathFUN','mathFUN_tutoring_N62.mat'))

[C,IA,IB] = intersect(behavior_mathfun(:,1),mathFUNtutoringN62(:,1))

% fluency - CCA similarity
[r p]=corr(GMV_index_mathfun(IB),behavior_mathfun.wjiii_math_fluency_std(IA))

% mathres - CCA similarity
[r p]=corr(GMV_index_mathfun(IB),behavior_mathfun.wjiii_applied_reasoning_std(IA))

% fluency - gene expression similarity
[r p]=corr(genetic_index_mathfun(IB,:),behavior_mathfun.wjiii_math_fluency_std(IA))

% mathres - gene expression similarity
[r p]=corr(genetic_index_mathfun(IB,:),behavior_mathfun.wjiii_applied_reasoning_std(IA))

% fluency - gene expression similarity + CCA similarity
y = behavior_mathfun.wjiii_math_fluency_std(IA)
[bb,dev,stats]=glmfit([genetic_index_mathfun(IB,:) GMV_index_mathfun(IB)],y)
y_predict = genetic_index_mathfun(IB,1)*bb(2) + genetic_index_mathfun(IB,2)*bb(3) + genetic_index_mathfun(IB,3)*bb(4) + GMV_index_mathfun(IB)*bb(5) + bb(1)
[r p]=corr(y,y_predict)

% mathres - gene expression similarity + CCA similarity
y = behavior_mathfun.wjiii_applied_reasoning_std(IA)
[bb,dev,stats]=glmfit([genetic_index_mathfun(IB,:) GMV_index_mathfun(IB)],y)
y_predict = genetic_index_mathfun(IB,1)*bb(2) + genetic_index_mathfun(IB,2)*bb(3) + genetic_index_mathfun(IB,3)*bb(4) + GMV_index_mathfun(IB)*bb(5) + bb(1)
[r p]=corr(y,y_predict)

% mathres - gene expression similarity 
[bb,dev,stats]=glmfit([genetic_index_mathfun(IB,:)],y)
y_predict = genetic_index_mathfun(IB,1)*bb(2) + genetic_index_mathfun(IB,2)*bb(3) + genetic_index_mathfun(IB,3)*bb(4) + bb(1)
[r p]=corr(y,y_predict)





