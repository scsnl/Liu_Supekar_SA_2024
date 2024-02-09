clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  This script is used to predict the learning outcomes of week 8
%  intervetion
%
%  Jin
%  8/17/2022  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% iMac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');
% oak_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','backups','2021_Longt_math_gene');

% Windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');
% % oak_path = fullfile('C:','Users','pbs','Box','backups','2021_Longt_math_gene');
% box_path = fullfile('C:','Users','jinliu5','Box','backups','2021_Longt_math_gene');

% path for code
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','figures_code')));

%% loading data
PLS1_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS1_brainloading.txt'));
PLS2_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS2_brainloading.txt'));
PLS3_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS3_brainloading.txt'));
mean_gene_expression = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','mean_gene_expression.txt'));

MAIP=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','CCA_math_brainmap_N219.txt'));


img_path = fullfile(oak_path,'results','smri','vbm','mathFUN','GMV_BN246_N62','ROISignals_GMV_BN246_N62.mat'); % the path of imaging features with a matrix of subject x features
Individual_gmv = importdata(img_path);

%% loading behavior data
beh_meas =readtable(fullfile(oak_path,'data','behavior','mathFUN','mathFUN_tutoring_N62.dat'));
final_subjectlist = beh_meas.PID;
cd(fullfile(oak_path,'data','behavior','mathFUN'))

beh_T=readtable('mathFUN_compnum_data_vbm.csv');
w4_subjlist = beh_T.PID;
RT_nde_diff = beh_T.num_meanrt_nde_gain;

[C,IA,IB] = intersect(final_subjectlist,w4_subjlist);

%% loading brain data
pls_T=[PLS1_brainloading,PLS2_brainloading,PLS3_brainloading];
valid_value_ind=find(pls_T(:,1)~=0);
genetic_index = corr(Individual_gmv(:,valid_value_ind)',pls_T(valid_value_ind,:));
GMV_index = corr(Individual_gmv',MAIP');

% output_filename = fullfile(oak_path,'results','smri','vbm','mathFUN','GMV_index.mat')
% save(output_filename,'GMV_index')

% output_filename = fullfile(oak_path,'results','smri','vbm','mathFUN','genetic_index.mat')
% save(output_filename,'genetic_index')

%% control analysis with mean gene expression
valid_value_ind=find(~isnan(mean_gene_expression(:,1)));
genetic_random_index = corr(Individual_gmv(:,valid_value_ind)',mean_gene_expression(valid_value_ind,:));

%% correlation between plasticity index and learning outcomes
% RT
x = [genetic_index(IA,:)];
y = RT_nde_diff(IB);
fluency = beh_meas.pre_fluency(IA,:);

behavior_mathfun = importdata(fullfile(oak_path,'data','behavior','mathFUN','mathfun_all_data.mat'))
[C,IC,ID] = intersect(final_subjectlist,behavior_mathfun.PID);
age = behavior_mathfun.Age(ID);

w4_subjlist(find(y>mean(y)+3*std(y)))
IB_new = IB;
IB_new(find(y>mean(y)+3*std(y)),:)=[];
x(find(y>mean(y)+3*std(y)),:)=[];
GMV_index(find(y>mean(y)+3*std(y)),:)=[];
fluency(find(y>mean(y)+3*std(y)),:)=[];
age(find(y>mean(y)+3*std(y)),:)=[];
y(find(y>mean(y)+3*std(y)),:)=[];

subjectlist_N61 = readtable(fullfile(oak_path,'data','behavior','mathFUN','week4_N61_demograph.csv'));
sex = subjectlist_N61.gender_1_F_2_M_;

output_name = fullfile(oak_path,'results','smri','vbm','mathFUN','learninggains_distribution.tiff');
histogram_tutoring(y,[200 200 200], ' ', ' ', ' ',output_name)

pre_RT=beh_T.num_meanrt_nde_1(IB_new);
pre_RT(find(y>mean(y)+3*std(y)),:)=[];
post_RT=beh_T.num_meanrt_nde_2(IB_new);
post_RT(find(y>mean(y)+3*std(y)),:)=[];

[h,p,ci,stats]=ttest2(pre_RT,post_RT)
output_name = fullfile(oak_path,'results','smri','vbm','mathFUN','learninggains.tiff');
line_pre_post(pre_RT,post_RT,[233 233 233]./255,[145 145 145]./255,[-400 600],[-400:200:600],[-400:200:600], {'Pre','Post'}, ' ','NDE of RT (ms)', output_name)

[r p]=corr(x,y)
% plot(x,y,'o')
output_name = fullfile(oak_path,'results','smri','vbm','mathFUN','corr_genePLS1_learninggains.tiff');
beh_corr_scatter(x(:,1),y','Transcriptome similarity index (PLS 1)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'NDE changes (ms)',[223,219,65]./255,[93 91 0]./255,' ',output_name)

[r p]=corr(GMV_index,y)
% [r p]=corr(x,GMV_index)
[r p]=corr(fluency,y)
% [r p]=corr(fluency,GMV_index)
[r p]=corr(genetic_random_index,y)

[bb,dev,stats]=glmfit(x,y)
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

% [bb,dev,stats]=glmfit([x pre_RT],y)
% y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + pre_RT*bb(5) + bb(1)
% [r p]=corr(y_predict,y)


% [bb,dev,stats]=glmfit([x GMV_index],y)
% y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) +  + GMV_index*bb(5) + bb(1)
% [r p]=corr(y_predict,y)
% output_name = fullfile(oak_path,'results','smri','vbm','mathFUN','corr_genePLSall_learninggains.tiff');
% beh_corr_scatter(y_predict,y','Predicted NDE changes (ms)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'NDE changes (ms)',[223,219,65]./255,[93 91 0]./255,' ',output_name)

% only structure prediction
[bb,dev,stats]=glmfit([GMV_index],y)
y_predict =  GMV_index*bb(2) + bb(1)
[r p]=corr(y_predict,y)

% age effect
[bb,dev,stats]=glmfit([x age],y);
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

[bb,dev,stats]=glmfit([x sex],y);
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

[bb,dev,stats]=glmfit([x age sex],y);
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)
% %% pls
% % z-score:
% X=zscore(x);
% Y=zscore(y);
% 
% %perform full PLS and plot variance in Y explained by top 15 components
% %typically top 2 or 3 components will explain a large part of the variance
% %(hopefully!)
% [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y);
% [r p]=corr(XS(:,1),Y)
% output_name = fullfile(oak_path,'results','smri','vbm','mathFUN','corr_gene_learninggains.tiff');
% beh_corr_scatter(XS(:,1),y','Transcriptome similarity index',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'NDE changes (ms)',[223,219,65]./255,[93 91 0]./255,' ',output_name)
% 
% num_permutation=1000;
% for dim=1:3
% [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
% temp=cumsum(100*PCTVAR(2,1:dim));
% Rsquared = temp(dim);
%     for j=1:num_permutation
%         j
%          order=randperm(size(Y,1));
%          Yp=Y(order,:);
%         Yp=zscore(Yp);
%         [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Yp,dim);
%         temp=cumsum(100*PCTVAR(2,1:dim));
%         Rsq(j) = temp(dim);
%     end
% dim
% R(dim)=Rsquared;
% p(dim)=length(find(Rsq>=Rsquared))/num_permutation
% end
