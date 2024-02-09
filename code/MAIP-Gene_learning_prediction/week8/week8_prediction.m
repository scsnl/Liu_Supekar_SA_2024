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
% 
% Windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');
% oak_path = fullfile('C:','Users','pbs','Box','backups','2021_Longt_math_gene');

% path for code
addpath(genpath(fullfile(oak_path,'scripts','smri','vbm','final','figures_code')))

%% loading data
PLS1_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS1_brainloading.txt'));
PLS2_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS2_brainloading.txt'));
PLS3_brainloading = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','PLS3_brainloading.txt'));
mean_gene_expression = importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','PLS_geneexpression_GMVmath_N219','mean_gene_expression.txt'));

MAIP=importdata(fullfile(oak_path,'results','smri','vbm','Stanford_cohort','CCA_GMV_math_N219','CCA_math_brainmap_N219.txt'));

img_path = fullfile(oak_path,'results','smri','vbm','week8','GMV_BN246_N25','ROISignals_GMV_BN246_N25.mat'); % the path of imaging features with a matrix of subject x features
Individual_gmv = importdata(img_path);
   
%% loading behavior data
beh_meas =importdata(fullfile(oak_path,'data','behavior','week8','week8_N25.mat'));
final_subjectlist = beh_meas.PID;
RT_T=readtable(fullfile(oak_path,'data','behavior','week8','20120114 behavior all_8week.xlsx'),'Sheet','RT graphs');
RT_subjlist = cellfun(@str2num,RT_T.Subject(1:39));
RT_T(40:42,:)=[];
age = beh_meas.age;
sex = beh_meas.gender;

for i=1:length(sex)
    if strcmp(sex{i},'Male')
        sex_ind(i,1) = 1;
    elseif strcmp(sex{i},'Female')
        sex_ind(i,1) = 2;
    end
end
[C,IA,IB] = intersect(final_subjectlist,RT_subjlist);

%% loading brain data
pls_T=[PLS1_brainloading,PLS2_brainloading,PLS3_brainloading];
valid_value_ind=find(pls_T(:,1)~=0);
genetic_index = corr(Individual_gmv(:,valid_value_ind)',pls_T(valid_value_ind,:));
GMV_index = corr(Individual_gmv',MAIP');


%% control analysis with mean gene expression
valid_value_ind=find(~isnan(mean_gene_expression(:,1)));
genetic_random_index = corr(Individual_gmv(:,valid_value_ind)',mean_gene_expression(valid_value_ind,:));

% output_filename = fullfile(oak_path,'results','smri','vbm','week8','GMV_index.mat')
% save(output_filename,'GMV_index')
% 
% output_filename = fullfile(oak_path,'results','smri','vbm','week8','genetic_index.mat')
% save(output_filename,'genetic_index')

%% correlation between plasticity index and learning outcomes
% corr 
x = [genetic_index(IA,:)];
w8_GMV_BN246_N24 = Individual_gmv(IA,:);

age = age(IA);
y = RT_T.additiondiffRT(IB);
mathres = beh_meas.mathresStd(IA,:);
numops = beh_meas.numopsStd(IA,:);
sex_ind = sex_ind(IA);

RT_subjlist_new=RT_subjlist(IB);
RT_subjlist_new(find(y>mean(y)+3*std(y)));
pre_RT=RT_T.PREaddadd(IB);
pre_RT(find(y>mean(y)+3*std(y)),:)=[];
post_RT=RT_T.POSTaddadd(IB);
post_RT(find(y>mean(y)+3*std(y)),:)=[];

x(find(y>mean(y)+3*std(y)),:)=[];
w8_GMV_BN246_N24(find(y>mean(y)+3*std(y)),:)=[];
GMV_index(find(y>mean(y)+3*std(y)),:)=[];
mathres(find(y>mean(y)+3*std(y)),:)=[];
numops(find(y>mean(y)+3*std(y)),:)=[];
genetic_random_index(find(y>mean(y)+3*std(y)),:)=[];
age(find(y>mean(y)+3*std(y)),:)=[];
sex_ind(find(y>mean(y)+3*std(y)),:)=[];
y(find(y>mean(y)+3*std(y)),:)=[];
% output_name = fullfile(oak_path,'results','smri','vbm','week8','learninggains_distribution.tiff');
% histogram_tutoring(y,[200 200 200], ' ', ' ', ' ',output_name)

[h,p,ci,stats]=ttest2(pre_RT,post_RT)
% output_name = fullfile(oak_path,'results','smri','vbm','week8','learninggains.tiff');
% line_pre_post(pre_RT,post_RT,[233 233 233]./255,[145 145 145]./255,[0 7000],[0:2000:6000],[0:2000:6000], {'Pre','Post'}, ' ','RT (ms)', output_name)

[r p]=corr(x,y)
% output_name = fullfile(oak_path,'results','smri','vbm','week8','corr_genePLS1_learninggains.tiff');
% beh_corr_scatter(x(:,1),y','Transcriptome similarity index (PLS 1)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'RT changes (ms)',[131,220,228]./255,[10 137 148]./255,' ',output_name)
% 
[r p]=corr(GMV_index,y)
% [r p]=corr(x,GMV_index)
% [r p]=corr(GMV_index(find(~isnan(mathres))),mathres(find(~isnan(mathres))))
% [r p]=corr(GMV_index(find(~isnan(numops))),numops(find(~isnan(numops))))

[r p]=corr(mathres(find(~isnan(mathres))),y(find(~isnan(mathres))))
[r p]=corr(numops(find(~isnan(numops))),y(find(~isnan(numops))))


% [bb,dev,stats]=glmfit([x GMV_index],y)
% y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) +  + GMV_index*bb(5) + bb(1)
% [r p]=corr(y_predict,y)
% output_name = fullfile(oak_path,'results','smri','vbm','week8','corr_genePLSall_learninggains.tiff');
% beh_corr_scatter(y_predict,y','Predicted RT changes (ms)',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'RT changes (ms)',[131,220,228]./255,[10 137 148]./255,' ',output_name)

% TSI prediction
[bb,dev,stats]=glmfit([x],y)
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

% structure only
[bb,dev,stats]=glmfit([GMV_index],y)
y_predict =  GMV_index*bb(2) + bb(1)
[r p]=corr(y_predict,y)

% random gene expression profile prediction
[r p]=corr(genetic_random_index,y)

% age effect
[bb,dev,stats]=glmfit([x age],y)
y_predict = x(:,1)*bb(2) + x(:,2)*bb(3) + x(:,3)*bb(4) + bb(1)
[r p]=corr(y_predict,y)

[bb,dev,stats]=glmfit([x sex_ind],y)
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
% % plot(XS(:,1),Y,'o')
% [r p]=corr(XS(:,1),y)
% output_name = fullfile(oak_path,'results','smri','vbm','week8','corr_gene_learninggains.tiff');
% beh_corr_scatter(XS(:,1),y','Transcriptome similarity index',['r = ' num2str(r(1),2)],['p  = ' num2str(p(1),3)],'RT changes (ms)',[131,220,228]./255,[10 137 148]./255,' ',output_name)
% 
% num_permutation=1000;
% dim=1
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
% 
