function [component_r,permutation_P,weight_img_r,weight_img_p,weight_beh_r,weight_beh_p,cca_output]=CCA_img_beh_jin(img_path,beh_path,component_num,cov_path,pca_ind,pca_percentage,perm_time,perm_pfwer,cca_output_path)

%% loading data
c = fix(clock);
disp('==================================================================');
fprintf('CCA analysis started at %d/%02d/%02d %02d:%02d:%02d \n',c');

% load imaging data
img_meas=importdata(img_path);
% load behavior data
beh_meas=importdata(beh_path);
% load covariates 
if isempty(cov_path)
    Cov=[];
else
    Cov=importdata(cov_path);
end

%% data preparation
% fill the missing data in behavior data by group mean value
for i=1:size(beh_meas,2)
    [beh_meas(find(isnan(beh_meas(:,i))),i)] = mean(beh_meas(~isnan(beh_meas(:,i)),i));
end

% regress of covariate of no interest
if ~isempty(Cov)
    clear R
    for j=1:size(img_meas,2)
        [B,BINT,R(:,j)] = regress(img_meas(:,j),[ones(size(img_meas,1),1) Cov]);
    end
    img_meas =R;
    
    clear S
    for u=1:size(beh_meas,2)
        [B,BINT,S(:,u)] = regress(beh_meas(:,u),[ones(size(beh_meas,1),1) Cov]);
    end
    beh_meas = S;
end

temp1 = img_meas;
temp2 = beh_meas;

% PCA    
if pca_ind ==1
    [COEFF,SCORE, latent,TSQUARED, EXPLAINED]=pca(img_meas);
    PCA_rank=cumsum(latent)./sum(latent);
    ind1=min(find(PCA_rank>pca_percentage));
    temp1=SCORE(:,1:ind1);
    
    cca_output.PCA_coeff = COEFF;
    cca_output.PCA_rank = PCA_rank;
end

%% CCA
[A,B,r,U,V,stats]=canoncorr(temp1,temp2);
cca_output.A=A;
cca_output.B=B;
cca_output.r=r;
cca_output.U=U;
cca_output.V=V;
cca_output.stats=stats;

component_r = r(component_num);

%% permutation test
% permutation method 1
for i=1:perm_time
    rand_id2=randperm(size(temp2,1));
    rand_sub_phen=temp2(rand_id2,:);
    [rand_A,rand_B,rand_r,rand_U,rand_V,rand_stats]=canoncorr(temp1,rand_sub_phen);
    rand_all(i)=rand_r(component_num);
end
permutation_P=length(find(rand_all>r(component_num)))/perm_time;
cca_output.permutation_p=permutation_P;

if perm_pfwer==1
% permutation method 2
[pfwer,r2,A2,B2,U2,V2] = permcca(temp2,temp1,perm_time);
% AM Winkler, O Renaud, SM Smith, TE Nichols NIH - Univ. of Geneva - Univ. of Oxford
cca_output.permutation_pfwer=pfwer;
end

%% feature weight calculation
[weight_img_r weight_img_p]=corr(U(:,component_num),img_meas);
[weight_beh_r weight_beh_p]=corr(V(:,component_num),beh_meas);
cca_output.weight_img_r=weight_img_r;
cca_output.weight_beh_r=weight_beh_r;
save(cca_output_path,'cca_output');

fprintf('CCA analysis finished at %d/%02d/%02d %02d:%02d:%02d \n',c);
disp('==================================================================');
end