clear,clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code was written by Dr Petra Vértes and is taken from Whitaker and Vértes, PNAS 2016 
% please cite that paper if you this code in your own work.
%
% Jin Liu modified to fit the data in this study
% 2022-3-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set the file path
current_path = pwd;
working_path = current_path(1:strfind(current_path,'code')-1);

addpath(genpath(fullfile(working_path,'code','subfunction')));
% addpath('spm12'));
% addpath('BrainNetViewer_20191031');

%% PLS 
% brain measures
math_related_brainmap=importdata(fullfile(working_path,'outputs','Stanford','CCA_math_brainmap_N219.txt'));

% gene expression
Table_gene=importdata(fullfile(working_path,'data','AllenBrainAtlas','BN246_geneexpression.csv'));

% the first column is index
gene_regional_expression=Table_gene.data(:,2:end);% this needs to be imported first

% gene name
genes=Table_gene.textdata(2:end);

% gene order
geneindex=1:size(gene_regional_expression,2);

% removed the region with NaN value of gene expression
num_regsion=find(~isnan(gene_regional_expression(:,1))); % only included region covered by donor samples

X=gene_regional_expression(num_regsion,:); % Predictors
Y=math_related_brainmap(num_regsion)'; % Response variable

% z-score:
clear gene_regional_expression Table_gene
X=zscore(X);
Y=zscore(Y);

%number of bootstrap iterations:
bootnum=1000;

% Do PLS in 3 dimensions (with 3 components):
dim=3;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions' IDs and weights in descending order of weight for both components:
[R,p]=corr([XS(:,1),XS(:,2),XS(:,3)],math_related_brainmap(num_regsion)');

%align PLS components with desired direction for interpretability 
if R(1,1)<0  %this is specific to the data shape we were using - will need ammending
    stats.W(:,1)=-1*stats.W(:,1);
    XS(:,1)=-1*XS(:,1);
end
if R(2,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,2)=-1*stats.W(:,2);
    XS(:,2)=-1*XS(:,2);
end
if R(3,1)<0 %this is specific to the data shape we were using - will need ammending
    stats.W(:,3)=-1*stats.W(:,3);
    XS(:,3)=-1*XS(:,3);
end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);
[PLS3w,x3] = sort(stats.W(:,3),'descend');
PLS3ids=genes(x3);
geneindex3=geneindex(x3);

% setting the output folder
output_path = fullfile(working_path,'outputs','GeneExpression')
mkdir(output_path)

%print out results
csvwrite(fullfile(output_path,'PLS1_ROIscores.csv'),XS(:,1));
csvwrite(fullfile(output_path,'PLS2_ROIscores.csv'),XS(:,2));
csvwrite(fullfile(output_path,'PLS3_ROIscores.csv'),XS(:,3));

%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];

%start bootstrap
for i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled subjects
    Yr=Y(myresample,:); % define X for resampled subjects
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
      
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run    
    
    temp=stats.W(:,3);%extract PLS3 weights
    newW=temp(x3); %order the newly obtained weights the same way as initial PLS 
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run   
end

%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');

%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';
temp3=PLS3w./PLS3sw';

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
[Z3 ind3]=sort(temp3,'descend');
PLS3=PLS3ids(ind3);
geneindex3=geneindex3(ind3);

%print out results
% later use first column of these csv files for pasting into GOrilla (for
% bootstrapped ordered list of genes) 
cd(output_path)
fid1 = fopen('PLS1_geneWeights_descend.csv','w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
end
fclose(fid1);

fid2 = fopen('PLS2_geneWeights_descend.csv','w');
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2);

fid3 = fopen('PLS3_geneWeights_descend.csv','w');
for i=1:length(genes)
  fprintf(fid3,'%s, %d, %f\n', PLS3{i},geneindex3(i), Z3(i));
end
fclose(fid3);

%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'ascend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);
[Z2 ind2]=sort(temp2,'ascend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
[Z3 ind3]=sort(temp3,'ascend');
PLS3=PLS3ids(ind3);
geneindex3=geneindex3(ind3);

% %print out results
% % later use first column of these csv files for pasting into GOrilla (for
% % bootstrapped ordered list of genes) 
% fid1 = fopen('PLS1_geneWeights_ascend.csv','w');
% for i=1:length(genes)
%   fprintf(fid1,'%s, %d, %f\n', PLS1{i}, geneindex1(i), Z1(i));
% end
% fclose(fid1);
% 
% fid2 = fopen('PLS2_geneWeights_ascend.csv','w');
% for i=1:length(genes)
%   fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
% end
% fclose(fid2);
% 
% fid3 = fopen('PLS3_geneWeights_ascend.csv','w');
% for i=1:length(genes)
%   fprintf(fid3,'%s, %d, %f\n', PLS3{i},geneindex3(i), Z3(i));
% end
% fclose(fid3);
