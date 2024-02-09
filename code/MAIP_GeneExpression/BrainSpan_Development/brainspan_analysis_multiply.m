function brainspan_analysis_multiply(PLSweighted_genelist,input_expression_matrix,input_gene_rows,input_sample_columns,output_path)


[C,IA,IB]=intersect(cell2table(input_gene_rows.gene_symbol),PLSweighted_genelist(:,1));

overlap_expression_matrix=input_expression_matrix(IA,2:end);
overlap_PLSweighted = PLSweighted_genelist(IB,3);

sample_age_summary=tabulate(input_sample_columns.age);

category_sample_new{1,1}={'8 pcw','9 pcw','12 pcw','13 pcw','16 pcw','17 pcw','19 pcw','21 pcw','24 pcw','25 pcw','26 pcw','35 pcw','37 pcw'}; %  prenatal
category_sample_new{2,1}={'4 mos','10 mos','1 yrs'}; % infancy
category_sample_new{3,1}={'2 yrs','3 yrs','4 yrs'}; % early childhood
category_sample_new{4,1}={'8 yrs','11 yrs','13 yrs'}; % late childhood
category_sample_new{5,1}={'15 yrs','18 yrs','19 yrs'}; % Adolescence
category_sample_new{6,1}={'21 yrs','23 yrs','30 yrs','36 yrs','37 yrs','40 yrs'}; % Adulthood

clear LIA top_mean_line bottom_mean_line category_sample

sample_region_summary_name=tabulate(input_sample_columns.structure_name);
sample_region_summary_index=tabulate(input_sample_columns.structure_acronym);

category_sample_region=[];
category_sample_region{1,1}='AMY'; % amygdala
category_sample_region{2,1}='MD'; % thalamus
category_sample_region{3,1}='STR'; % striatum
category_sample_region{4,1}='HIP'; % hippocampus
category_sample_region{5,1}='MFC'; % frontal
category_sample_region{6,1}='OFC'; % frontal
category_sample_region{7,1}='DFC'; % frontal
category_sample_region{8,1}='VFC'; % frontal
category_sample_region{9,1}='M1C'; % frontal-motor
category_sample_region{10,1}='S1C'; % parietal
category_sample_region{11,1}='IPC'; % parietal
category_sample_region{12,1}='A1C'; % temporal
category_sample_region{13,1}='STC'; % temporal
category_sample_region{14,1}='ITC'; % temporal
category_sample_region{15,1}='V1C'; % occipital

category_sample_region_name = {'Amgdala';'Thalamus';'Striatum';'Hippocampus';...
    'Medial prefrontal cortex';'Orbital frontal cortex';'Dorsolateral prefrontal cortex';'Ventrolateral prefrontal cortex';...
    'Primary motor cortex';'Primary somatosensory cortex';'Inferior parietal cortex';...
    'Primary auditory cortex';'Posterior superior temporal cortex';'Inferolateral temporal cortex';...
    'Primary visual cortex'};

for j=1:length(category_sample_new)
    for l=1:length(category_sample_region)
        [LIA]=find(ismember(input_sample_columns.age,category_sample_new{j}) & ismember(input_sample_columns.structure_acronym,category_sample_region{l,1}));
        mean_line(j,l) = sum(mean(table2array(overlap_expression_matrix(:,LIA)),2).*table2array(overlap_PLSweighted));
        category_sample(LIA,1) = j;
        category_sample(LIA,2) = l;
    end
end

colors = [178 91 38;...
    180 160 155;...
    105,63,153;...
    203,179,214;...
    247,128,31;...
    252,244,156;...
    228,31,38;...
    244,153,154;...
    55,160,72;...
    181,216,140;...
    96,185,168;...
    33,121,180;...
    166,206,227;...
    201,201,201;...
    122,122,122]./255

lineall=plot(mean_line,'LineWidth',3);
for i=1:length(colors)
    lineall(i).Color=colors(i,:);
end
set(gca,'XTick',[1:1:length(category_sample_new)]);
set(gca,'YLim',[0 max(max(mean_line))+max(max(mean_line)).*0.1]);
set(gca,'XTickLabel',{'Prenatal','Infancy','Early childhood',' Late childhood','Adolescence','Adulthood'});
ylabel('Gene Score');
set(gca,'FontSize',17);
set(gca,'LineWidth',2.5);
lgd=legend(category_sample_region_name,'Location','eastoutside');
lgd.Box = 'off';
set(gcf,'Position',[10,100,1180,350]);
print(gcf,'-dtiff',output_path);

% imagesc(top_mean_line)
% set(gca,'XTick',[1:1:length(category_sample_region)])
% set(gca,'XTickLabel',category_sample_region)
% set(gca,'YTick',[1:1:length(category_sample_new)])
% set(gca,'YTickLabel',{'Prenatal','Infancy','Early childhood','Late childhood','Adolescence','Adulthood'})
% set(gca,'FontSize',15)
% set(gcf,'position',[100,100,770,300]);
% print(gcf,'-dtiff',fullfile(output_path,'PLS1top_mean_brainspan.tif'));





