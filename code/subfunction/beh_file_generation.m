function beh_file_generation(beh_input_path,N,output_path)

%% loading data
disp('==================================================================');
fprintf('covariate files preparation\n');

% load behavior data
beh_meas=readtable(beh_input_path,'sheet',['Stanford_cohort_' N]);

%% data preparation
numopt_mathres_var=[beh_meas.numopsStd beh_meas.mathresStd];
save(fullfile(output_path,['numopt_mathres_' N '.mat']),'numopt_mathres_var');

% age
age_var=[beh_meas.age];
save(fullfile(output_path,['age_' N '.mat']),'age_var');

% gender
for i=1:length(beh_meas.gender)
    if strcmp(beh_meas.gender{i},'Male')
        gender_ind(i,1)=1;
    elseif strcmp(beh_meas.gender{i},'Female')
        gender_ind(i,1)=2;
    end
end
gender_var=[gender_ind];
save(fullfile(output_path,['gender_' N '.mat']),'gender_var');

% study
for i=1:length(beh_meas.study)
    if strcmp(beh_meas.study{i},'20: Math Longitudinal Development') | strcmp(beh_meas.study{i},'20: Math Longitudinal Development (20)')
        study_ind(i,1)=1;
    elseif strcmp(beh_meas.study{i},'24: Math 8-week Intervention') | strcmp(beh_meas.study{i},'24: Math 8-week Intervention (24)')
        study_ind(i,1)=2;
    elseif strcmp(beh_meas.study{i},'262: Math Whiz ASD/1-week Retrieval Intervention') | strcmp(beh_meas.study{i}, '262: Math Whiz ASD/1-week Retrieval Intervention (262)')
        study_ind(i,1)=3;
    elseif strcmp(beh_meas.study{i},'31.1: Auditory Autism Voice Identification') | strcmp(beh_meas.study{i},'31.2: Auditory Autism Voice Identification 1st Follow-Up')
        study_ind(i,1)=4;
    elseif strcmp(beh_meas.study{i},'38: Cross-sect Auditory Autism Voice Identification')
        study_ind(i,1)=5;
    elseif strcmp(beh_meas.study{i},'42: ASD Prosody') | strcmp(beh_meas.study{i},'42: ASD Prosody (42)')
        study_ind(i,1)=6;
    elseif strcmp(beh_meas.study{i},'43: ADHD Study') | strcmp(beh_meas.study{i},'43: ADHD Study (43)')
        study_ind(i,1)=7;
    end
end
study_var=[study_ind];
save(fullfile(output_path,['study_' N '.mat']),'study_var');

age_stu_var=[beh_meas.age study_ind];
save(fullfile(output_path,['age_stu_' N '.mat']),'age_stu_var');

fprintf('covariate files done\n');
disp('==================================================================');

end