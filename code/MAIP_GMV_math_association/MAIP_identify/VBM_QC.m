clear,clc
%% QC for VBM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script is used to copy pdf from VBM
%
% Jin Liu
% 6/9/2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting path
% imac
box_path = fullfile(filesep,'Users','jinjin','Library','CloudStorage','Box-Box','Jin Liu','2021 Longt math gene');
oak_path = fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene');

% windows
% box_path = fullfile('C:','Users','jinliu5','Box','Jin Liu','2021 Longt math gene');
% oak_path = fullfile('Y:','projects','jinliu5','2021_Longt_math_gene');

%% copy pdf
output_folder = fullfile(oak_path,'results','smri','vbm','Stanford','wjiii','pdf_vbm');
mkdir(output_folder)
subjectlist = readtable(fullfile(oak_path,'data','subjectlist','wjiii','subjectlist_multiple_times_N213.csv'));
missing_list =[];
for i=1:size(subjectlist)
    input_path = fullfile(oak_path,'data','imaging','participants',num2str(subjectlist.PID(i),'%04.f'),['visit' num2str(subjectlist.Visit(i))],['session' num2str(subjectlist.Session(i))],'anatomical','VBM_spm12','report',['catreport_' subjectlist.Var4{i} '.pdf']);
    if exist(input_path)
        output_path = fullfile(output_folder,[num2str(subjectlist.PID(i),'%04.f'),'_visit' num2str(subjectlist.Visit(i)), '_catreport_' subjectlist.Var4{i} '.pdf'])
        copyfile(input_path,output_path)
    else
        temp = subjectlist.PID(i);
        missing_list = [missing_list temp];
    end
end

% 582 and 8008 fail to pass VBM pipeline


cd(output_folder)