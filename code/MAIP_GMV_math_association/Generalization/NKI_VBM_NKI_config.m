%-Configfile for preprocessmri.m
%__________________________________________________________________________

%-SPM version
paralist.spmversion = 'spm12';

%-Please specify parallel or nonparallel
%-e.g. for preprocessing and individualstats, set to 1 (parallel)
%-for groupstats, set to 0 (nonparallel)
paralist.parallel = '1';


%-Subject list (full path to the csv file)
paralist.subjectlist = '/oak/stanford/groups/menon/projects/jinliu5/2021_Longt_math_gene/data/subjectlist/subjectlist_NKI267.csv';
%---- example ----
%- PID, visit, session, spgrfilename
%- 7014, 1 ,1 , spgr

%- List of smri images 
% (no .nii or .img extensions)
% i.e.--> 
%      spgr (name of spgr file for the 1st subject in paralist.subjectlist)
%	   spgr (name of spgr file for the 2nd subject in paralist.subjectlist)
%	   spgr (name of spgr file for the 3rd subject in paralist.subjectlist)
%	   spgr_1 (name of spgr file for the 4th subject in paralist.subjectlist)
%      ....
paralist.smoothwidth= [12 12 12];

% I/O parameters
% - Raw data directory
paralist.rawdatadir = '/oak/stanford/groups/menon/rawdata/public/NKI-RS/';

% - Project directory - output of the preprocessing will be saved in the
% data/imaging folder of the project directory
paralist.projectdir = '/oak/stanford/groups/menon/projects/jinliu5/2021_Longt_math_gene/';

% 

