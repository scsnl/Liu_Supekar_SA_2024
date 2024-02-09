function BN2brainmap(brainmeasure,mask_path,output_name)

% V = spm_vol(fullfile(filesep,'Volumes','menon','projects','jinliu5','2021_Longt_math_gene','data','imaging','roi','Reslice_BN_Atlas_246_3mm.nii'));
% V = spm_vol(fullfile(filesep,'oak','stanford','groups','menon','projects','jinliu5','2021_Longt_math_gene','data','imaging','roi','BN_Atlas_246_3mm.nii'));
% V = spm_vol(fullfile('Y:','projects','jinliu5','2021_Longt_math_gene','data','imaging','roi','BN_Atlas_246_3mm.nii'));
 V = spm_vol(fullfile(mask_path));

[Y,XYZmm] = spm_read_vols(V);
Y_new=Y;
for i=1:246
Y_new(find(Y==i))=brainmeasure(i);
end

V.fname = output_name;
V.dt = [64 0];
V.private.dat.fname = V.fname;  
spm_write_vol(V,Y_new); 
% y_Write(Y_new,V,output_name);

end
