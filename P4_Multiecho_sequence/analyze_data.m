%% 
% Load data
twix =mapVBVD('C:\Users\salom\OneDrive - epfl.ch\Bureau\phd\pulseq_workshop\github\pulseq-ekin\Pulseq_Workshop\P4_Multiecho_sequence\meas_MID00058_FID319207_FatSeperators.dat');
%%
% Look at the data
twix{:};
% Keep the second struct (no noise there)
twix=twix{2};
% image info
twix.image
