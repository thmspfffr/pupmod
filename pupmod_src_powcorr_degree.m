%% pupmod_src_powcorr_degree
% COMPUTES THE NUMBER OF ALTERED CORRELATIONS AS A FUNCTION OF
% CARRIER FREQUENCY.

clear

addpath /home/gnolte/meg_toolbox/toolbox/
addpath /home/gnolte/meg_toolbox/fieldtrip_utilities/
addpath /home/gnolte/meg_toolbox/toolbox_nightly/
addpath /home/gnolte/meg_toolbox/meg/
addpath ~/pconn/matlab/
addpath ~/Documents/MATLAB/Colormaps/Colormaps' (5)'/Colormaps/

v = 23;


outdir   = '/home/tpfeffer/pupmod/proc/conn/';

addpath ~/pconn/matlab/
SUBJLIST  = [4 5 6 7 8 9 10 11 12 13 15 16 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34];

load(sprintf('~/pupmod/proc/conn/pupmod_src_powcorr_cleaned_v%d.mat',v));

f = fopen('~/Documents/MATLAB/aal_symm.nii.txt','rt');
aal_labels = textscan(f,'%d %s %d','headerlines',0);

load aalmask_grid_coarse.mat


