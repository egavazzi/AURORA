% This script is here to control and call all the make_all_* scripts.

results_dir = "/mnt/data/etienne/Julia/AURORA.jl/data/Visions2/"
% RunDirs = {"Alfven_536s_correct_msis_and_scattering"}
RunDirs = {'Alfven_536s_correct_msis'}
%%
make_all_Ie_top
make_all_Ie_top_MI
%%
make_all_Q_lambda
make_all_I_lambda
%%
make_all_IQ_lambda_plots
%%
cxmax = 8
make_all_animations
