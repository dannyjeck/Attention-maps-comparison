%main_script.m runs the analysis of the fixation points, interest 
%points, tap points, and computed salience (following the itti et al.
%(1998) model) on the same set of natural scenes. This code is not 
%particularly clean but was used to generate Figures 3 and 5 in Jeck et
%al. (2017). 

%Code by Daniel Jeck 2015


FM_by_time = true; %Weight fixation points by fixation duration
firstInterestPoint = true; %Use only the first selection of the interest points or not.
filtflag = false; %Convolve all maps with a Gaussian (done in Masccioci et al. 2009 and others)
sig = 27; %Width of the Gaussian if applied
binsize = [64]; %initial map resolution (used in gen_maps). They are downsampled later
rem_list = [];
FM_firstsec = false; %flag to restrict fixation analysis to the first second

%%

gen_maps %Create the map files and saves them as .mat files. These may be useful for visualization.

% for binsize = [4 8 16 32 64 128 256] %binsize must be one of these values
for binsize = [64] %64 used in Figure 3, 256 used in figure 5
    map_correlation(binsize,filtflag,sig,firstInterestPoint,rem_list) % Computes all the relevant correlation valuse after downsampling the maps
    plot_corr
    h = figure(1);
    % pause;
    % saveas(h,['./Results/corrplot_' num2str(binsize) '_sal_FMbytime.fig']);
    % saveas(h,['./Results/corrplot_' num2str(binsize) '_sal_FMbytime.png']);


end
