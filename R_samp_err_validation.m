% Script to validate the sample error hypothesis test (Figure S2)
%
% By Daniel Jeck 2017

%% Parameters
binsize = 64; %must be a factor of both dimensions of the images (power of 2)
filtflag = false;
sig = 27;
firstInterestPoint = true;
rem_list = [];

%% load data

% load('fixation_maps.mat') %variable: fixmaps
% load('interest_maps_firstclick.mat') %variable: interest_maps
% load('tap_maps.mat')      %variable: tap_maps
load('salmaps.mat');      %variable: salmaps


load('fixation_points.mat');


load('interest_points.mat');


load('tap_points.mat');



%% select appropriate natural scenes
tap_maps = tap_maps_points(31:78); 
interest_maps = interest_maps_points([1:25 76:98]);
fixmaps = fixmaps_points([1:25 76:98]);

% tap_maps = tap_maps(31:78); 
% interest_maps = interest_maps([1:25 76:98]);
% fixmaps = fixmaps([1:25 76:98]);


fix_points_x = fix_points_x([1:25 76:98]);
tap_points = tap_points(31:78);
all_x_points_int = all_x_points_int([1:25 76:98]);

list = 1:48;
list([rem_list]-30) = [];

tap_maps = tap_maps(list);
interest_maps = interest_maps(list);
fixmaps = fixmaps(list);
fix_points_x = fix_points_x(list);
tap_points = tap_points(list);
all_x_points_int = all_x_points_int(list);
salmaps = salmaps(list);

%% Downsample maps

remove_center_bias = false;
shift_bias_correct = false;

for pic = 1:length(fixmaps)
    fixmaps{pic} = downsize_map(fixmaps{pic},binsize);
    interest_maps{pic} = downsize_map(interest_maps{pic},binsize);
    tap_maps{pic} = downsize_map(tap_maps{pic},binsize);
    salmaps{pic} = downsize_map(salmaps{pic},binsize);
end

tap_mean = zeros(size(tap_maps{1}));
int_mean = zeros(size(interest_maps{1}));
fix_mean = zeros(size(fixmaps{1}));
for pic = 1:length(fixmaps)
    tap_mean = tap_mean + tap_maps{pic};
    fix_mean = fix_mean + fixmaps{pic};
    int_mean = int_mean + interest_maps{pic};
end
tap_mean = tap_mean/length(fixmaps);
fix_mean = fix_mean/length(fixmaps);
int_mean = int_mean/length(fixmaps);


if remove_center_bias
    for pic = 1:length(fixmaps)
        fixmaps{pic} = fixmaps{pic} - fix_mean;
        interest_maps{pic} = interest_maps{pic} - int_mean;
        tap_maps{pic} = tap_maps{pic} - tap_mean;
    end
end

if shift_bias_correct
    load('bias_calc.mat');
    for pic = 1:length(fixmaps)
        tap_maps{pic} = circshift(tap_maps{pic},round([-bias_all(2) bias_all(1)]));
    end
end


%% Upper bound analysis params
Nresamp=1000; %number of times to resample the appropriate map

Npoints_tap = zeros(size(tap_points));
Npoints_int = zeros(size(tap_points));
Npoints_fix = zeros(size(tap_points));
for pic= 1:length(fixmaps)
    Npoints_fix(pic) = length(fix_points_x{pic});
    Npoints_int(pic) = length(all_x_points_int{pic});
    Npoints_tap(pic) = length(tap_points(pic).Xdata);
    
end


%%
Ntrials = 1000;
p_Rfixtap_samp_err = zeros(Ntrials,1);
p_Rfixtap_samp_err2 = zeros(Ntrials,1);

for k = 1:Ntrials
    R_samp_err_fixtap = R_samp_err(fixmaps,fixmaps,Npoints_tap,Nresamp,sig,filtflag);
    Rfixtap = R_samp_err(fixmaps,fixmaps,Npoints_tap,1,sig,filtflag); %only one sample for the data to test against the null
    Rfixtap = Rfixtap(:);
    
    [~, p_Rfixtap_samp_err(k)] = ztest(mean(Rfixtap),mean(R_samp_err_fixtap(:)), ...
        sqrt(var(Rfixtap(:))/length(Rfixtap) + var(mean(R_samp_err_fixtap))/length(mean(R_samp_err_fixtap))),0.95,'both');

    
    p_Rfixtap_samp_err2(k) = (sum(mean(Rfixtap)>mean(R_samp_err_fixtap,2))+1)/(Nresamp+1);

    
    disp(['k = ' num2str(k) '; pfixint = ' num2str(p_Rfixtap_samp_err(k)) '/' num2str(p_Rfixtap_samp_err2(k))]);
    
end

figure(1)
hist(p_Rfixtap_samp_err,20);
figure(2)
hist(p_Rfixtap_samp_err2,20);