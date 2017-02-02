function map_correlation(binsize,filtflag,sig,firstInterestPoint, rem_list)
% Find the correlations between different types of maps

%% Parameters
if nargin==0
    binsize = 32; %must be a factor of both dimensions of the images (power of 2)
    filtflag = false;
    sig = 27;
    firstInterestPoint = true;
end

%% load data

% load('fixation_maps.mat') %variable: fixmaps
% load('interest_maps_firstclick.mat') %variable: interest_maps
% load('tap_maps.mat')      %variable: tap_maps
load('salmaps.mat');      %variable: salmaps
load('russmaps.mat'); %variable: russmaps


load('fixation_points.mat');

if firstInterestPoint
    load('interest_points_firstclick.mat');
else
    load('interest_points.mat');
end

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
russmaps = russmaps(list);

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

[X, Y] = meshgrid(-3*sig:3*sig,-3*sig:3*sig);
gauss = 1*exp(-((X.^2)+(Y.^2))./((2*sig)^2));



%% Remove center bias (average over all maps of a given type);
remove_center_bias = false;
shift_bias_correct = false;

for pic = 1:length(fixmaps)
    fixmaps{pic} = downsize_map(fixmaps{pic},binsize);
    interest_maps{pic} = downsize_map(interest_maps{pic},binsize);
    tap_maps{pic} = downsize_map(tap_maps{pic},binsize);
    salmaps{pic} = downsize_map(salmaps{pic},binsize);
    russmaps{pic} = downsize_map(russmaps{pic},binsize);
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

%% Show maps of an example image
inx = 71;
inx = find(list==(inx-30),1);
figure(1)
subplot(231)
imagesc(fixmaps{inx}); colormap('gray');
title('Fixation map');
subplot(232)
imagesc(interest_maps{inx}); colormap('gray');
title('Interest map');
subplot(233)
imagesc(tap_maps{inx}); colormap('gray');
title('Tap map');
subplot(234)
imagesc(salmaps{inx}); colormap('gray');
title('Saliency map (Itti 1998)');
subplot(235)
imagesc(russmaps{inx}); colormap('gray');
title('Saliency map (Russell et al., 2014)'); 
boldify

%% Compute True correlations

R_true = zeros(5,5,length(fixmaps));
for pic = 1:length(fixmaps)
    fix = fixmaps{pic}(:);
    int = interest_maps{pic}(:);
    tap = tap_maps{pic}(:);
    sal = salmaps{pic}(:);
    russ = russmaps{pic}(:);
    R_true(:,:,pic) = corrcoef([fix int tap sal russ]);
end

Rfixint = squeeze(R_true(1,2,:));
Rfixtap = squeeze(R_true(1,3,:));
Rinttap = squeeze(R_true(2,3,:));
Rtapsal = squeeze(R_true(3,4,:));
Rfixsal = squeeze(R_true(1,4,:));
Rfixruss = squeeze(R_true(1,5,:));
Rintsal = squeeze(R_true(2,4,:));


save('corr_true','R_true','Rfixint','Rfixtap','Rinttap','Rtapsal','Rfixsal','Rfixruss','Rintsal','remove_center_bias','shift_bias_correct');

Ravg = mean(R_true,3);

%% Plot all correlations

% for pic = 1:length(fixmaps)
%     figure(1);
%     plot(fixmaps{pic}(:),tap_maps{pic}(:),'b.','MarkerSize',0.5)
% 
%     pause
% end

%% Compute correlations from mismatched images
Rfixint_rand = zeros(length(fixmaps),length(fixmaps));
Rfixtap_rand = zeros(length(fixmaps),length(fixmaps));
Rinttap_rand = zeros(length(fixmaps),length(fixmaps));
Rfixsal_rand = zeros(length(fixmaps),length(fixmaps));
Rtapsal_rand = zeros(length(fixmaps),length(fixmaps));
Rfixruss_rand = zeros(length(fixmaps),length(fixmaps));
Rintsal_rand = zeros(length(fixmaps),length(fixmaps));

for pic1 = 1:length(fixmaps)
    for pic2 = 1:length(fixmaps)
        fix = fixmaps{pic1}(:);
        int = interest_maps{pic2}(:);
        R = corrcoef(fix,int);
        Rfixint_rand(pic1,pic2) = R(1,2);
        
        fix = fixmaps{pic1}(:);
        tap = tap_maps{pic2}(:);
        R = corrcoef(fix,tap);
        Rfixtap_rand(pic1,pic2) = R(1,2);
        
        int = interest_maps{pic1}(:);
        tap = tap_maps{pic2}(:);
        R = corrcoef(int,tap);
        Rinttap_rand(pic1,pic2) = R(1,2);
        
        tap = tap_maps{pic1}(:);
        sal = salmaps{pic2}(:);
        R = corrcoef(tap, sal);
        Rtapsal_rand(pic1,pic2) = R(1,2);
        
        fix = fixmaps{pic1}(:);
        sal = salmaps{pic2}(:);
        R = corrcoef(fix,sal);
        Rfixsal_rand(pic1,pic2) = R(1,2);
        
        fix = fixmaps{pic1}(:);
        russ = russmaps{pic2}(:);
        R = corrcoef(fix,russ);
        Rfixruss_rand(pic1,pic2) = R(1,2);

        int = interest_maps{pic1}(:);
        sal = salmaps{pic2}(:);
        R = corrcoef(int,sal);
        Rintsal_rand(pic1,pic2) = R(1,2);

%         fix1 = fixmaps{pic1}(:);
%         fix2 = fixmaps{pic2}(:);
%         R = corrcoef(fix1,fix2);
%         Rfixfix_rand(pic1,pic2) = R(1,2);
    end
end

save('corr_rand','Rfixint_rand','Rfixtap_rand','Rinttap_rand','Rtapsal_rand', 'Rfixsal_rand', 'Rfixruss_rand','Rintsal_rand');

%% Generate "upper bound" correllation of interest and fixation based on resampling from interest maps

R_upper_intfix = R_upper(interest_maps,interest_maps,Npoints_fix,Nresamp,sig,filtflag);
save('R_upper_intfix','R_upper_intfix');


%% Generate "upper bound" correllation of taps and fixations based on resampling from fixation maps

R_upper_fixtap = R_upper(fixmaps,fixmaps,Npoints_tap,Nresamp,sig,filtflag);
save('R_upper_fixtap','R_upper_fixtap');


%% Generate "upper bound" correllation of interest and taps based on resampling from interest maps

R_upper_inttap = R_upper(interest_maps,interest_maps,Npoints_tap,Nresamp,sig,filtflag);
save('R_upper_inttap','R_upper_inttap');

%% Correllate Interest and Fixation maps using only Ntaps number of interest points

R_intNtaps_fix = R_upper(interest_maps,fixmaps,Npoints_tap,Nresamp,sig,filtflag);
save('R_intNtaps_fix2','R_intNtaps_fix');

%% Generate "upper bound" correlation of taps to salience based on resampling from saliency maps

R_upper_tapsal = R_upper(salmaps,salmaps,Npoints_tap,Nresamp,sig,filtflag);
save('R_upper_tapsal','R_upper_tapsal');

%% Generate "upper bound" correlation of fixations to salience based on resampling from saliency maps
R_upper_fixsal = R_upper(salmaps,salmaps,Npoints_fix,Nresamp,sig,filtflag);
save('R_upper_fixsal','R_upper_fixsal');

%% Generate "upper bound" correlation of fixations to russell salience based on resampling from saliency maps
R_upper_fixruss = R_upper(russmaps,russmaps,Npoints_fix,Nresamp,sig,filtflag);
save('R_upper_fixruss','R_upper_fixruss');

%%
R_upper_intsal = R_upper(salmaps,salmaps,Npoints_int,Nresamp,sig,filtflag);
save('R_upper_intsal','R_upper_intsal');