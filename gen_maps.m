% script to turn the various point data sets into (fixation/interest/tap) maps

%% Generate Fixation maps from Chris's study (Masciocchi et al. 2009)
sig = 27; %The standard deviation of the gaussian
fixmaps = cell(100,1);
fixmaps_points = cell(100,1);
imsize = [768 1024];

%load in the eye data
load_eye_data=load('revised_eye_data.mat');
eye_data=load_eye_data.revised_final_data;
clear load_eye_data

%remove bad data
for picture=1:100
    for ss=1:21
        %I'm not sure what this does besides remove NaNs (copied from Chris code)
        %Based on his analysis of fixation data by fixation number, it
        %looks like the NaNs are just junk (rather than lost data orsomething)
        actual_data = eye_data{picture}{ss}(:,1)>0; 
        eye_data{picture}{ss}=eye_data{picture}{ss}(actual_data,:);
        
        if FM_firstsec & ~isempty(eye_data{picture}{ss})
            tot_time = cumsum(eye_data{picture}{ss}(:,1));
            sel = tot_time < 1; %save first second of data
            lastinx = find(~sel,1);
            lastdata = eye_data{picture}{ss}(lastinx,:);
            if lastinx ~=1
                lastdata(1) = 1-tot_time(lastinx-1); %modify last entry so only first second counts
            elseif ~isempty(lastdata)
                lastdata(1) = 1;
            end
            
            eye_data{picture}{ss} = [eye_data{picture}{ss}(sel,:); lastdata];
        end
    end
end

[X, Y] = meshgrid(-3*sig:3*sig,-3*sig:3*sig);
gauss = 1*exp(-((X.^2)+(Y.^2))./((2*sig)^2));

fix_points_x = cell(100,1);
fix_points_y = cell(100,1);
fix_time = cell(100,1);

%generate fixation maps
for pic = 1:100
    pic_data = vertcat(eye_data{pic}{:}); %concatenate subject data
    y = round(pic_data(:,3));
	x = round(pic_data(:,2));
    sel = x>=1 & y>=1 & x<=imsize(2) & y<=imsize(1); %select only valid fixations
    fix_points_x{pic} = x(sel);
    fix_points_y{pic} = y(sel);
    fix_time{pic} = pic_data(sel,1);
    if FM_by_time
        t = fix_time{pic};
    else
        t = 1;
    end
    fixmaps_points{pic} = accumarray([y(sel) x(sel)],t,imsize);
    fixmaps{pic} = conv2(fixmaps_points{pic},gauss,'same');
end

save('fixation_maps','fixmaps');
save('fixation_points','fix_points_x','fix_points_y','fix_time','fixmaps_points');

%% Generate Interest maps from Chris's study (Masciocchi et al. 2009)
% sig = 27; %The standard deviation of the gaussian
interest_maps = cell(100,1);
interest_maps_points = cell(100,1);
imsize = [768 1024];

load_int_data=load('int_data1.mat');

% Around 33 pix=1 degree visual angle
% Chris code to get interest points out of the data file
all_x_points_int = cell(100,1);
all_y_points_int = cell(100,1);
for pic=1:100
    int_data_begin=load_int_data.int_data{1,pic};
    p=length(int_data_begin);
    
    all_x_points_int{pic}(1:p,1)=round(int_data_begin(:,6)*1.6);
    if ~firstInterestPoint
        all_x_points_int{pic}(p+1:2*p,1)=round(int_data_begin(:,9)*1.6);
        all_x_points_int{pic}((2*p)+1:3*p,1)=round(int_data_begin(:,12)*1.6);
        all_x_points_int{pic}((3*p)+1:4*p,1)=round(int_data_begin(:,15)*1.6);
        all_x_points_int{pic}((4*p)+1:5*p,1)=round(int_data_begin(:,18)*1.6);
    end
    all_y_points_int{pic}(1:p,1)=round(int_data_begin(:,7)*1.6);
    if ~firstInterestPoint
        all_y_points_int{pic}(p+1:2*p,1)=round(int_data_begin(:,10)*1.6);
        all_y_points_int{pic}((2*p)+1:3*p,1)=round(int_data_begin(:,13)*1.6);
        all_y_points_int{pic}((3*p)+1:4*p,1)=round(int_data_begin(:,16)*1.6);
        all_y_points_int{pic}((4*p)+1:5*p,1)=round(int_data_begin(:,19)*1.6);
    end
end

clear load_int_data

[X, Y] = meshgrid(-3*sig:3*sig,-3*sig:3*sig);
gauss = 1*exp(-((X.^2)+(Y.^2))./((2*sig)^2));

%generate interest maps
for pic = 1:100
    y = double(all_y_points_int{pic});
	x = double(all_x_points_int{pic});
    sel = x>=1 & y>=1 & x<=imsize(2) & y<=imsize(1); %select only valid fixations
    all_x_points_int{pic} = x(sel);
    all_y_points_int{pic} = y(sel);
    interest_maps_points{pic} = accumarray([y(sel) x(sel)],1,imsize);
    interest_maps{pic} = conv2(interest_maps_points{pic},gauss,'same');
end

save('interest_maps','interest_maps');
save('interest_points','all_x_points_int','all_y_points_int','interest_maps_points');

%% Generate tap maps Jeck et al. (submitted Sept 2016 to Vision Research)

data_dir = './tap_data252/';
% sig = 27; %The standard deviation of the gaussian
imsize = [768 1024];

tap_maps = cell(80,1);
tap_maps_points = cell(80,1);


% Read all the data from text files
Nsubj = 252;
subj = cell(Nsubj,1);
for k = 1:Nsubj
    if k>=100
        fname = [data_dir 'subject_0' num2str(k) '.txt'];
    elseif k>=10
        fname = [data_dir 'subject_00' num2str(k) '.txt'];
    else
        fname = [data_dir 'subject_000' num2str(k) '.txt'];
    end
    fid = fopen(fname);
    data = fread(fid, '*char')'; %read all contents into data as a char array (don't forget the `'` to make it a row rather than a column).
    fclose(fid);
    subj{k} = regexp(data, ',', 'split'); %This will return a cell array with the individual entries for each string you have between the commas.
    subj{k} = reshape(subj{k}(2:end),7,[]);
    
end

% Organize onto the images
first_tap = false;
img_list = 31:80; %skip the 30 simple scenes
% Nimg = length(img_list);
image_data = repmat(struct('Xdata',[],'Ydata',[],'RT',[]),80,1);
for k = 1:Nsubj
    if first_tap;
        m = 2;
    else
        m = size(subj{k},2);
    end
    for l = 1:m
        im_idx = str2double(subj{k}{3,l});
        image_data(im_idx).Ydata(end+1) = str2double(subj{k}{4,l});
        image_data(im_idx).Xdata(end+1) = str2double(subj{k}{5,l});
        image_data(im_idx).RT(end+1)    = str2double(subj{k}{6,l});
    end
end

[X, Y] = meshgrid(-3*sig:3*sig,-3*sig:3*sig);
gauss = 1*exp(-((X.^2)+(Y.^2))./((2*sig)^2));

%generate tap maps
for pic = img_list
    %images were presented rotated 90 degrees, so x,y becomes (1024-y,x)
    y = image_data(pic).Xdata(:);
	x = imsize(2)-image_data(pic).Ydata(:);
    
    sel = x>=1 & y>=1 & x<=imsize(2) & y<=imsize(1); %select only valid fixations
    tap_maps_points{pic} = accumarray([y(sel) x(sel)],1,imsize);
    tap_maps{pic} = conv2(tap_maps_points{pic},gauss,'same');
end

tap_points = image_data;
save('tap_maps','tap_maps');
save('tap_points','tap_points','tap_maps_points');

%% Format Saliency maps generated using the virtualbox version of Itti's saliency
%http://ilab.usc.edu/toolkit/downloads-virtualbox.shtml

salmaps = cell(48,1);
% fixmaps_points = cell(100,1);
imsize = [768 1024];

for k = 1:length(salmaps)
    
    im = imread(['./sal_data/salMap' num2str(30+k) '-SM00.png']);
    salmaps{k} = double(imresize(imrotate(im,-90),imsize));
end

save('salmaps','salmaps')