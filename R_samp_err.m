function R_se = R_samp_err(maps1,maps2,Npoints_list,Nresamp,sig,filtflag)
% R_se = R_samp_err(maps1,maps2,Npoints_list,Nresamp,sig,filtflag)
% Generates new maps by resampling with replacement from maps1 and computes
% the correllation with maps2.
%
%
% maps1 - The maps to be resampled
% maps2 - The maps that the resampled maps will be compared against
% Npoints_list - number of resampling points to be drawn for each map
% Nresamp - number of times to resample each map (usually 1000)
% sig - standard deviation of gaussian kernel for filtering
% filtflag - true if Gaussian kernel is to be applied to resampled maps
%
% If you want to compare this to the notation in the paper (latex expression below),
% then maps1 and maps2 are equal (P), and Npoints_list is the number of
% points sampled in Q.
%
% $R(\hat{P},\tilde{P}^{Q})$
%
% example: R_samp_err_intfix = R_samp_err(interest_maps,interest_maps,Npoints_fix,1000,27,False);
%
% gives $R(\hat{I},\tilde{I}^{F})$

%
% By Daniel Jeck 2015.

%% handle missing parameters
% Set filter flag
if nargin < 5
    filtflag = false;
end
if nargin ==5
    filtflag = true;
end

%convert singleton to list
if length(Npoints_list)==1
    Npoints_list = Npoints_list*ones(size(maps1));
end


if filtflag %make the gaussian kernel if you're using it
    [X, Y] = meshgrid(-3*sig:3*sig,-3*sig:3*sig);
    gauss = 1*exp(-((X.^2)+(Y.^2))./((2*sig)^2));
end

%% Compute distribution of correlations under the sample error hyothesis for each image
R_se= zeros(Nresamp,length(maps1));
for pic = 1:length(maps1)
    n = Npoints_list(pic); %number of points for this image
    map1 = maps1{pic};
    map2 = maps2{pic};
    pdf = map1/sum(map1(:));
    R = corrcoef(map1(:),map2(:)); 
    Rtrue = R(1,2);
    
    rand_inx = discretesample(pdf,n*Nresamp);
    rand_inx = reshape(rand_inx,[n,Nresamp]);
    [Y, X] = ind2sub(size(map1),rand_inx);
    
    
    
    for k = 1:Nresamp
        im = accumarray([Y(:,k) X(:,k)],1,size(map2));
        if filtflag
            im = conv2(im,gauss,'same');
        end
%         if pic == 64-30
%             keyboard;
%         end
        R = corrcoef(im(:),map2(:));
        R_se(k,pic) = R(1,2);
    end
%     disp([pic Rtrue mean(R_se(:,pic)) Rtrue-mean(R_se(:,pic))]);
end
