function R_upper = R_upper(maps1,maps2,Npoints_list,Nresamp,sig,filtflag)
% R_upper = R_upper(maps1,maps2,Npoints_list,Nresamp,sig,filtflag)
% Resamples maps1 and computes the correllation with maps2
% maps1 - The maps to be resampled
% maps2 - The maps that the resampled maps will be compared against
% Npoints_list - number of resampling points to be drawn for each map
% Nresamp - number of times to resample each map
% sig - standard deviation of gaussian kernel for filtering
% filtflag - true if gaussian kernel is to be applied to resampled maps

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


R_upper= zeros(Nresamp,length(maps1));
for pic = 1:length(maps1)
    n = Npoints_list(pic); %number of tap points for this image
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
        R_upper(k,pic) = R(1,2);
    end
    disp([pic Rtrue mean(R_upper(:,pic)) Rtrue-mean(R_upper(:,pic))]);
end
