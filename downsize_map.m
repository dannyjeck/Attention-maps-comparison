function result = downsize_map(map,binsize)
%result = downsize_map(map,binsize)
%Takes in a map and downsamples it in binsize x binsize squares. Output matrix is
%composed of sums of the binsize x binsize blocks in map.


n = size(map,1)/binsize;
m = size(map,2)/binsize;

subs = kron(reshape(1:(n*m), n, m),ones(binsize)); 
result = reshape(accumarray(subs(:), map(:)), n, m);
