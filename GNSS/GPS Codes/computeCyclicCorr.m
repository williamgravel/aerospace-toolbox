function [corr] = computeCyclicCorr(x,y,opts)
arguments
    x (1,:) double
    y (1,:) double
    opts.Normalized (1,1) logical = false
end

len = length(x);
corr = zeros(1,len);

for n = 0:len-1
    corr(n+1) = sum(x.*circshift(y,n));
end

if opts.Normalized
    corr = corr/len;
end

end
