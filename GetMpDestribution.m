function [GaussianMean,GaussianVariance] = GetMpDestribution(series)
    global seriesLength;
    global seriesCount;
    global dimension;
    distributions = 0;
    GaussianMean = zeros(1,seriesCount);
    GaussianVariance = zeros(1,seriesCount);
    series = cell2mat(series);
    for i = 1:seriesLength
        GaussianMean(1,i) = mean(series(:,i));
        GaussianVariance(1,i) = var(series(:,i));
    end
    
end

