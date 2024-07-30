function distMatrix = AggregateDist(dist, opt, th)
% This code returns an aggregate function of two distance matrices.
% dist1 and dist2 must be of the same size.
% opt : aggregate function option
%       (1:Min, 2:Max, 3:Threshold-count, 4:Sum (or Average), 5:Geometric mean, 6:Harmonic mean)
% th: threshold
% When opt is outside the range, this function returns dist{1} as dist.
% 
% [numOfMetrics,~] = size(dist);
numOfMetrics = length(dist);
distMatrix = dist{1};
if numOfMetrics == 1
    opt = 0;
end
if opt == 1
    for i = 2:numOfMetrics
        distMatrix = min(distMatrix,dist{i});
    end
elseif opt == 2
    for i = 2:numOfMetrics
        distMatrix = max(distMatrix,dist{i});
    end
elseif opt == 3
    distMatrix = zeros(length(distMatrix));
    for i = 1:numOfMetrics
        distMatrix = distMatrix + (dist{i} >= th).*dist{i};
    end
elseif opt == 4
    for i = 2:numOfMetrics
        distMatrix = (distMatrix + dist{i});
    end
elseif opt == 5
    for i = 2:numOfMetrics
        distMatrix = distMatrix.*dist{i};
    end
    distMatrix = distMatrix.^(1/numOfMetrics);
elseif opt == 6
    distMatrix = zeros(length(distMatrix));
    for i = 1:numOfMetrics
        distMatrix = distMatrix + 1./(dist{i});
    end
    distMatrix = numOfMetrics./distMatrix;
end



