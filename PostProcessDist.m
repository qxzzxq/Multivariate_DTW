function newDist = PostProcessDist(dist, post, pathLength, CF) 
% This function normalize a distance matrix.
% post_processing = 0; % Post-processing option (0: None, 1: Length norm, 2: Complexity norm)
%         dist : distance matrices in form of cell where each cell contains
%         a symmetric distance matrix.
%         post : post-processing option 
%               (1: None, 2: Length normalization, 3:Complexity normalization)
%--------------------------------------------------------------------------

if post == 2
    numOfJobs = length(dist);
    newDist = dist./sqrt(pathLength);
    for i = 1:numOfJobs
        newDist(i,i) = 0;
    end
elseif post == 3
    newDist = dist.*CF;
else
    newDist = dist;
end


