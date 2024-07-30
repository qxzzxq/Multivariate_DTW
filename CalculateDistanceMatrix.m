function [pathLength, CF, dist] = CalculateDistanceMatrix(dat, distanceType, windowOption, windowSize, Minkowski_p) 
% This function calculates DTW distance matrices.
%         dat  : data set in form of cell to be used which contains all the
%                time series data of an attribute. 
% distanceType : Method of calculating a Distance matrix
%               (1:DTW, 2:Minkowski)
% windowOption : DTW windowSize option. 
%               (1: percentage, 2: Absolute time range in the same unit to that of time series data)
% windowSize   : 
% Minkowski_p  : Mincowski parameter


if ~(distanceType==1 || distanceType==2)
    distanceType = 1;
end

numOfJobs = length(dat);

p_rec = 1/Minkowski_p; % reciprocal of p
CE = cellfun(@(x) sqrt(sum(diff(x).^2)), dat, 'UniformOutput',false); % Complexity Estimate
    
% Calculate dist
dist = zeros(numOfJobs);
pathLength = zeros(numOfJobs);
CF = zeros(numOfJobs); % Correction Factor
if distanceType == 1
    if windowOption == 1
        for ii = 1:numOfJobs
            %disp(ii)
            for jj = ii+1:numOfJobs
                m = min(length(dat{ii}),length(dat{jj}));
                A = [];
                B = [];

                A = dat{ii};
                B = dat{jj};            

                out = dtwCPP(A,B,ceil(m*windowSize/100));

                dist(ii,jj) = out(1);
                dist(jj,ii) = dist(ii,jj);
                %dist(ii,jj) = 1-((out(1))^2/(2*out(2))); 
                pathLength(ii,jj) = out(2);
                pathLength(jj,ii) = pathLength(ii,jj);

                CF(ii,jj) = max(CE{ii},CE{jj})/min(CE{ii},CE{jj});
                CF(jj,ii) = CF(ii,jj);
            end

            dist(ii,ii) = 0;    
        end
    else
        for ii = 1:numOfJobs
            %%disp(ii)
            for jj = ii+1:numOfJobs
%                 m = min(length(dat{ii}),length(dat{jj}));
%                 
%                 if windowSize > m
%                     windowSize = m;
%                 end
                
                A = [];
                B = [];

                A = dat{ii};
                B = dat{jj};            

                out = dtwCPP(A,B,windowSize);

                dist(ii,jj) = out(1);
                dist(jj,ii) = dist(ii,jj);
                %dist(ii,jj) = 1-((out(1))^2/(2*out(2)));        
                pathLength(ii,jj) = out(2);
                pathLength(jj,ii) = pathLength(ii,jj);

                CF(ii,jj) = max(CE{ii},CE{jj})/min(CE{ii},CE{jj});
                CF(jj,ii) = CF(ii,jj);
            end

            dist(ii,ii) = 0;
        end
    end
else
    for ii = 1:numOfJobs
        %%disp(ii)
        for jj = ii+1:numOfJobs
            %Take the smaller series as reference
            m = min(length(dat{ii}) , length(dat{jj}));
            A = [];
            B = [];
            
            %Get downsampled versions of the time series   
            A = dat{ii}(1:m);
            B = dat{jj}(1:m);

            %Euclidean distance
            dist(ii,jj) = (sum((abs(A - B)).^Minkowski_p)).^p_rec;
            dist(jj,ii) = dist(ii,jj);
            
            CF(ii,jj) = max(CE{ii},CE{jj})/min(CE{ii},CE{jj});
            CF(jj,ii) = CF(ii,jj);
        end
        dist(ii,ii) = 0;
    end
end


