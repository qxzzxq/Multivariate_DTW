format long
%
% Parameters
%---------------------------------------------------------------------
%adjustRI = 1; % 1: use adjusted rand index, 0: use (regular) rand index. 
cutPercentile = 60; % Cut percentile used in hierarchical clustering
numberDesiredClusters = 6; % Number of clusters to output
%----------------------------------------------------------------------
loadDist = 1; % Distance matrices option (0: Calculate new distance matrices, 1: Load pre-calculated)
%fileNameDist = 'dist_xxxxxx_xxxxxx.mat'; %this file should contain cell
%arrays.

DataName = 'Data_BU';

if loadDist
    fileNameDist = fullfile([pwd '/' DataName '/dist_All_metrics_trimmed_062023_185124_w5min.mat']);
else
    fileNameData = strcat([DataName '/dat_All_trimmed_041723.mat']);
    fileNameDist = fullfile([pwd '/' DataName '/dist_All_trimmed_' datestr(now,'mmddyy_HHMMSS')]);
end

fileNameRItable = fullfile([pwd '/' DataName '/RI_table_' datestr(now,'mmddyy')]);
distName = ''; % name of the variable in which distance matrices are stored. The default name is 'dist'.
fileNameOut = 'RI'; % Output file name
%----------------------------------------------------------------------
pre_processing = 1; % Pre-processing option (1: None, 2: zNorm, 3: 0/1 norm, 4: filter)
mono_norm = 1; % 1: Unify normalization option for all metric, 0: find all possible combinations
%----------------------------------------------------------------------
distanceType = 1; % Distance matrix method option (1: Dynamic Time Warping, 2: Minkowski)
WindowSize = 120; % Either absolute or percentage. (EX: 5(%), 120 (abs))
WindowOption = 2; % DTW windowSize option. (1: percentage, 2: Absolute time range)
Minkowski_p = 2; % The order of Minkowski metric. (p >= 1)
%----------------------------------------------------------------------
post_processing = 1; % Post-processing option (1: None, 2: Length norm, 3: Complexity norm)
%----------------------------------------------------------------------
agg = 0; % Aggregate function (1: Min, 2: Max, 3: Threshold-count, 4: Average, 5: Sum, 6: Geometric mean, 7: Harmonic mean)
th = 0.5; % Threshold for count in aggregate function
%----DB SCAN -----------------------------------------------------------
e = 0.5; % epsilon: neighborhood search radius 
minpts = 5; % minimum number of objects considered as a cluster
%----------------------------------------------------------------------
% Setting number of preprocesses, distance methods, postprocesses, and
% aggregation methods. 
numPre = 7;
numDist = 2;
numPost = 3;
numAgg = 6;

clear titleString;
titleString.pre = ["None" "Z-norm" "0/1 norm" "-1/1 norm" "Decimal norm" "Tanh norm" "Sigmoid norm"];
titleString.dist = ["DTW" "Minkowski"];
titleString.post = ["None" "PathLength norm" "Complexity norm"];
titleString.agg = ["Min" "Max" "Threshold Count" "Sum" "Geometric Mean" "Harmonic Mean"];

%-----------------------------------------------------
% Calculate or Load distance matrices
%
if loadDist == 1
    load(fileNameDist);
%{
    if isempty(distName)
        dist = load(fileNameDist).dist;
        correctionFactor = load
    else
        dist = load(fileNameDist).(distName);
    end
    %}
    [~, ~, numOfMetrics] = size(dist); 
else
    % ----------Calculate Distance matrices of all the metrics-------------  
    load(fileNameData);
    dat_orig = dat;
    numOfMetrics = length(dat);
    %current_date_time = datestr(now,'mmddyy_HHMMSS');
    % Format: dist = CalculateDistanceMatrix(dat, pre, post, windowOption, windowSize);

    postdat = cell(1,numPre);
    dist = cell(numPre, numDist, numOfMetrics);
    pathLength = cell(numPre, numDist, numOfMetrics);
    correctionFactor = cell(numPre, numDist, numOfMetrics);
    
    
    distTime = zeros(numPre, numDist)-1;

    for pre = 1:numPre
        %------------- Pre-Processing -----------------------------------------
        for i = 1:numOfMetrics
            dat{i} = PreProcessDat(dat{i}, pre);
        end
        %----------------------------------------------------------------------
        postdat{pre} = dat;
        dat = dat_orig;

        for distOpt = 1:numDist

            %-------------Calculate Distance matrices ----------------------------
            tic

            for i = 1:numOfMetrics
                display(['pre ' num2str(pre), ', distOpt ' num2str(distOpt) ', metric ' num2str(i)])
                [pathLength{pre,distOpt,i}, correctionFactor{pre,distOpt,i}, dist{pre,distOpt,i}] = CalculateDistanceMatrix(postdat{pre}{i}, distOpt, WindowOption, WindowSize, Minkowski_p);
                %{
                if (distOpt ~= 1)
                    pathLength{pre,distOpt,i} = [];
                end
                %}
                pl = pathLength{pre,distOpt,i};
                cf = correctionFactor{pre,distOpt,i};
                distMat = dist{pre,distOpt,i};
                clear pl cf distMat
            end
            distTime(pre, distOpt) = toc;
        end
    end    
    %save(fileNameDist, 'pathLength','correctionFactor','dist','labels','classes', '-v7.3');
end

%--------------------------------------------------------------------------    
%%
% Post processing
distPost = cell(numPre, numDist, numPost,numOfMetrics);
for pre = 1:numPre
    for distOpt = 1:numDist
        for post = 1:numPost
            %--------------- Post-Processing ------------------------------------------
            for i = 1:numOfMetrics
                if (distOpt ~= 1) && (post == 2) 
                    continue
                end
                distPost{pre,distOpt,post,i} = PostProcessDist(dist{pre,distOpt,i}, post, pathLength{pre,distOpt,i}, correctionFactor{pre,distOpt,i});
            end
        end
    end
end

%%
% Aggregation

if mono_norm
    RI_table = zeros(6*numPre, numAgg);
    numComb = numPre;
    combList = zeros(numPre,numOfMetrics);
    for i = 1:numPre
        combList(i,:) = i;
    end
else 
    numComb = numPre^numOfMetrics;
    combList = combs(1:numPre,numOfMetrics);
end

% Initialize tables of Rand Index (RI), Adjusted RI (ARI), and Adjusted
% Mutual Information (AMI).
RI = zeros(numDist, numPost, numAgg)-1;
RI_max = -1;
ARI = zeros(numDist, numPost, numAgg)-1;
ARI_max = -1;
AMI = zeros(numDist, numPost, numAgg)-1;
AMI_max = -1;

ind_temp = cell(1,3);
for pre = 1:numComb
    fprintf('pre: %d\n', pre)
    for distOpt = 1:numDist
        for post = 1:numPost

            pairDistance = cell(numOfMetrics,1);
            for i = 1:numOfMetrics
                pairDistance{i} = distPost{combList(pre,i),distOpt,post,i};
            end
        %--------------- Aggrgate function ------------------------------------------
            for agg = 1:numAgg

                if (~isempty(find(combList(pre,:) == 1, 1)) && (agg == 3))||((distOpt ~= 1) && (post == 2))
                    continue
                end
                
                %{
                try
                    aggDist = AggregateDist(pairDistance, agg, th);
                    predictedClusters = Cluster_Hierarchical(labels, aggDist, numberDesiredClusters);
                    %predictedClusters = CalculateRI_DBscan(labels, aggDist, e, minpts);
                    %predictedClusters = Cluster_DBscan_exhaustive(labels, aggDist, minpts, numberDesiredClusters);
                    
                    [~,RI(distOpt,post,agg)] = rand_index(labels, predictedClusters);
                    [~,ARI(distOpt,post,agg)] = rand_index(labels, predictedClusters, 'adjusted');
                    AMI(distOpt,post,agg) = ami(labels, predictedClusters);
                    %[~ ,RI(distOpt,post,agg)] = CalculateRI_DBscan(labels, aggDist, e, minpts);                    
                catch
                    warning('Aggregate function - (%s) is not valid (pre:%d, dist:%d, post:%d)', titleString.agg(agg), pre, distOpt, post);
                end
                %}
                
                aggDist = AggregateDist(pairDistance, agg, th);
                predictedClusters = Cluster_Hierarchical(labels, aggDist, numberDesiredClusters);
                %predictedClusters = CalculateRI_DBscan(labels, aggDist, e, minpts);
                %predictedClusters = Cluster_DBscan_exhaustive(labels, aggDist, minpts, numberDesiredClusters);
                
                % Calculate evaluations (RI, ARI, AMI).
                [~,RI(distOpt,post,agg)] = rand_index(labels, predictedClusters);
                [~,ARI(distOpt,post,agg)] = rand_index(labels, predictedClusters, 'adjusted');
                AMI(distOpt,post,agg) = ami(labels, predictedClusters);
                %[~ ,RI(distOpt,post,agg)] = CalculateRI_DBscan(labels, aggDist, e, minpts);                  

            end
        end
    end
    j = (numPost*numDist)*pre-numPost;
    for i = 1:numDist
        RI_table(j-numPost+1:j,:) = squeeze(RI(i,:,:));
        ARI_table(j-numPost+1:j,:) = squeeze(ARI(i,:,:));
        AMI_table(j-numPost+1:j,:) = squeeze(AMI(i,:,:));
        
        j = j+numPost;
    end

    % Update Rand Index
    RI_test = reshape(RI,[],1);
    [RI_max_temp, RI_ind] = max(RI_test);
    [ind_temp{:}] = ind2sub(size(RI), RI_ind);
    if RI(ind_temp{:}) == RI(RI_ind) && RI(ind_temp{:}) > RI_max
        RI_max = RI(ind_temp{:});
        ind1 = pre;
        [ind2, ind3, ind4] = ind_temp{:};
    end
    RI = zeros(numDist, numPost, numAgg)-1;
    
    % Update Adjusted Rand Index
    ARI_test = reshape(ARI,[],1);
    [ARI_max_temp, ARI_ind] = max(ARI_test);
    [ind_temp{:}] = ind2sub(size(ARI), ARI_ind);
    if ARI(ind_temp{:}) == ARI(ARI_ind) && ARI(ind_temp{:}) > ARI_max
        ARI_max = ARI(ind_temp{:});
        ind11 = pre;
        [ind22, ind33, ind44] = ind_temp{:};
    end
    ARI = zeros(numDist, numPost, numAgg)-1;
    
    % Update Adjusted Mutual Information
    AMI_test = reshape(AMI,[],1);
    [AMI_max_temp, AMI_ind] = max(AMI_test);
    [ind_temp{:}] = ind2sub(size(AMI), AMI_ind);
    if AMI(ind_temp{:}) == AMI(AMI_ind) && AMI(ind_temp{:}) > AMI_max
        AMI_max = AMI(ind_temp{:});
        ind111 = pre;
        [ind222, ind333, ind444] = ind_temp{:};
    end
    AMI = zeros(numDist, numPost, numAgg)-1;

end

RI_max_index = [ind1, ind2, ind3, ind4];
pre_metric_RI = combList(ind1,:);

ARI_max_index = [ind11, ind22, ind33, ind44];
pre_metric_ARI = combList(ind11,:);

AMI_max_index = [ind111, ind222, ind333, ind444];
pre_metric_AMI = combList(ind111,:);
%%
if mono_norm  
    RI_cell = num2cell(RI_table);
    ind_na = find(RI_table < 0);
    [RI_cell{ind_na}] = deal('N/A');
    
    ARI_cell = num2cell(ARI_table);
    ind_na = find(ARI_table < 0);
    [ARI_cell{ind_na}] = deal('N/A');   
    
    AMI_cell = num2cell(AMI_table);
    ind_na = find(AMI_table < 0);
    [AMI_cell{ind_na}] = deal('N/A');   
    
    fileNameRItable = append(fileNameRItable,'_mono');
    save(fileNameRItable, 'RI_cell', 'RI_max_index', 'RI_max', 'pre_metric_RI',...
        'ARI_cell', 'ARI_max_index', 'ARI_max', 'pre_metric_ARI',...
        'AMI_cell', 'AMI_max_index', 'AMI_max', 'pre_metric_AMI');
else
    save(fileNameRItable, 'RI_max_index', 'RI_max', 'pre_metric_RI',...
        'ARI_max_index', 'ARI_max', 'pre_metric_ARI',...
        'AMI_max_index', 'AMI_max', 'pre_metric_AMI');
end




%%
disp('>>>>>>  Start tabulating')


for i = 1:numPre
    for ii = 1:numPre
        for j = 1:numDist
            RI_table = [RI_table; squeeze(RI(i,ii,j,:,:))];
        end
    end
end
%}
RI_cell = num2cell(RI_table);
ind_na = find(RI_table < 0);
[RI_cell{ind_na}] = deal('N/A');
[row, col] = find(RI_table == max(max(RI_table)));


%%
%save(fileNameRItable, 'RI_cell', 'RI_table');
%save(fullfile([pwd '/' DataName '/Time_table_' datestr(now,'mmddyy')]), 'distTime', 'time_cell', 'time_table');
