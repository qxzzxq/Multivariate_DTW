clear
clc
close all


DataName = 'Data_BU'; app_type = 'sw4'; ano_type = 'membw'; fileNameData = fullfile([pwd '\' DataName '\' app_type '\' ano_type '\dat_All_' app_type '_samples_120623.mat']);
% DataName = 'Data_EMPIRE'; fileNameData = fullfile('C:\Users\Solji\Dropbox\Research\ModularCode\Data_EMPIRE\slidingW_5s\dat_All_slidingW_5s_36ts_102423.mat');

distanceType = 1; % Distance matrix method option (1: Dynamic Time Warping, 2: Minkowski)
WindowSize = 120; % Either absolute or percentage. (EX: 5(%), 120 (abs))
WindowOption = 2; % DTW windowSize option. (1: percentage, 2: Absolute time range)
Minkowski_p = 2; 
numberDesiredClusters = 6;

%%
load(fileNameData);
%----------------------------- 
% set timeseries length: t_length 
t_length = 600; 

% Select one anomaly(index:5) and the rest normal(6:end) 
%index_wrt_full = 1:15;
index_wrt_full = [3,6:15];
%index_wrt_full = 5:15;

for i=1:length(dat)
    dat{i} = dat{i}(index_wrt_full);
    for j=1:length(index_wrt_full)
        if length(dat{i}{j}) > t_length
            dat{i}{j}=dat{i}{j}(1:t_length);
        end
    end
end
%dat{1} = dat{1}(index_wrt_full);
jobID = jobID(index_wrt_full);
labels = labels(index_wrt_full);
labels_num = labels_num(index_wrt_full);

%-----------------------------
dat_orig = dat;

% sample - kernel function 268 
% dat = dat(37);

% sample - kernel function 18 
% dat = dat(17);



%
%--------------------------single metric
dat = cell(1);
%dat{1} = dat_orig{1};  metric_type = "Active Mem"; %set which metric to investigate
dat{1} = dat_orig{2};  metric_type = "CPU Util"; 
dat_orig = dat;
%----------------------------- 
%}

numOfMetrics = length(dat);



%%



%%%
th = 0.05;
mono_norm = 1;

numPre = 7;
numDist = 1;
numPost = 3;
numAgg = 6;


postdat = cell(1,numPre);
dist = cell(numPre, numDist, numOfMetrics);
pL = cell(numPre, numDist, numOfMetrics);
cF = cell(numPre, numDist, numOfMetrics);
max_cost = cell(numPre, numDist, numOfMetrics);
t_test = cell(numPre, numDist, numOfMetrics);
t_ref = cell(numPre, numDist, numOfMetrics);

%%%
for pre = 1:numPre
    %------------- Pre-Processing -----------------------------------------
    for i = 1:numOfMetrics
        %dat{i} = PreProcessDat(dat{i}, pre_processing);
        dat{i} = PreProcessDat(dat{i}, pre);
    end
    %----------------------------------------------------------------------
    postdat{pre} = dat;
    dat = dat_orig;

    for distOpt = 1:numDist

        %-------------C alculate Distance matrices ----------------------------
        tic

        for i = 1:numOfMetrics
            display(['pre ' num2str(pre), ', distOpt ' num2str(distOpt) ', metric ' num2str(i)])
            [pL{pre,distOpt,i}, cF{pre,distOpt,i}, dist{pre,distOpt,i},t_test{pre,distOpt,i}, t_ref{pre,distOpt,i}, max_cost{pre,distOpt,i}] = CalculateDistanceMatrix_max_cost(postdat{pre}{i}, distOpt, WindowOption, WindowSize,Minkowski_p);
            %{
            if (distOpt ~= 1)
                pathLength{pre,distOpt,i} = [];
            end
            %}
            pl = pL{pre,distOpt,i};
            cf = cF{pre,distOpt,i};
            distMat = dist{pre,distOpt,i};
            %save(fullfile([pwd '\distDir\' 'dist_pre_' num2str(pre) '_distOpt_' num2str(distOpt) '_metric_' num2str(i) '_' datestr(now,'mmddyy')]), 'pl','cf','distMat');
            clear pl cf distMat
        end
    end
end    

%%% Post processing
distPost = cell(numPre, numDist, numPost,numOfMetrics);
costPost = cell(numPre, numDist, numPost,numOfMetrics);
for pre = 1:numPre
    for distOpt = 1:numDist
        for post = 1:numPost
            %--------------- Post-Processing ------------------------------------------
            for i = 1:numOfMetrics
                if (distOpt ~= 1) && (post == 2) 
                    continue
                end
                distPost{pre,distOpt,post,i} = PostProcessDist(dist{pre,distOpt,i}, post, pL{pre,distOpt,i}, cF{pre,distOpt,i});
                costPost{pre,distOpt,post,i} = PostProcessDist(dist{pre,distOpt,i}, post, pL{pre,distOpt,i}, cF{pre,distOpt,i});
            end
        end
    end
end

%%% Aggregation
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

RI = zeros(numDist, numPost, numAgg)-1;
RI_max = -1;
ARI = zeros(numDist, numPost, numAgg)-1;
ARI_max = -1;
AMI = zeros(numDist, numPost, numAgg)-1;
AMI_max = -1;

ind_temp = cell(1,3); 
%pre = 2; distOpt = 1; post = 1; agg = 1;
%% Single combination test
pre = 1; distOpt = 1; post = 1; agg = 1;
pairDistance = cell(numOfMetrics,1);
pairCost = cell(numOfMetrics,1)
for i = 1:numOfMetrics
    pairDistance{i} = distPost{combList(pre,i),distOpt,post,i};
    pairCost{i} = distPost{combList(pre,i),distOpt,post,i};
end
aggDist = AggregateDist(pairDistance, agg, th);
                  
pre_type=["none", "Z-norm", "MinMax", "[-1/1]", "Decimal", "Tanh", "Sigmoid"];
post_type=["none", "length norm", "complexity"];

%%
%------------------------------------------------------------------
% labels - Ground truth labels 
% dist - n by n distance matrix
% nC - number of desired clusters
%------------------------------------------------------------------

distMatrix = aggDist; nC = numberDesiredClusters;
%distMatrix_orig = distMatrix;
numOfJobs = length(distMatrix);
%jobIndex = 1:numOfJobs;


if nC > numOfJobs
    nC = numOfJobs;
end



numOfJobs = length(distMatrix);
z = linkage(squareform(distMatrix));
%inconsist = inconsistent(z);

if max(z(:,3)) == 0
    cutoff = linspace(0,1,10);
else
    %cutoff = linspace(0,max(z(:,3)),10000);
    cutoff = [0;nonzeros(z(:,3))];
    numOfIndividualClusters = zeros(length(cutoff),1);
    numOfIndividualClusters(1) = numOfJobs;
    numOfClusters = zeros(length(cutoff),1);
end

%{
if max(z(:,3))==0
    continue
end
%}

for i = 2:length(cutoff)
    %Cluster the matrix
    if isnan(cutoff(i))
        cutoff(i)=Inf;
    end
    C = cluster(z,'Cutoff',cutoff(i),'Criterion','distance');

    %Process the clusters and print them
    uv = unique(C);
    edges = uv(1)-0.5:uv(end)+0.5;
    n = histcounts(C,edges);
    numOfIndividualClusters(i) = length(find(n == 1));
    % find the number of clusters that have 5 <= jobs. 
    numOfClusters(i) = length(find(n >= 5));
    %[~, I] = sort(n,'descend');
end


%{
figure
subplot(1,2,1)
%plot(cutoff, numOfIndividualClusters,'.b')
%plot(numOfIndividualClusters,cutoff,'.b')
plot(numOfClusters,cutoff,'.b')
%title('Kernel function ', num2str(kernelList(ii)));
ylabel('Cutoff height');
xlabel('Number of groups');
ylim([0,max(z(:,3))])
%}

%{
if max(numOfClusters) >= nC
    optimalIndex = find(numOfClusters == nC, 1, 'last');
    if isempty(optimalIndex)
        optimalIndex = find(numOfClusters == min(numOfClusters(numOfClusters > nC)), 1, 'last');
    end
    optimalCutoff = cutoff(optimalIndex);
else
    optimalCutoff = cutoff(find(numOfClusters == max(numOfClusters), 1, 'last'));
end
%}


nC_tempList = numOfClusters(numOfClusters >= nC);
if isempty(nC_tempList)
    optimalIndex = find(numOfClusters == max(numOfClusters), 1, 'last');
    %optimalCutoff = cutoff(find(numOfClusters == max(numOfClusters), 1, 'last'));
    optimalCutoff = cutoff(optimalIndex);
    optimalFlag = 0;
else 
    nC_temp = min(nC_tempList);
    %optimalCutoff = cutoff(find(numOfClusters == nC_temp, 1, 'last'));
    optimalFlag = 1;
end

%
while optimalFlag
    optimalIndex = find(numOfClusters == nC_temp, 1, 'last');
    if ~isempty(optimalIndex) && numOfIndividualClusters(optimalIndex) <= (0.2)*numOfIndividualClusters(1)
        optimalCutoff = cutoff(optimalIndex);
        optimalFlag = 0;
    elseif isempty(nC_tempList)
        optimalCutoff = cutoff(find(numOfClusters == max(numOfClusters(numOfClusters<=nC)), 1, 'last'));
        optimalFlag = 0;
    else
        nC_tempList(find(nC_tempList == nC_temp, 1, 'last')) = [];
        nC_temp = min(nC_tempList);
        if isempty(nC_tempList)
            optimalCutoff = cutoff(find(numOfClusters == max(numOfClusters(numOfClusters<=nC)), 1, 'last'));
            optimalFlag = 0;
        end
    end
end
%}


%Cluster the matrix with the largest 
C = cluster(z,'Cutoff',optimalCutoff,'Criterion','distance');
%{
subplot(1,2,2)
%H = dendrogram(z,0,'ColorThreshold',optimalCutoff);
H = dendrogram(z,0,'ColorThreshold',cutoff(optimalIndex+1));
title('Clusters found (predicted) with optimal cutoff');
xlabel('Index');
ylabel('Distance');
ylim([0,max(z(:,3))])

titleText = DataName;
%titleText = 'Data_BU (8nodes)';
t = sgtitle(sprintf(titleText));
set(t,'Interpreter','none')
t.FontSize = 16;

%}

uv = unique(C);
edges = uv(1)-0.5:uv(end)+0.5;
countC = histcounts(C,edges);
countC(countC < 5) = NaN;

numClustersfound = length(countC(~isnan(countC)));
Cnew = 1:length(C);
Cnew = Cnew'+ length(countC)*ones(length(C),1);  
clusters = cell(1, numClustersfound);
for i = 1:numClustersfound
    ind = find(C == find(countC == max(countC),1,'first'));
    Cnew(ind) = i;
    clusters{i} = ind;
    countC(find(countC == max(countC),1,'first')) = NaN;
end
C = Cnew;



%%
labelList = unique(labels);
indexGT = cell(1,length(labelList));

labels_str = strings(size(labels));
for i = 1:length(labelList)
    ind = find(labels == labelList(i));
    labels_str(ind) = labelList(i);
end

for i = 1:length(labelList)
    indexGT{i} = find(labels == labelList(i));
end

% Change the order of classes and the color assignment(red:anomaly)
indexGT = flip(indexGT); 


newcolor=[0 0 0];
colororder(newcolor)
indexClusters = clusters;

z = linkage(squareform(distMatrix));
prc = prctile(z(:,3),60);
[H, ~, leaf] = dendrogram(z,0);
set(H, 'Color', newcolor);

ClusterColoring(indexGT, length(labelList), labels_str, 1);
set(gca, 'XTickLabel', index_wrt_full(leaf));

%[t,s] = title(sprintf('BU: Active Mem (pre:%s/post:%s) %d min(s)', pre_type(pre), post_type(post), t_length/60),'Ground Truth');
[t,s] = title(sprintf('BU: %s (pre:%s/post:%s) %d min(s)', metric_type, pre_type(pre), post_type(post), t_length/60),'Ground Truth');
%[t,s] = title(sprintf('BU: Aggregate (pre:%s/post:%s) %d min(s)', pre_type(pre), post_type(post), t_length/60),'Ground Truth');

set(t,'Interpreter','none')
t.FontSize = 16;


%% plot samples
%{
pre = 1;
figure
hold on
for i = 1:length(index_wrt_full)
    if (i < 6)
        plot(postdat{pre}{1}{i}, 'r')
    %{
    elseif (i == 10)
        plot(postdat{pre}{1}{i}, 'b')
    %}
    else
        plot(postdat{pre}{1}{i}, 'g')
    end
end

%title(sprintf('ExaMiniMD - Active Mem (pre:%s/post:%s)', pre_type(pre), post_type(post)), sprintf('normal(green), %s(red)', ano_type));
title(sprintf('%s - %s (pre:%s/post:%s)', app_type, metric_type, pre_type(pre), post_type(post)), sprintf('normal(green), %s(red)', ano_type));

%}
