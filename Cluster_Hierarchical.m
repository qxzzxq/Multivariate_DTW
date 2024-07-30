function C = Cluster_Hierarchical(labels, distMatrix, nC)
%------------------------------------------------------------------
% labels - Ground truth labels 
% dist - n by n distance matrix
% nC - number of desired clusters
%------------------------------------------------------------------

%distMatrix = aggDist; nC = numberDesiredClusters;
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



%[~,RI] = rand_index(labels, C);
%%

%{
labelList = unique(labels);
indexGT = cell(1,length(labelList));



label_text = ["ExaMiniMD" "LAMMPS" "sw4lite" "sw4" "SWFFT" "HACC"]';
%
%label_text = ["normal" "anomaly"]';
labels_str = strings(size(labels));
for i = 1:length(labelList)
    ind = find(labels == labelList(i));
    %labels_str(ind) = classes(i);
    labels_str(ind) = labelList(i);
end




%labels_str = num2str(labels);
for i = 1:length(labelList)
    indexGT{i} = find(labels == labelList(i));
end
newcolor=[0 0 0];
colororder(newcolor)
indexClusters = clusters;

z = linkage(squareform(distMatrix));
prc = prctile(z(:,3),60);
figure
subplot(2,1,1)
%H = dendrogram(z,0,'ColorThreshold',prc);
H = dendrogram(z,0);
set(H, 'Color', newcolor);
ClusterColoring(indexGT, length(labelList), labels_str, 1);
%ClusterColoring(indexGT, length(labelList), labels, 1);
t = title('Ground Truth');
set(t,'Interpreter','none')

subplot(2,1,2)
%{
H = dendrogram(z,0,'ColorThreshold',optimalCutoff);
H = dendrogram(z,0,'ColorThreshold',cutoff(optimalIndex+1)); % same as default coloring
%}
%H = dendrogram(z,0,'ColorThreshold',min(cutoff(end), optimalIndex+1));
H = dendrogram(z,0); set(H, 'Color', newcolor);
%H = dendrogram(z,0); set(H, 'Color', newcolor); ClusterColoring(indexClusters, numClustersfound, C, 0);
%
t = title('Clusters found');
set(t,'Interpreter','none')

titleText = DataName;
%titleText = 'Data_BU (8nodes)';
%titleText = [DataName '  (var1: 0/1 norm, var2: Z-norm)'] ;
t = sgtitle(sprintf(titleText));
set(t,'Interpreter','none')
t.FontSize = 16;

%}

%%
%plot(cutoff, numOfIndividualClusters,'.b')
%{
figure
histogram(aggregate_anomaly_index)
xlabel('index(job)')
ylabel('counts')
title('kernel anomaly counts')
%}