function C = Cluster_DBscan(labels, distMatrix, e, minpts)
%------------------------------------------------------------------
% labels - Ground truth labels 
% dist - n by n distance matrix
% e - epsilon: threshold for a neighborhood search radius
% minpts -  minimum number of neighbors
%------------------------------------------------------------------

%distMatrix_orig = distMatrix;
numOfJobs = length(distMatrix);
%jobIndex = 1:numOfJobs;

%e = 0.75;
idx = dbscan(distMatrix,e,minpts,'Distance','precomputed');

uv = unique(idx);
if length(unique(uv)) == 1
    C = 1:length(idx);
else
    if uv(1) == -1        
        uv_outliers = uv(1);
        uv = uv(2:end); 
    end
    edges = uv(1)-0.5:uv(end)+0.5;
    countC = histcounts(idx,edges);
    [~,I] = sort(countC, 'descend');
    numClustersfound = length(I);

    Cnew = 1:length(idx);
    Cnew = Cnew'+ length(countC)*ones(length(idx),1);  

    clusters = cell(1, numClustersfound);
    for i = 1:numClustersfound
        ind = find(idx == uv(I(i)));
        Cnew(ind) = i;
        clusters{i} = ind;
    end
    C = Cnew;
end
    
   
%[~,RI] = rand_index(labels, C);



%{

labelList = unique(labels);
indexGT = cell(1,length(labelList));

labels_str = strings(size(labels));
for i = 1:length(labelList)
    ind = find(labels == labelList(i));
    %labels_str(ind) = label_text(i);
    labels_str(ind) = labelList(i);
end



labels_str = num2str(labels);
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
H = dendrogram(z,0);
set(H, 'Color', newcolor);
ClusterColoring(indexClusters, numClustersfound, C, 0);
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







