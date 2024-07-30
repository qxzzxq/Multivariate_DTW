function C = Cluster_DBscan_exhaustive(labels, distMatrix, k, nC)
%------------------------------------------------------------------
% labels - Ground truth labels 
% dist - n by n distance matrix
% k -  minimum number of objects considered as a cluster
%------------------------------------------------------------------

if k < 0
    minpts = 1;
else
    minpts = k;
end

if nC < 1
    nC = 5; % (default 5) when number of expected cluster is not provided.
end
    
%distMatrix = aggDist;
numOfJobs = length(distMatrix);
%jobIndex = 1:numOfJobs;


[RD,~,order]=optics_dist(minpts,distMatrix);
order_label = cell(numOfJobs,1);
for i = 1:numOfJobs
    order_label{i} = num2str(order(i));
end
RD_ordered = RD(order);
deltaRD = RD_ordered(1:end-1)-RD_ordered(2:end);
[~, ind] = sort(deltaRD,'descend');



%  Plot reachability of Ground Truth clusters
%{
load('colorOptions.mat')
colorOpt = [];
for i = 1:numOfJobs
    colorOpt(i,:) = colorOptions(labels(i)+1,:);
end



figure
subplot(2,1,1)
for i = 1:numOfJobs
    hold on
    %h(i)= plot(i,RD(order(i)),'o','Color',colorOpt(order(i),:));
    h(i)= line([i,i],[0,RD(order(i))],'Color',colorOpt(order(i),:),'Linewidth', 1);
end
%hold off
xlim([0, numOfJobs])
ylim([0,max(RD)])
xticks(1:numOfJobs);
xticklabels(order_label);
ylabel('Epsilon distance')

%[t,s] = title('Reachability', sprintf('Normalized, WindowSize=%dmin, k=%d',WindowSize/60, k));
t = title('Ground Truth', sprintf('minPts=%d', k));
%t.FontSize = 16;


labelList = unique(labels);
indexGT = cell(1,length(labelList));
for i = 1:length(labelList)
    indexGT{i} = find(labels == labelList(i));
end


% Add lengend to the reachability figure.
ind_leg = zeros(1,length(indexGT));
for i = 1:length(indexGT)
    ind_leg(i) = find(order==indexGT{i}(1));
end

%legend(h(ind_leg), string(labels(order(ind_leg))), 'Location','northwest')
legend(h(ind_leg), classes, 'Location','northwest')

hold off 

%}

% Find optimal epsilon distance e
%[pks, loc] = findpeaks(RD(order));
e = RD_ordered(ind(nC));
[RD_sorted, ~] = sort(RD_ordered,'descend');
j = find(RD_sorted == e,1);

optimalFlag = 1;

while optimalFlag
    idx = dbscan(distMatrix,e,minpts,'Distance','precomputed');

    uv = unique(idx);
    if length(unique(uv)) == 1 || e == min(RD)
        C = 1:length(idx);
        return
    else
        if uv(1) == -1        
            uv_outliers = uv(1);
            uv = uv(2:end); 
        end
        edges = uv(1)-0.5:uv(end)+0.5;
        countC = histcounts(idx,edges);
        [~,I] = sort(countC, 'descend');
        numClustersfound = length(I);
        if numClustersfound < nC
            j = j+1;
            e = RD_sorted(j);
            continue;
        end

        Cnew = 1:length(idx);
        %Cnew = Cnew'+ length(countC)*ones(length(idx),1);  % singletons 
        Cnew = Cnew'*(-1); 

        clusters = cell(1, numClustersfound);
        for i = 1:numClustersfound
            ind = find(idx == uv(I(i)));
            Cnew(ind) = i;
            clusters{i} = ind;
        end
        C = Cnew;
        optimalFlag = 0;
    end
end




% Plot reachability of predicted clusters
%{

colorOpt = zeros(numOfJobs,3);
%colorOptions = rand(numClustersfound,3);
colorOptions = hsv(numClustersfound);
for i = 1:numClustersfound
    ind = find(C == i);
    colorOpt(ind,1) = colorOptions(i,1);
    colorOpt(ind,2) = colorOptions(i,2);
    colorOpt(ind,3) = colorOptions(i,3);
end


subplot(2,1,2)
for i = 1:numOfJobs
    hold on
    %h(i)= plot(i,RD(order(i)),'o','Color',colorOpt(order(i),:));
    h(i)= line([i,i],[0,RD(order(i))],'Color',colorOpt(order(i),:),'Linewidth', 1);
end
%hold off
xlim([0, numOfJobs])
ylim([0,max(RD)])
xticks(1:numOfJobs);
xticklabels(order_label);
ylabel('Epsilon distance')
hold on
line([0,numOfJobs],[e,e],'LineStyle','--','Color','k');
t = title('Predicted clusters', sprintf('minPts=%d, and optimal eps=%f', k,e));
%t.FontSize = 16;

titleText = 'Reachability';
%titleText = 'Data_BU (8nodes)';
%titleText = [DataName '  (var1: 0/1 norm, var2: Z-norm)'] ;
t = sgtitle(sprintf(titleText));
set(t,'Interpreter','none')
t.FontSize = 16;

%}


end


