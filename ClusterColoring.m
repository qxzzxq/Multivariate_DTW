function ClusterColoring(indexCluster, numC, labels, legendOpt)
% indexCluster is a cell array. Each cell of the array is a cluster
% containing the index for the cluster.
% numC is the number of clusters to be colored in the dendrogram. 
% For example, numC = 3, only first 3 clusters will be colored. 
% labels should be given an array whose element means the label for the
% index. The labels are used to generate a legend in the dendrogram.
colorOptions = load('colorOptions.mat').colorOptions;

numOfClusters = length(indexCluster);

if (numC == 0)
    numOfColors = 3;
elseif (numC > numOfClusters)
    numOfColors = numOfClusters;
else
    numOfColors = numC;
end

ax = gca; % get the axes handle

X = get(ax.Children,'XData'); % Get x values of all lines
Y = get(ax.Children,'YData'); % Get y values of all lines

lab = ax.XAxis.TickLabels; % get all the labels
[row, ~] = size(lab);

lab_int = [];

for i = 1:row
    lab_int(i) = str2num(lab(i,:));
end


X = cell2mat(X);
Y = cell2mat(Y);
indexToDelete=[];
for i = 1:length(Y)
    if sum(Y(i,:)) == 0
        indexToDelete=[indexToDelete, i];
    end
end

XY_ori = [X, Y];
XY = XY_ori;
XY = [XY(:,1), XY(:,4), XY(:,5), XY(:,6), XY(:,8)];
X_tailor=[];
Y_tailor=[];
XY = XY';
for i = 1:length(XY)
    if (XY(3,i) == 0)
        %{
        if(XY(5,i) == 0)
            X_tailor = horzcat(X_tailor, XY(1,i), XY(2,i));
            Y_tailor = horzcat(Y_tailor, XY(4,i), XY(4,i));
        else
            X_tailor = horzcat(X_tailor, XY(1,i));
            Y_tailor = horzcat(Y_tailor, XY(4,i));
        end
        %}
        if(XY(5,i) == 0 && XY(4,i) ~= 0)
            X_tailor = horzcat(X_tailor, XY(1,i), XY(2,i));
            Y_tailor = horzcat(Y_tailor, XY(4,i), XY(4,i));
        elseif (XY(5,i) ~= 0)
            X_tailor = horzcat(X_tailor, XY(1,i));
            Y_tailor = horzcat(Y_tailor, XY(4,i));
        end
    elseif (XY(5,i) == 0)
        X_tailor = horzcat(X_tailor, XY(2,i));
        Y_tailor = horzcat(Y_tailor, XY(4,i));
    end
end


row = min(row, length(X_tailor));
XY = [X_tailor',X_tailor',Y_tailor',zeros(row,1)];

%XY_sort = sortrows(XY, 1);
XY = sortrows(XY,1);




indexToCluster = zeros(1,row);
intermediateIndex = find(rem(XY(:,1),1) ~=0);
lab_ind_intermediate = lab_int;
if ~isempty(intermediateIndex)    
    lab_ind_intermediate(intermediateIndex+1) = [];  
end

for i = 1:numOfColors
    ind = find(ismember(lab_ind_intermediate,indexCluster{i}));
    indexToCluster(ind) = i;
end

nc = numOfColors+1;
%nc = numOfColors;
clr = cell(nc,1);
[numColor, ~] = size(colorOptions);

if (nc > numColor)
    for i = 1:numColor
        clr{i} = colorOptions(i,:);
    end
    for i = numColor:numOfColors
        clr{i+1} = rand(1,3);
    end
    
else 
    for i = 1:nc
        clr{i} = colorOptions(i,:);
    end
end
        
for i = 1:row
    h(i) = line(ax,XY(i,1:2),XY(i,3:4),'Color',clr{indexToCluster(i)+1},'Linewidth', 1);
    %line(ax,XY(i,1:2),XY(i,3:4),'Color',clr{indexToCluster(i)},'Linewidth', 1.2);
end
xlabel('Index');
%% Add lengend to the figure.
ind_leg = [];
for i = 1:max(indexToCluster)
    ind_leg = [ind_leg find(indexToCluster==i, 1)];
end

if legendOpt
    %legend(h, labels(lab_int))
    if ~isempty(intermediateIndex)    
        lab_int(intermediateIndex+1) = [];
        labels(intermediateIndex+1,:) = [];   
    end
    legend(h(ind_leg), labels(lab_int(ind_leg),:),'Location','northwest')
    %legend(h(ind_leg), string(indexToCluster(ind_leg)),'Location','northwest')
end


