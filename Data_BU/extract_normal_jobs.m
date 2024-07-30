load('BU_data_trimmed_032822_w_abs2min.mat')
dat_raw = readtable('\BUdata_trimmed_normal\metadata.csv');
ind_normal = find('None'==anomaly);

nodes = dat_raw.nodes(ind_normal);


%%
dat_normal = cell(1,2);
dat_normal{1} = dat_trimmed{1}(ind_normal);
dat_normal{2} = dat_trimmed{2}(ind_normal);
% dat_normal = cell(1,2);
% dat_normal{1} = dat{1}(ind_normal);
% dat_normal{2} = dat{2}(ind_normal);

labels_normal = labels(ind_normal);
jobID_normal = jobID(ind_normal);
%%
classes = unique(labels_normal);
labels_normal_num = zeros(length(labels_normal),1);
for ii = 1:length(classes)
    ind = find(labels_normal==classes(ii));
    labels_normal_num(ind) = ii;
end

%{
label_text = ["ExaMiniMD" "LAMMPS" "sw4lite" "sw4" "SWFFT" "HACC"]';
for ii = 1:length(label_text)
    ind = find(labels_normal==label_text(ii));
    labels_normal_num(ind) = ii;
end
%}

dat = dat_normal;
labels = labels_normal_num;

%------ Select 170 samples(jobs) for training set
%{
r = randsample(693, 170);
training = unique(r);
testing = setdiff(1:693, training);

job_training = jobID_normal(training);
job_testing = jobID_normal(testing);
label_training = labels_normal(training);
label_testing = labels_normal(testing);

appname = [label_training; label_testing];
jobID = [job_training; job_testing];

timeseries_id = [1:length(appname)]';
label_table = table(timeseries_id, appname, jobID);
writetable(label_table, 'labels.csv','Delimiter',',','QuoteStrings',1)

% Save training data
for i = 1:length(training)
    activeMem = dat{1}{training(i)};
    cpuUtil = dat{2}{training(i)};
    metric_table = table(activeMem, cpuUtil);
    %save(fullfile([pwd '\data\' num2str(i)]), 'activeMem', 'cpuUtil',
    %'-v7.3');
    writetable(metric_table, fullfile([pwd '\data\' num2str(i) '.csv']),'Delimiter',',','QuoteStrings',1)
end

% Save testing data
for j = 1:length(testing)
    activeMem = dat{1}{testing(j)};
    cpuUtil = dat{2}{testing(j)};
    metric_table = table(activeMem, cpuUtil);
    writetable(metric_table, fullfile([pwd '\data\' num2str(j+length(training)) '.csv']),'Delimiter',',','QuoteStrings',1)
end
    
%}

%------ delete short jobs-------
%{
dat_delete = sort([268,530,675],'descend');
for i = 1:numel(dat_delete)
    dat{1}{dat_delete(i)} = [];
    dat{2}{dat_delete(i)} = [];
end
dat{1} = dat{1}(~cellfun(@isempty, dat{1}));
dat{2} = dat{2}(~cellfun(@isempty, dat{2}));

labels(dat_delete) = [];
%}
%save(fullfile([pwd '\dat_All_' datestr(now,'mmddyy')]), 'dat','classes','labels', '-v7.3');