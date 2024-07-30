format long
dat_raw = readtable('EMPIREanon.csv');
load('kernelList.mat');
load('funcNum082922.mat');
%%
timestamp = dat_raw.timestamp;
nodes = dat_raw.node_name;
job_id = dat_raw.job_id;
ranks = dat_raw.rank;
name = dat_raw.name;
ckc = dat_raw.current_kernel_count;
ckt = dat_raw.current_kernel_time;
tkc = dat_raw.total_kernel_count;

%save original jobs
jobs = job_id;

%% Eliminate NaN values of jobID

job_id = jobs;

ind_Nan = find(isnan(job_id));

timestamp(ind_Nan) = [];
nodes(ind_Nan) = [];
job_id(ind_Nan) = [];
ranks(ind_Nan) = [];
name(ind_Nan) = [];
ckc(ind_Nan) = [];
ckt(ind_Nan) = [];

%{
funcNum = zeros(length(name),1)-1;
for j = 1:length(name)
    temp = sscanf(name{j},'%13c%d');
    funcNum(j) = temp(end);
end
%}
%{
D = regexp(name,'\d+','match');
funcNum2 = str2double(vertcat(D{:}));
%}
jobID = unique(job_id);
%functionNames = unique(name);

% jobID: 12497415,12497416,12497440 don't exist in file.
%jobID_bad = [12421212,12421220,12423332,12423334,12497415,12497416,12497440,12562163];
jobID_bad = [12421212,12421220,12423332,12423333,12423334,12562163];
jobID_good = setdiff(jobID, jobID_bad);

%{
ind_max = find(ckc == max(ckc)); %find(tkc == max(tkc));

max_node = nodes(ind_max);
max_rank = ranks(ind_max);
max_name = name(ind_max);
%}


%%

titleString = "all"; jobNum = jobID; 
%titleString = "normal"; jobNum = jobID_good; jobNum = jobID_good([1,3,5]); 
%titleString = "bad"; jobNum = jobID_bad; jobNum = jobID_bad([2,3,5]);

%jobNum = jobID_bad(1:3); %jobID;
%jobNum = [jobNum(17); jobNum(32:end)];

%functionName = 'kernel_funct_268';
%kernelList=276;

Flag = 0;
for ii = 1:numel(kernelList)
    functionName = strcat('kernel_function_',num2str(kernelList(ii)));
    %rankNum = unique(ranks); %rankNum = 6;
    %nodeNum = 4;


    % ind_func = find(contains(name, functionName)); %% wrong!!
    ind_func = find(funcNum == kernelList(ii));
    Stimestamp = timestamp(ind_func);
    Snodes = nodes(ind_func);
    Sranks = ranks(ind_func);
    Sname = name(ind_func);
    Sckc = ckc(ind_func);
    Sckt = ckt(ind_func);
    Stkc = tkc(ind_func);
    Sjob_id = job_id(ind_func);


    %test{ii} = unique(Sname);
    %{
    dat = cell(1,1);
    dat_s = cell(1,1);
    dat_j = cell(1,1);
    tim = cell(1,1);
    tim_s = cell(1,1);
    tim_j = cell(1,1);
    %}
    dat = cell(1,1);
    tim = cell(1,1);
    %
    figure
    t = tiledlayout(2,1);
    ax1 = nexttile; title('raw data')
    ax2 = nexttile; title('average in 10s window')

    figure
    t2 = tiledlayout(2,1);
    ax21 = nexttile; title('raw data')
    ax22 = nexttile; title('average in 10s window')
    %}

    rankList = unique(Sranks);
    for i = 1:length(jobNum)

        Ltimestamp = Stimestamp;
        Lnodes = Snodes;
        Lranks = Sranks;
        Lname = Sname;
        Lckc = Sckc;
        Lckt = Sckt;
        Ltkc = Stkc;

        ind = find(Sjob_id == jobNum(i)); 
        Lranks = Lranks(ind);
        Lckt = Lckt(ind);
        Ltimestamp = Ltimestamp(ind);
        Lnodes = Lnodes(ind);
        Lckc = Lckc(ind);
        Lname = Lname(ind);



        %{
        temp_timestamp = cell(length(rankList),1);
        temp_ckt = cell(length(rankList),1);
        for j = 1:length(rankList)
            ind = find(Lranks == rankList(j));
            temp_time = Ltimestamp(ind);
            temp_timestamp{j} = temp_time - min(temp_time); % offset timestamp
            temp_ckt{j} = Lckt(ind);
            [temp_timestamp{j}, I] = sort(temp_timestamp{j});
            temp_ckt{j} = temp_ckt{j}(I);
        end
        %}


        % Get rid of duplicate data tuple
        [uniqueTime k l] = unique(Ltimestamp,'first');
        indexToDupes = find(not(ismember(1:numel(Ltimestamp),k)));
        if ~isempty(indexToDupes)
            Lranks(indexToDupes) = [];
            Lckt(indexToDupes) = [];
            Ltimestamp(indexToDupes) = [];
            Lnodes(indexToDupes) = [];
            Lckc(indexToDupes) = [];
            Lname(indexToDupes) = [];
        end



        % Adjust timestamp to overlay
        timeOffset = min(Ltimestamp);
        Ltimestamp = Ltimestamp - timeOffset;

        % sort timestamp in ascending order and ckt accordingly
        [Ltimestamp, I] = sort(Ltimestamp);
        Lckt = Lckt(I);
        Lckc = Lckc(I);
        %D_Lckc = Lckc(2:end) - Lckc(1:end-1);
        %tim{1}{i} = Ltimestamp;
        %dat{1}{i} = Lckt;

        %
        hold([ax1 ax2],'on')
        plot(ax1,Ltimestamp,Lckt)
        %scatter(ax2,Ltimestamp,Lckt,'filled')
        %}
        % Uniform moving window for time
        windowWidth = 5;
        %
        tim_base = 0;
        next = 0;
        j = 1;
        k = 1;
        testTime = {};
        testCkt = [];
        testCkc = [];
        if ~isempty(Ltimestamp)
            while tim_base <= Ltimestamp(end)
            %while next <= 300 
                next = tim_base + windowWidth;
                temp = [];
                testCkt(k) = 0;
                testCkc(k) = 0;
                while j <= length(Ltimestamp) && Ltimestamp(j)< next
                    temp = [temp, Ltimestamp(j)];
                    testCkt(k) = testCkt(k) + Lckt(j);
                    testCkc(k) = testCkc(k) + Lckc(j);
                    %testCkt(k) = testCkt(k) + D_Lckc(j);
                    j = j + 1;
                end
                testTime{k} = temp;
                k = k + 1;
                tim_base = next;
            end
            testCkt = testCkt/length(rankList);
            testCkc = testCkc/length(rankList);
            testTimestamp = 0:windowWidth:next-windowWidth;
            %
            if length(testTimestamp)<2
                Flag = 1;
                break;
            end
                
            tim{1}{i} = testTimestamp;
            dat{1}{i} = testCkt;
            
            plot(ax2,testTimestamp,testCkt)
            %
            hold([ax21 ax22],'on')
            %;
            plot(ax21,Ltimestamp,Lckc)
            plot(ax22,testTimestamp,testCkc)
            %scatter(ax22,testTimestamp,testCkt,'filled')
            %}
        end

    end
    %hold off

    %{
    % For the "1st" tiled layout
    xlabel(t,'Time (s)')
    ylabel(t,'Current kernel time')
    t.TileSpacing = 'compact';
    t.Padding = 'compact';
    % add legend
    legendCell = strcat('jobID:',string(num2cell(jobNum)));
    lg  = legend(ax2,legendCell,'Orientation','Horizontal','NumColumns',10); 
    lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
    % add title
    [upT, lowT] = title(t, sprintf('"%s" jobs',titleString), sprintf('%s', functionName));
    upT.FontSize = 16;
    set(lowT ,'Interpreter','none');


    % For the "2nd" tiled layout
    xlabel(t2,'Time (s)')
    ylabel(t2,'Current kernel count')
    %ylabel(t2,'\Delta Current kernel count')
    t2.TileSpacing = 'compact';
    t2.Padding = 'compact';
    % -- add legend
    legendCell = strcat('jobID:',string(num2cell(jobNum)));
    lg  = legend(ax22,legendCell,'Orientation','Horizontal','NumColumns',10); 
    lg.Layout.Tile = 'North'; % <-- Legend placement with tiled layout
    % -- add title
    [upT, lowT] = title(t2, sprintf('"%s" jobs',titleString), sprintf('%s', functionName));
    upT.FontSize = 16;
    set(lowT ,'Interpreter','none');
    %}
    %}

    if Flag == 1
        Flag = 0 ;
        kernelList{ii} = [];
        continue;
    end

    labels = cell(length(jobID),1);

    jobID_bad = [12421212,12421220,12423332,12423333,12423334,12562163];
    jobID_bad_ind = ismember(jobID_bad, jobID);
    for i = 1:length(jobID)
        labels{i} = "normal";
    end
    labels = string(labels);

    labels_num = zeros(length(jobID),1);
    for i = 1:length(jobID_bad)
        if jobID_bad_ind(i)==1
            ind_bad = find(jobID == jobID_bad(i));
            labels(ind_bad) = "bad";
            labels_num(ind_bad) = 1;
        end
    end
    %save('EMPIRE_dat_window_1s_first5Min_kf_7.mat', 'labels', 'dat', 'tim', 'jobID');
    %save('EMPIRE_dat_window_10s_kf_306.mat', 'labels', 'dat', 'tim', 'jobID');
    %save(fullfile(['Data_EMPIRE\slidingW_5s\EMPIRE_dat_window_5s_' functionName '_' datestr(now,'mmddyy')]), 'labels', 'labels_num', 'dat', 'tim', 'jobID');
    
end

%kernelList = kernelList(~cellfun('isempty',kernelList));
%kernelList = cell2mat(kernelList);
%save(fullfile(['new_kernelList_' datestr(now,'mmddyy')]), 'kernelList');
