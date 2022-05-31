clear
load('xdata.mat'); 
load('2019_305_cityflow.mat'); 
load('city_label.mat'); 
day=load('2019_month_day.mat'); %the number of day in specific month of 2019

%% parameter
% ---number of cities---
cityName = City.name;
theta.n = 305; % total number of cities
theta.iCenter=1; % Wuhan
cityName = cellstr(str2mat(cityName{1:theta.n}));

% ---transmission rate---
ve = 0.4;
npi = 0.18;
theta.R0 = 10; %basic reproduction rate
theta.betan = (1-ve)*(1-npi);

% ---flow---
theta.cityflow=cityflow;
theta.day=day.day;

% ---testing and tracing---
theta.k = 0.165;
theta.lag0 = 3; %time lag in outbreak city
theta.lag1 = 3; %time lag in other cities
theta.break0 = 0; %break interval in outbreak city
theta.break1 = 0; %break interval in other cities

% ---initial state---
theta.E0 = 50;
theta.A0 = 50;
theta.P0 = 50;
theta.I0 = 50;

% ---days of simulation---
t = 1:1000;

%% Range
intervalRange = [1:1:5];
lagRange = [3:1:30];

%% results variable
res=cell(1,1); % save primary data of model
% ---total case---
totalCaseOther=zeros(length(intervalRange),length(lagRange)); % total cases of other cities
% ---duration---
avgDurOther=zeros(length(intervalRange),length(lagRange)); % average duration of other cities
% ---testing---
avgRoundOther=zeros(length(intervalRange),length(lagRange)); % average round of testing of other cities
avgNumTestOther = zeros(length(intervalRange),length(lagRange)); % average number of testing individuals of other cities
% ---outbreak cities---
nOutbreak=zeros(length(intervalRange),length(lagRange)); % number of outbreak cities (except for epicenter)
% ---peak size---
avgPeakSizeOther=zeros(length(intervalRange),length(lagRange)); % average max peak size of other cities

% scenario simulation
for interval = intervalRange
    for lag =lagRange
       %% setting
        % ---testing and tracing---
        theta.lag0 = lag;
        theta.lag1 = lag; %time lag in other cities
        theta.testing0 = 1/interval; %testing interval in outbreak city
        theta.testing1 = 1/interval; %testing interval in other cities
    
        % ---index---
        iLag = find(abs(lag-lagRange)<=eps);
        iInterval=find(abs(interval-intervalRange)<=eps);
        
       %% run model
        ydot = f4(t, theta, xdata); %daily reported cases

       %% case
        yint=round(ydot);  % The daily number of cases of each city is rounded as integer
        sumYint=sum(yint,2); % the daily number of cases of nation
        
        % ---save case matrix---
        res{1,1}{iLag, iInterval}=yint;
        totalCaseNation(iInterval, iLag)=sum(sumYint);
        totalCaseEpicenter(iInterval, iLag) = sum(yint(:,theta.iCenter));
        totalCaseOther(iInterval, iLag) = totalCaseNation(iInterval, iLag)-totalCaseEpicenter(iInterval, iLag);

        tDailyCases = array2table(yint,'VariableNames',cityName);
        writetable(tDailyCases, strcat(".\simu_result\lag_duration_interval_",num2str(interval),"_dailycases.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");      
        
       %% duration
        % ---average duration of other cities---
        [~,starind] = max((yint>=1),[],1); %the first day of epidemic 
        [~,endind]=max((flipud(yint)>0),[],1);%the last day of epidemic,sort from the end of the data
        fendind=size(yint,1)-endind+1;%the last day of epidemic,positive sequence
        duration=fendind-starind+1;%duration
        duration(find(sum(yint,1)==0))=0;%for city of no case,duration set to be 0;
        durationOther=duration;
        durationOther(theta.iCenter)=[];
        avgDurOther(iInterval, iLag)=mean(durationOther);%average duration of cities 
        
        %% testing
        % ---average round of testing of other cities---
        iother0 = find(durationOther==0);
        iother1 = find(durationOther~=0);
        avgRoundOther(iInterval, iLag)=mean([durationOther(iother0),round((durationOther(iother1)-lag)./interval)]);
        
        % average number of testing individuals of other cities
        avgNumTestOther(iInterval, iLag)=mean([durationOther(iother0).*xdata.pop(iother0),round((durationOther(iother1)-lag)./interval).*xdata.pop(iother1)]);
        
        %% outbreak cities
         % ---number of outbreak cities(except for epicenter)---
        nOutbreak(iInterval, iLag)=sum(durationOther~=0);

        %% peak size
        % ---average peak size of other cities---
        MaxinciOther=Maxinci;
        MaxinciOther(theta.iCenter)=[];
        avgPeakSizeOther(iInterval,iLag)=mean(MaxinciOther);%average max peak size of cities
    end
end
%% save simulation results 
% ---total case---
tTotalCaseOther = array2table(totalCaseOther,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tTotalCaseOther, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("total_case_other"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---duration---
tAvgDurOther = array2table(avgDurOther,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tAvgDurOther, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("avg_dur_other"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---testing---
tAvgRoundOther=array2table(avgRoundOther,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
tAvgNumTestOther=array2table(avgNumTestOther,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tAvgRoundOther, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("avg_round_other"),"WriteRowNames",true,"WriteMode","overwritesheet");
writetable(tAvgNumTestOther, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("avg_num_test_other"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---outbreak cities---
tNOutbreak = array2table(nOutbreak,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tNOutbreak, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("n_outbreak"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---peak size---
tAvgPeakSizeOther = array2table(avgPeakSizeOther,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tAvgPeakSizeOther, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("avg_peak_size_other"),"WriteRowNames",true,"WriteMode","overwritesheet");

save simulation_result