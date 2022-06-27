clear
load('xdata.mat'); 
load('2019_366_cityflow_avg.mat'); 
load('2019_366_city_label.mat'); 

%% parameter
% ---number of cities---
cityName = City.name;
theta.n = 366; % total number of cities
theta.iCenter=168; % Wuhan
cityName = cellstr(str2mat(cityName{1:theta.n}));

% ---vaccination---
ve=0.4;
theta.ve = ve; % vaccine protectoin against infection

% ---transmission rate---
npi = 0.18; %npi on reduction of transmission rate  
theta.betan = 1-npi;

% ---flow---
theta.cityflow=cityflow;

% ---testing and tracing---
theta.k = 0.165; %the fraction of contacts that were successfully traced
lag=21; % ***response lag***
theta.lag0 = lag; %time lag in outbreak city
theta.lag1 = lag; %time lag in other cities
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
coverageRange = [0:0.05:1];
intervalRange  = [1:1:4];

%% results variable
res=cell(1,1); % save primary data of model
% ---total case---
totalCaseNation=zeros(length(intervalRange),length(coverageRange)); % total cases of nation
% ---duration---
avgDurNation=zeros(length(intervalRange),length(coverageRange)); % average duration of nation
% ---testing---
avgRoundNation=zeros(length(intervalRange),length(coverageRange)); % average round of testing of nation
% ---outbreak cities---
nOutbreak=zeros(length(intervalRange),length(coverageRange)); % number of outbreak cities (except for epicenter)

for R0 = [10,15]
   %% scenario simulation
    for perV = coverageRange
        for interval = intervalRange
            %% setting
            theta.R0 = R0; %basic reproduction rate
            % ---vaccination---
            theta.perV = perV; %vaccine coverage
            % ---testing and tracing---
            if interval == 0
                theta.testing0 = 0; %testing interval in outbreak city
                theta.testing1 = 0; %testing interval in other cities
            else
                theta.testing0 = 1/interval; %testing interval in outbreak city
                theta.testing1 = 1/interval; %testing interval in other cities
            end
            % ---index---
            iInterval = find(abs(interval-intervalRange)<=eps);
            iCoverage=find(abs(perV-coverageRange)<=eps);

            %% run model
            ydot = f4(t, theta, xdata); %daily reported cases

            %% case
            yint=round(ydot.case);  % The daily number of cases of each city is rounded as integer
            sumYint=sum(yint,2); % the daily number of cases of nation

            % ---save case matrix---
            totalCaseNation(iInterval, iCoverage)=sum(sumYint);

            tDailyCases = array2table(yint,'VariableNames',cityName);
            writetable(tDailyCases, strcat(".\simu_result\R0_",num2str(R0),"_Coverage_",num2str(perV),"_dailycases.xlsx"), ...
                'Sheet',num2str(interval),"WriteMode","overwritesheet"); 

            %% lockdown
            bLockdown = ydot.bLockdown;
            tBLockdown = array2table(bLockdown,'VariableNames',cityName);
            writetable(tBLockdown, strcat(".\simu_result\lag_duration_interval_",num2str(interval),"_lockdown.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");      

            %% duration
            % ---average duration of nation---
            [~,starind] = max((yint>=1),[],1); %the first day of epidemic 
            [~,endind]=max((flipud(yint)>0),[],1);%the last day of epidemic,sort from the end of the data
            fendind=size(yint,1)-endind+1;%the last day of epidemic,positive sequence
            duration=fendind-starind+1;%duration
            duration(find(sum(yint,1)==0))=0;%for city of no case,duration set to be 0;
            duration(find(duration==0))=[];
            avgDurNation(iInterval, iCoverage)=mean(duration);

            %% testing
            % ---average round of testing of nation---
            avgRoundNation(iInterval, iCoverage)=mean(round((duration-lag)./interval));
            
            %% outbreak cities
            % ---number of outbreak cities(except for epicenter)---
            nOutbreak(iInterval, iCoverage)=length(duration);

        end
    end
    %% save simulation results 
    % ---total case---
    tTotalCaseNation = array2table(totalCaseNation,'VariableNames',sprintfc('%d',coverageRange),"RowNames",sprintfc('%d',intervalRange));
    writetable(tTotalCaseNation, strcat(".\simu_result\ve_",num2str(ve),"_R0_",num2str(R0),"_interval_duration.xlsx"), 'Sheet',num2str("total_case_nation"),"WriteRowNames",true,"WriteMode","overwritesheet");
    % ---duration---
    tAvgDurNation = array2table(avgDurNation,'VariableNames',sprintfc('%d',coverageRange),"RowNames",sprintfc('%d',intervalRange));
    writetable(tAvgDurNation, strcat(".\simu_result\ve_",num2str(ve),"_R0_",num2str(R0),"_interval_duration.xlsx"), 'Sheet',num2str("avg_dur_nation"),"WriteRowNames",true,"WriteMode","overwritesheet");
    % ---testing---
    tAvgRoundNation=array2table(avgRoundNation,'VariableNames',sprintfc('%d',coverageRange),"RowNames",sprintfc('%d',intervalRange));
    writetable(tAvgRoundNation, strcat(".\simu_result\ve_",num2str(ve),"_R0_",num2str(R0),"_interval_duration.xlsx"), 'Sheet',num2str("avg_round_nation"),"WriteRowNames",true,"WriteMode","overwritesheet");
    % ---outbreak cities---
    tNOutbreak = array2table(nOutbreak,'VariableNames',sprintfc('%d',coverageRange),"RowNames",sprintfc('%d',intervalRange));
    writetable(tNOutbreak, strcat(".\simu_result\ve_",num2str(ve),"_R0_",num2str(R0),"_interval_duration.xlsx"), 'Sheet',num2str("n_outbreak"),"WriteRowNames",true,"WriteMode","overwritesheet");

end
save simulation_result