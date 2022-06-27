clear
load ('result.mat');
load('data.mat'); 

%% parameter
R0=3.2; %***************************

n=304;  % the number of cities
t=1:66; % period 1/1-3/6
tSimu=1:1000;% days of simulation
xdata=data.xdata;
cityName=data.xdata.cityName;
cityName = cellstr(str2mat(cityName{1:end},"武汉市"));

theta=results2.mean; % results of f1.m

%% Range
controlRange = [0:0.05:0.95];
travelBanRange  = [0:0.05:1];

%% results variable
res=cell(1,1); % save primary data of model
totalCaseNation=zeros(length(travelBanRange),length(controlRange)); % total cases of nation
avgDurOther=zeros(length(travelBanRange),length(controlRange)); % average duration of other cities
nOutbreak=zeros(length(travelBanRange),length(controlRange)); % number of outbreak cities (except for epicenter)

%% model
for perC=controlRange % control intensity of social distancing 
    for perTravelBan=travelBanRange % control intensity of inter-city flow
        %% setting
        % ---index---
        iflow=find(abs(perTravelBan-travelBanRange)<=eps);
        iC=find(abs(perC-controlRange)<=eps);
        
       %% run model
        ydot = f4(t,theta,xdata,perTravelBan,perC,R0); %daily reported cases

       %% case
        % ---case---
        yint=round(ydot);  % The daily number of cases of each city is rounded as integer
        sumYint=sum(yint,2); % the daily number of cases of nation
        
        % ---save case matrix---
        res{1,1}{iflow,iC}=yint; % the daily number of cases of each city
        totalCaseNation(iflow,iC)=sum(sumYint); % the total number of cases of nation
        
%         tDailyCases = array2table(yint,'VariableNames',cityName);
%         writetable(tDailyCases, strcat(".\simu_result\Control_",num2str(perC),"_dailycases.xlsx"), ...
%             'Sheet',num2str(perTravelBan),"WriteMode","overwritesheet"); 
        
        %% duration
        % ---average duration of other cities---
        [~,starind] = max((yint>=1),[],1); %the first day of epidemic 
        [~,endind]=max((flipud(yint)>0),[],1);%the last day of epidemic,sort from the end of the data
        fendind=size(yint,1)-endind+1; % the last day of epidemic,positive sequence
        duration=fendind-starind+1;%duration
        duration(find(sum(yint,1)==0))=0;%for city of no case,duration set to be 0;
        avgDurOther(iflow,iC)=mean(duration(1:n));
        
        %% outbreak cities
         % ---number of outbreak cities(except for epicenter)---
        nOutbreak(iflow, iC)=sum(duration(1:n)~=0);
   
    end
end
%% save simulation results 
tDurOther = array2table(avgDurOther,'VariableNames',sprintfc('%d',controlRange),"RowNames",sprintfc('%d',travelBanRange));
tNOutbreak = array2table(nOutbreak,'VariableNames',sprintfc('%d',controlRange),"RowNames",sprintfc('%d',travelBanRange));
writetable(tDurOther, strcat('.\simu_result\average_duration_R0_',num2str(R0),'.xlsx'),"Sheet","avg_dur_other","WriteRowNames",true,"WriteMode","overwritesheet");
writetable(tNOutbreak, strcat('.\simu_result\average_duration_R0_',num2str(R0),'.xlsx'),"Sheet","n_outbreak","WriteRowNames",true,"WriteMode","overwritesheet");
save result_sinulation;