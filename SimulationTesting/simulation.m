clear
load('xdata.mat'); 

t = 1:(365*1); %days of simulation

%% get data
pop=xdata.pop;
cityName = xdata.cityName;
cityName = cellstr(char(cityName{1:length(cityName)}));

%% set parameter
% ---transmission---
theta.npi = 0.18; %npi on reduction of transmission rate  
theta.ve = 0.4; %vaccine protectoin against infection

% ---testing and tracing---
theta.k = 0.165; %the fraction of contacts that were successfully traced

% ----dynamic parameters----
theta.R0 = 10;
theta.pa = [0.194285714	0.224056604	0.189701897	0.136882129];%asymptomatic proportion for age group [ref. Beijing data]
theta.phs = [0.108233333 0.08595 0.1254 0.43025]; %hospitalization rate for symptomatic Omicron-infection in unvaccinated individuals
theta.phv = [0.03247 0.025785 0.03762 0.129075]; %hospitalization rate for symptomatic Omicron-infection in vaccinated individuals
theta.pus = [0.017233333 0.07605 0.17075 0.38945]; %ICU admission rate for symptomatic Omicron-infection in unvaccinated hospitalized patients
theta.puv = [0.00517 0.022815 0.051225 0.116835]; %ICU admission rate for symptomatic Omicron-infection in vaccinated hospitalized patients
theta.rE = 1./1.2; % 1/latent period
theta.rA = 1./3.5; % 1/recovery period for asymptomatic [ref. Jantien A Backer, et al. t, Eurosurveillance, 2022]
theta.rI = 1./3.5; % 1/recovery period for symptomatic
theta.rP = 1./(3.2 - 1.2); % 1/latent period but for pre-symptomaic [ref. Jantien A Backer, Eurosurveillance, 2022]
theta.rI2H = 1./3.5; % 1/interval between symptom onset and admission [ref. Shao, Jiasheng, et al. Vaccines， 2022]
theta.rI2U = 1./3.5; % 1/interval between symptom onset and severe illness [ref. Shao, Jiasheng, et al. Vaccines， 2022]
theta.rH = 1/6; % 1/hospital length of stay [ref. Cai, Jun, et al. Nature Medicine, 2022]
theta.rU = 1/8; % 1/icu length of stay [ref. Cai, Jun, et al. Nature Medicine, 2022]

% ---initial state---
theta.perV = 0.89; %vaccine coverage in China
theta.E0 = 50;
theta.A0 = 50;
theta.P0 = 50;
theta.I0 = 50;

%% Range
n = size(cityName,1);
intervalRange = [0,1,2,3,4,5];
lagRange = [7,14,21,28];

%% results variable
caseTotal=zeros(length(intervalRange),length(lagRange)); % total cases of nation
durAvg=zeros(length(intervalRange),length(lagRange)); % average duration of nation
roundAvg=zeros(length(intervalRange),length(lagRange)); % average round of testing of nation
testAvg =zeros(length(intervalRange),length(lagRange)); % average number of testing individuals of nation
nOutbreak=zeros(length(intervalRange),length(lagRange)); % number of outbreak cities

%% scenario simulation
for interval = intervalRange
    for lag =lagRange
       %% setting
        % ---testing and tracing---
        theta.lag = lag; %time lag
        theta.interval = interval; %testing interval
        
        % ---index---
        iLag = find(abs(lag-lagRange)<=eps);
        iInterval=find(abs(interval-intervalRange)<=eps);
        
        %% run model
        [tsCaseAge, tsCaseTotal, tsQuarAge, tsQuarTotal, tsHosAge, tsHosTotal,tsICUAge, tsICUTotal, cumHAge, cumHTotal, cumUAge, cumUTotal, bLockdown] = f4(t, theta, xdata); %simulation
        
        %% save simulation results 
        % ---save lockdown and wave matrix
        tBLockdown = array2table(bLockdown,'VariableNames',cityName);
        writetable(tBLockdown, strcat(".\simu_result\interval_",num2str(interval),"_lockdown.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");  
        
        % ---save case matrix---
        tTsCaseAge1 = array2table(reshape(tsCaseAge(:,1,:),size(tsCaseAge,1),size(tsCaseAge,3)),'VariableNames',cityName);
        tTsCaseAge2 = array2table(reshape(tsCaseAge(:,2,:),size(tsCaseAge,1),size(tsCaseAge,3)),'VariableNames',cityName);
        tTsCaseAge3 = array2table(reshape(tsCaseAge(:,3,:),size(tsCaseAge,1),size(tsCaseAge,3)),'VariableNames',cityName);
        tTsCaseAge4 = array2table(reshape(tsCaseAge(:,4,:),size(tsCaseAge,1),size(tsCaseAge,3)),'VariableNames',cityName);
        tTsCaseTotal = array2table(tsCaseTotal(:,:),'VariableNames',cityName);

        writetable(tTsCaseAge1, strcat(".\simu_result\interval_",num2str(interval),"_dailycases_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsCaseAge2, strcat(".\simu_result\interval_",num2str(interval),"_dailycases_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsCaseAge3, strcat(".\simu_result\interval_",num2str(interval),"_dailycases_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsCaseAge4, strcat(".\simu_result\interval_",num2str(interval),"_dailycases_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsCaseTotal, strcat(".\simu_result\interval_",num2str(interval),"_dailycases.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");

        % ---save qurantined matrix---
        tTsQuarAge1 = array2table(reshape(tsQuarAge(:,1,:),size(tsQuarAge,1),size(tsQuarAge,3)),'VariableNames',cityName);
        tTsQuarAge2 = array2table(reshape(tsQuarAge(:,2,:),size(tsQuarAge,1),size(tsQuarAge,3)),'VariableNames',cityName);
        tTsQuarAge3 = array2table(reshape(tsQuarAge(:,3,:),size(tsQuarAge,1),size(tsQuarAge,3)),'VariableNames',cityName);
        tTsQuarAge4 = array2table(reshape(tsQuarAge(:,4,:),size(tsQuarAge,1),size(tsQuarAge,3)),'VariableNames',cityName);
        tTsQuarTotal = array2table(tsQuarTotal(:,:),'VariableNames',cityName);

        writetable(tTsQuarAge1, strcat(".\simu_result\interval_",num2str(interval),"_quarantined_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsQuarAge2, strcat(".\simu_result\interval_",num2str(interval),"_quarantined_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsQuarAge3, strcat(".\simu_result\interval_",num2str(interval),"_quarantined_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsQuarAge4, strcat(".\simu_result\interval_",num2str(interval),"_quarantined_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsQuarTotal, strcat(".\simu_result\interval_",num2str(interval),"_quarantined.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");   

        % ---save Hos matrix---
        tTsHosAge1 = array2table(reshape(tsHosAge(:,1,:),size(tsHosAge,1),size(tsHosAge,3)),'VariableNames',cityName);
        tTsHosAge2 = array2table(reshape(tsHosAge(:,2,:),size(tsHosAge,1),size(tsHosAge,3)),'VariableNames',cityName);
        tTsHosAge3 = array2table(reshape(tsHosAge(:,3,:),size(tsHosAge,1),size(tsHosAge,3)),'VariableNames',cityName);
        tTsHosAge4 = array2table(reshape(tsHosAge(:,4,:),size(tsHosAge,1),size(tsHosAge,3)),'VariableNames',cityName);
        tTsHosTotal = array2table(tsHosTotal(:,:),'VariableNames',cityName);
        
        tCumHAge1 = array2table(reshape(cumHAge(:,1,:),size(cumHAge,1),size(cumHAge,3)),'VariableNames',cityName);
        tCumHAge2 = array2table(reshape(cumHAge(:,2,:),size(cumHAge,1),size(cumHAge,3)),'VariableNames',cityName);
        tCumHAge3 = array2table(reshape(cumHAge(:,3,:),size(cumHAge,1),size(cumHAge,3)),'VariableNames',cityName);
        tCumHAge4 = array2table(reshape(cumHAge(:,4,:),size(cumHAge,1),size(cumHAge,3)),'VariableNames',cityName);
        tCumHTotal = array2table(cumHTotal(:,:),'VariableNames',cityName);

        writetable(tTsHosAge1, strcat(".\simu_result\interval_",num2str(interval),"_hospitalization_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsHosAge2, strcat(".\simu_result\interval_",num2str(interval),"_hospitalization_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsHosAge3, strcat(".\simu_result\interval_",num2str(interval),"_hospitalization_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsHosAge4, strcat(".\simu_result\interval_",num2str(interval),"_hospitalization_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsHosTotal, strcat(".\simu_result\interval_",num2str(interval),"_hospitalization.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        
        writetable(tCumHAge1, strcat(".\simu_result\interval_",num2str(interval),"_cum_hospitalization_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumHAge2, strcat(".\simu_result\interval_",num2str(interval),"_cum_hospitalization_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumHAge3, strcat(".\simu_result\interval_",num2str(interval),"_cum_hospitalization_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumHAge4, strcat(".\simu_result\interval_",num2str(interval),"_cum_hospitalization_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
         writetable(tCumHTotal, strcat(".\simu_result\interval_",num2str(interval),"_cum_hospitalization.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        
        % ---save ICU matrix---
        tTsICUAge1 = array2table(reshape(tsICUAge(:,1,:),size(tsICUAge,1),size(tsICUAge,3)),'VariableNames',cityName);
        tTsICUAge2 = array2table(reshape(tsICUAge(:,2,:),size(tsICUAge,1),size(tsICUAge,3)),'VariableNames',cityName);
        tTsICUAge3 = array2table(reshape(tsICUAge(:,3,:),size(tsICUAge,1),size(tsICUAge,3)),'VariableNames',cityName);
        tTsICUAge4 = array2table(reshape(tsICUAge(:,4,:),size(tsICUAge,1),size(tsICUAge,3)),'VariableNames',cityName);
        tTsICUTotal = array2table(tsICUTotal(:,:),'VariableNames',cityName);
        
        tCumUAge1 = array2table(reshape(cumUAge(:,1,:),size(cumUAge,1),size(cumUAge,3)),'VariableNames',cityName);
        tCumUAge2 = array2table(reshape(cumUAge(:,2,:),size(cumUAge,1),size(cumUAge,3)),'VariableNames',cityName);
        tCumUAge3 = array2table(reshape(cumUAge(:,3,:),size(cumUAge,1),size(cumUAge,3)),'VariableNames',cityName);
        tCumUAge4 = array2table(reshape(cumUAge(:,4,:),size(cumUAge,1),size(cumUAge,3)),'VariableNames',cityName);
        tCumUTotal = array2table(cumUTotal(:,:),'VariableNames',cityName);
        
        writetable(tTsICUAge1, strcat(".\simu_result\interval_",num2str(interval),"_icu_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsICUAge2, strcat(".\simu_result\interval_",num2str(interval),"_icu_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsICUAge3, strcat(".\simu_result\interval_",num2str(interval),"_icu_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsICUAge4, strcat(".\simu_result\interval_",num2str(interval),"_icu_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tTsICUTotal, strcat(".\simu_result\interval_",num2str(interval),"_icu.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        
        writetable(tCumUAge1, strcat(".\simu_result\interval_",num2str(interval),"_cum_icu_20.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumUAge2, strcat(".\simu_result\interval_",num2str(interval),"_cum_icu_20_40.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumUAge3, strcat(".\simu_result\interval_",num2str(interval),"_cum_icu_40_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumUAge4, strcat(".\simu_result\interval_",num2str(interval),"_cum_icu_60.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tCumUTotal, strcat(".\simu_result\interval_",num2str(interval),"_cum_icu.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet");

        % ---summary---
        caseTotal(iInterval,iLag) = sum(tsCaseTotal,"all");
        
       %% duration
        [~,starind] = max((tsCaseTotal>=1),[],1); %the first day of epidemic 
        [~,endind]=max((flipud(tsCaseTotal)>0),[],1);%the last day of epidemic,sort from the end of the data
        fendind=size(tsCaseTotal,1)-endind+1;%the last day of epidemic,positive sequence
        duration=fendind-starind+1;%duration
        duration(sum(tsCaseTotal,1)==0)=0;%for city of no case,duration set to be 0;
        duration(duration<lag)=0; 
        
        tDuration = array2table(duration,'VariableNames',cityName);
        writetable(tDuration, strcat(".\simu_result\interval_",num2str(interval),"_duration.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 

        % ---summary---
        durAvg(iInterval,iLag)=mean(duration(duration~=0&duration>lag)); %average duration of nation

        %% testing
        testRound = zeros(1,n);
        testRound(duration~=0&duration>lag) =round((duration((duration~=0&duration>lag))-lag)./interval); %number of rounds of testing
        test = zeros(1,n);
        test(duration~=0&duration>lag) = round((duration(duration~=0&duration>lag)-lag)./interval).*pop(duration~=0&duration>lag); %number of test
        tRound = array2table(testRound,'VariableNames',cityName);
        tTest = array2table(test,'VariableNames',cityName);
        writetable(tRound, strcat(".\simu_result\interval_",num2str(interval),"_round.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tTest, strcat(".\simu_result\interval_",num2str(interval),"_test.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 

        % ---summary---
        roundAvg(iInterval, iLag)=mean(round((duration(duration~=0&duration>lag)-lag)./interval)); %average round of testing of nation
        testAvg(iInterval,iLag)=mean(round((duration(duration~=0&duration>lag)-lag)./interval).*pop(duration~=0&duration>lag)); %average number of testing individuals of nation
        
        %% outbreak cities
        % ---summary---
        nOutbreak(iInterval, iLag)=length(duration(duration~=0&duration>lag)); %number of outbreak cities
    end
end
%% save simulation results 
% ---total case---
tCaseTotal = array2table(caseTotal,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tCaseTotal, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("case"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---duration---
tDurAvg = array2table(durAvg,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tDurAvg, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("dur"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---testing---
tRoundAvg=array2table(roundAvg,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
tTestAvg=array2table(testAvg,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tRoundAvg, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("round"),"WriteRowNames",true,"WriteMode","overwritesheet");
writetable(tTestAvg, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("test"),"WriteRowNames",true,"WriteMode","overwritesheet");

% ---outbreak cities---
tNOutbreak = array2table(nOutbreak,'VariableNames',sprintfc('%d',lagRange),"RowNames",sprintfc('%d',intervalRange));
writetable(tNOutbreak, strcat(".\simu_result\duration_speed.xlsx"), 'Sheet',num2str("n_outbreak"),"WriteRowNames",true,"WriteMode","overwritesheet");

save simulation_result