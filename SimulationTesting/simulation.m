clear
load('xdata.mat'); 

t = 1:(365); %days of simulation
lims = [0.025,0.5,0.975];

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
theta.phs = [0.108233333 0.08595 0.1254 0.191862]; %adjusted hospitalization rate for symptomatic Omicron-infection in unvaccinated individuals
theta.phv = [0.03247 0.025785 0.03762 0.0575586]; %adjusted hospitalization rate for symptomatic Omicron-infection in vaccinated individuals
theta.pus = [0.017233333 0.07605 0.17075 0.2851525]; %adjusted ICU admission rate for symptomatic Omicron-infection in unvaccinated hospitalized patients
theta.puv = [0.00517 0.022815 0.051225 0.08554575]; %adjusted ICU admission rate for symptomatic Omicron-infection in vaccinated hospitalized patients
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
nC = size(cityName,1);
intervalRange = [1];
lagRange = [7,14,21,28];
nSample = 100;

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
        %% tmpVariable
        tsCaseAge1Res = zeros(nSample,max(t),nC);
        tsCaseAge2Res = zeros(nSample,max(t),nC);
        tsCaseAge3Res = zeros(nSample,max(t),nC);
        tsCaseAge4Res = zeros(nSample,max(t),nC);
        tsQuarAge1Res = zeros(nSample,max(t),nC);
        tsQuarAge2Res = zeros(nSample,max(t),nC);
        tsQuarAge3Res = zeros(nSample,max(t),nC);
        tsQuarAge4Res = zeros(nSample,max(t),nC);
        tsHosAge1Res = zeros(nSample,max(t),nC);
        tsHosAge2Res = zeros(nSample,max(t),nC);
        tsHosAge3Res = zeros(nSample,max(t),nC);
        tsHosAge4Res = zeros(nSample,max(t),nC);
        tsHosDisAge1Res = zeros(nSample,max(t),nC);
        tsHosDisAge2Res = zeros(nSample,max(t),nC);
        tsHosDisAge3Res = zeros(nSample,max(t),nC);
        tsHosDisAge4Res = zeros(nSample,max(t),nC);
        tsICUAge1Res = zeros(nSample,max(t),nC);
        tsICUAge2Res = zeros(nSample,max(t),nC);
        tsICUAge3Res = zeros(nSample,max(t),nC);
        tsICUAge4Res = zeros(nSample,max(t),nC);
        tsICUDisAge1Res = zeros(nSample,max(t),nC);
        tsICUDisAge2Res = zeros(nSample,max(t),nC);
        tsICUDisAge3Res = zeros(nSample,max(t),nC);
        tsICUDisAge4Res = zeros(nSample,max(t),nC);
        tsCumHosAge1Res = zeros(nSample,max(t),nC);
        tsCumHosAge2Res = zeros(nSample,max(t),nC);
        tsCumHosAge3Res = zeros(nSample,max(t),nC);
        tsCumHosAge4Res = zeros(nSample,max(t),nC);
        cumICUAge1Res = zeros(nSample,max(t),nC);
        cumICUAge2Res = zeros(nSample,max(t),nC);
        cumICUAge3Res = zeros(nSample,max(t),nC);
        cumICUAge4Res = zeros(nSample,max(t),nC);
        durationRes = zeros(nSample,1,nC);
        % 100 simulation
        for isample = 1:nSample 
            [tsCaseAge, tsQuarAge, tsHosAge, tsHosDisAge, tsICUAge, tsICUDisAge, tsCumHosAge, tsCumICUAge, duration] = f4_stochastic(t, theta, xdata); %simulation
            % age-stratified case
            tsCaseAge1Res(isample,:,:) = reshape(tsCaseAge(:,1,:),size(tsCaseAge,1),size(tsCaseAge,3));
            tsCaseAge2Res(isample,:,:) = reshape(tsCaseAge(:,2,:),size(tsCaseAge,1),size(tsCaseAge,3));
            tsCaseAge3Res(isample,:,:) = reshape(tsCaseAge(:,3,:),size(tsCaseAge,1),size(tsCaseAge,3));
            tsCaseAge4Res(isample,:,:) = reshape(tsCaseAge(:,4,:),size(tsCaseAge,1),size(tsCaseAge,3));
            

            % age-stratified quarantined individuals
            tsQuarAge1Res(isample,:,:) = reshape(tsQuarAge(:,1,:),size(tsQuarAge,1),size(tsQuarAge,3));
            tsQuarAge2Res(isample,:,:) = reshape(tsQuarAge(:,2,:),size(tsQuarAge,1),size(tsQuarAge,3));
            tsQuarAge3Res(isample,:,:) = reshape(tsQuarAge(:,3,:),size(tsQuarAge,1),size(tsQuarAge,3));
            tsQuarAge4Res(isample,:,:) = reshape(tsQuarAge(:,4,:),size(tsQuarAge,1),size(tsQuarAge,3));

            % age-stratified hospitalization (current)
            tsHosAge1Res(isample,:,:) = reshape(tsHosAge(:,1,:),size(tsHosAge,1),size(tsHosAge,3));
            tsHosAge2Res(isample,:,:) = reshape(tsHosAge(:,2,:),size(tsHosAge,1),size(tsHosAge,3));
            tsHosAge3Res(isample,:,:) = reshape(tsHosAge(:,3,:),size(tsHosAge,1),size(tsHosAge,3));
            tsHosAge4Res(isample,:,:) = reshape(tsHosAge(:,4,:),size(tsHosAge,1),size(tsHosAge,3));

            % age-stratified hospitalization considering primaty disease (current)
            tsHosDisAge1Res(isample,:,:) = reshape(tsHosDisAge(:,1,:),size(tsHosDisAge,1),size(tsHosDisAge,3));
            tsHosDisAge2Res(isample,:,:) = reshape(tsHosDisAge(:,2,:),size(tsHosDisAge,1),size(tsHosDisAge,3));
            tsHosDisAge3Res(isample,:,:) = reshape(tsHosDisAge(:,3,:),size(tsHosDisAge,1),size(tsHosDisAge,3));
            tsHosDisAge4Res(isample,:,:) = reshape(tsHosDisAge(:,4,:),size(tsHosDisAge,1),size(tsHosDisAge,3));

            % age-stratified ICU (current)
            tsICUAge1Res(isample,:,:) = reshape(tsICUAge(:,1,:),size(tsICUAge,1),size(tsICUAge,3));
            tsICUAge2Res(isample,:,:) = reshape(tsICUAge(:,2,:),size(tsICUAge,1),size(tsICUAge,3));
            tsICUAge3Res(isample,:,:) = reshape(tsICUAge(:,3,:),size(tsICUAge,1),size(tsICUAge,3));
            tsICUAge4Res(isample,:,:) = reshape(tsICUAge(:,4,:),size(tsICUAge,1),size(tsICUAge,3));

            % age-stratified ICU considering primaty disease (current)
            tsICUDisAge1Res(isample,:,:) = reshape(tsICUDisAge(:,1,:),size(tsICUDisAge,1),size(tsICUDisAge,3));
            tsICUDisAge2Res(isample,:,:) = reshape(tsICUDisAge(:,2,:),size(tsICUDisAge,1),size(tsICUDisAge,3));
            tsICUDisAge3Res(isample,:,:) = reshape(tsICUDisAge(:,3,:),size(tsICUDisAge,1),size(tsICUDisAge,3));
            tsICUDisAge4Res(isample,:,:) = reshape(tsICUDisAge(:,4,:),size(tsICUDisAge,1),size(tsICUDisAge,3));

            % age-stratified hospitalization (cumulative)
            tsCumHosAge1Res(isample,:,:) = reshape(tsCumHosAge(:,1,:),size(tsCumHosAge,1),size(tsCumHosAge,3));
            tsCumHosAge2Res(isample,:,:) = reshape(tsCumHosAge(:,2,:),size(tsCumHosAge,1),size(tsCumHosAge,3));
            tsCumHosAge3Res(isample,:,:) = reshape(tsCumHosAge(:,3,:),size(tsCumHosAge,1),size(tsCumHosAge,3));
            tsCumHosAge4Res(isample,:,:) = reshape(tsCumHosAge(:,4,:),size(tsCumHosAge,1),size(tsCumHosAge,3));

            % age-stratified ICU (cumulative)
            cumICUAge1Res(isample,:,:) = reshape(tsCumICUAge(:,1,:),size(tsCumICUAge,1),size(tsCumICUAge,3));
            cumICUAge2Res(isample,:,:) = reshape(tsCumICUAge(:,2,:),size(tsCumICUAge,1),size(tsCumICUAge,3));
            cumICUAge3Res(isample,:,:) = reshape(tsCumICUAge(:,3,:),size(tsCumICUAge,1),size(tsCumICUAge,3));
            cumICUAge4Res(isample,:,:) = reshape(tsCumICUAge(:,4,:),size(tsCumICUAge,1),size(tsCumICUAge,3));

            durationRes(isample,1,:) = duration;
        end

       %% save simulation results        
        [tsCaseTotalMean,tsCaseTotalUp,tsCaseTotalDown] = saveAgeResult(tsCaseAge1Res,tsCaseAge2Res,tsCaseAge3Res,tsCaseAge4Res,interval,lag,1,cityName,xdata);
        [tsQuarTotalMean,tsQuarTotalUp,tsQuarTotalDown] = saveAgeResult(tsQuarAge1Res,tsQuarAge2Res,tsQuarAge3Res,tsQuarAge4Res,interval,lag,2,cityName,xdata);
        [tsHosTotalMean,tsHosTotalUp,tsHosTotalDown] = saveAgeResult(tsHosAge1Res,tsHosAge2Res,tsHosAge3Res,tsHosAge4Res,interval,lag,3,cityName,xdata);
        [tsHosDisTotalMean,tsHosDisTotalUp,tsHosDisTotalDown] = saveAgeResult(tsHosDisAge1Res,tsHosDisAge2Res,tsHosDisAge3Res,tsHosDisAge4Res,interval,lag,4,cityName,xdata);
        [tsICUTotalMean,tsICUTotalUp,tsICUTotalDown] = saveAgeResult(tsICUAge1Res,tsICUAge2Res,tsICUAge3Res,tsICUAge4Res,interval,lag,5,cityName,xdata);
        [tsICUDisTotalMean,tsICUDisTotalUp,tsICUDisTotalDown] = saveAgeResult(tsICUDisAge1Res,tsICUDisAge2Res,tsICUDisAge3Res,tsICUDisAge4Res,interval,lag,6,cityName,xdata);
        [tsCumHosTotalMean,tsCumHosTotalUp,tsCumHosTotalDown] = saveAgeResult(tsCumHosAge1Res,tsCumHosAge2Res,tsCumHosAge3Res,tsCumHosAge4Res,interval,lag,7,cityName,xdata);
        [tsCumHosDisTotalMean,tsCumHosDisTotalUp,tsCumHosDisTotalDown] = saveAgeResult(tsCumHosAge1Res,tsCumHosAge2Res,tsCumHosAge3Res,tsCumHosAge4Res,interval,lag,8,cityName,xdata);
        [tsCumICUTotalMean,tsCumICUTotalUp,tsCumICUTotalDown] = saveAgeResult(cumICUAge1Res,cumICUAge2Res,cumICUAge3Res,cumICUAge4Res,interval,lag,9,cityName,xdata);
        [tsCumICUDisTotalMean,tsCumICUDisTotalUp,tsCumICUDisTotalDown] = saveAgeResult(cumICUAge1Res,cumICUAge2Res,cumICUAge3Res,cumICUAge4Res,interval,lag,10,cityName,xdata);

        % ---summary---
        caseTotal(iInterval,iLag) = sum(tsCaseTotalMean,"all");
        
        %% other
        durationCI = zeros(3,1,nC);
        testRoundMean = zeros(1,nC);
        testRoundUp = zeros(1,nC);
        testRoundDown = zeros(1,nC);
        testMean = zeros(1,nC);
        testUp = zeros(1,nC);
        testDown = zeros(1,nC);

        for iC=1:nC
            durationCI(:,1,iC) = plims(durationRes(:,:,iC),lims);
        end
        %duration
        durationMean = round(reshape(durationCI(2,:,:),1,nC));
        durationUp = round(reshape(durationCI(3,:,:),1,nC));
        durationDown = round(reshape(durationCI(1,:,:),1,nC));

        tDurationMean = array2table(durationMean,'VariableNames',cityName);
        tDurationUp = array2table(durationUp,'VariableNames',cityName);
        tDurationDown = array2table(durationDown,'VariableNames',cityName);
        writetable(tDurationMean, strcat(".\simu_result\interval_",num2str(interval),"_duration_mean.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tDurationUp, strcat(".\simu_result\interval_",num2str(interval),"_duration_up.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tDurationDown, strcat(".\simu_result\interval_",num2str(interval),"_duration_down.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 

        %test
        testRoundMean(durationMean~=0&durationMean>lag) =round((durationMean((durationMean~=0&durationMean>lag))-lag)./interval); %number of rounds of testing
        testRoundUp(durationUp~=0&durationUp>lag) =round((durationUp((durationUp~=0&durationUp>lag))-lag)./interval); %number of rounds of testing
        testRoundDown(durationDown~=0&durationDown>lag) =round((durationDown((durationDown~=0&durationDown>lag))-lag)./interval); %number of rounds of testing
    
        testMean(durationMean~=0&durationMean>lag) = round((durationMean(durationMean~=0&durationMean>lag)-lag)./interval).*pop(durationMean~=0&durationMean>lag); %number of test
        testUp(durationUp~=0&durationUp>lag) = round((durationUp(durationUp~=0&durationUp>lag)-lag)./interval).*pop(durationUp~=0&durationUp>lag); %number of test
        testDown(durationDown~=0&durationDown>lag) = round((durationDown(durationDown~=0&durationDown>lag)-lag)./interval).*pop(durationDown~=0&durationDown>lag); %number of test
        
        tRoundMean = array2table(testRoundMean,'VariableNames',cityName);
        tRoundUp = array2table(testRoundUp,'VariableNames',cityName);
        tRoundDown = array2table(testRoundDown,'VariableNames',cityName);
        writetable(tRoundMean, strcat(".\simu_result\interval_",num2str(interval),"_round_mean.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tRoundUp, strcat(".\simu_result\interval_",num2str(interval),"_round_up.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tRoundDown, strcat(".\simu_result\interval_",num2str(interval),"_round_down.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 

        tTestMean = array2table(testMean,'VariableNames',cityName);
        tTestUp = array2table(testUp,'VariableNames',cityName);
        tTestDown = array2table(testDown,'VariableNames',cityName);
        writetable(tTestMean, strcat(".\simu_result\interval_",num2str(interval),"_test_mean.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tTestUp, strcat(".\simu_result\interval_",num2str(interval),"_test_up.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 
        writetable(tTestDown, strcat(".\simu_result\interval_",num2str(interval),"_test_down.xlsx"), ...
            'Sheet',num2str(lag),"WriteMode","overwritesheet"); 


        % ---summary---
        durAvg(iInterval,iLag)=mean(durationMean(durationMean~=0&durationMean>lag)); %average duration of nation
        roundAvg(iInterval, iLag)=mean(round((durationMean(durationMean~=0&durationMean>lag)-lag)./interval)); %average round of testing of nation
        testAvg(iInterval,iLag)=mean(round((durationMean(durationMean~=0&durationMean>lag)-lag)./interval).*pop(durationMean~=0&durationMean>lag)); %average number of testing individuals of nation
        
        %% outbreak cities
        % ---summary---
        nOutbreak(iInterval, iLag)=length(durationMean(durationMean~=0&durationMean>lag)); %number of outbreak cities
        
        save(strcat("simulation_result_interval_",num2str(interval),"_lag_",num2str(lag)))
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