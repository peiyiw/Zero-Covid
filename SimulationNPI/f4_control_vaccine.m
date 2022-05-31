function ydot = f4(t,theta,xdata,perC,perV,R0)

load('2019_304_cityflow.mat') % cityflow of 2019-12 months
day=load('2019_month_day.mat'); %the number of day in specific month of 2019

tSimu=1:1000;%days of simulation
n=304;  % the number of cities

%% param
% ---NPI---
perTravelBan(1:max(tSimu),1)=0;
perTravelBan(1:25,1)=1;
perControlT(1:n)=(1-perC); % intensity of social distancing

% --dynamic parameters----
pa = 0.25;% proportion of asymptomatic
rE = 1./2.9; % latent period
rA = 1./5; % recovery period for asymptomatic
rI = 1./5; % recovery period for symptomatic
rP = 1./2.3; % latent period but for pre-symptomaic
r1 = 1./3.84; % Relative infectiousness of asymptomatic infections vs. symptomatic infections
r2 = 0.15; % Relative infectiousness of pre-symptomatic infections vs. symptomatic infections
betaBasic = R0./(pa.*r1./rA+(1-pa).*(r2./rP+1./rI)); % Transmission rate for symptomatic infections
ve = 0.59; % vaccine effectiveness
betaBasicWh = betaBasic;
betaBasicWhV = betaBasic.* (1-ve);
betaWh = betaBasicWh;
betaWhV = betaBasicWhV;
betaBasicN(1:n) = betaBasic;
betaBasicNV(1:n) = betaBasic.* (1-ve);
betaN = betaBasicN;
betaNV = betaBasicNV;

% --unknown parameters----
rReportN(tSimu,1:n)=theta(2); %report rate in China
rReportN(26:end,1:n)=theta(3); %report rate in China
rReportWh(1:max(tSimu))=theta(4); %report rate in Wuhan
rReportWh(26:end)=theta(5); %report rate in Wuhan
rReportHb = theta(6); %Relative report rate in Hubei vs. other provinces
rReportN(tSimu,164:175)=rReportN(tSimu,164:175)*rReportHb;

%% data
% ----data of wuhan----
popWh=xdata.popWh; % population of Wuhan
flowWh=xdata.flowWh;% daily flow from Wuhan to other city
flowWh(26:end+max(tSimu)-max(t),:)=repmat(mean(flowWh([1:25],:)),size(flowWh,1)+max(tSimu)-max(t)-26+1,1);
travelWh=sum(flowWh,2);
time2Wh=xdata.time2Wh; % time of 2nd control of Wuhan
time2Wh(end:end+max(tSimu)-max(t),:)=1;

% ----data of others----
pop=xdata.pop; % population of other cities
time2=xdata.time2; % time of 2nd control of other cities
time2(end:end+max(tSimu)-max(t),:)=1;

%% intial state
% ----vaccination state----
HVWh(1) = perV.*popWh; % vaccinated
HV(1,1:n)=perV.*pop;

% ----initial state of wuhan----
HEWh(1) = theta(7); % latent asymptomatic 
HAWh(1) =  theta(8); % infectious asymptomatic
HPWh(1) =  theta(9); % presymptomatic
HIWh(1) = theta(10); % infectious symptomatic 
HRWh(1) = 0; % removed infectious indivial
IWh(1) = 0; % daily reported cases
HNWh(1) = popWh;%all population
HSWh(1) = HNWh(1)-HEWh(1)-HAWh(1)-HPWh(1)-HIWh(1)-HRWh(1)-HVWh(1); %susceptible

% ----initial state of others----
HE(1,1:n) = 0; %latent asymptomatic 
HA(1,1:n) = 0; % infectious asymptomatic
HP(1,1:n) = 0; %presymptomatic
HI(1,1:n) = 0; %infectious symptomatic 
HR(1,1:n) = 0; %removed infectious indivial
I(1,1:n)=0; %daily reported cases
HN(1,1:n)= pop;%all population
HS(1,1:n) = pop-HV(1,1:n); %susceptible

roundI(1,1:n+1)=round([I(1,:),IWh(1)]);

Cutoff = 1;
for i = 2:max(tSimu)
    % ---other cities---
    if sum(HA(i-1,1:n)<0)>0
        HA(i-1,find(HA(i-1,1:n)<0))=0;
    end
    if sum(HI(i-1,1:n)<0)>0
        HI(i-1,find(HI(i-1,1:n)<0))=0;
    end
    if sum(HP(i-1,1:n)<0)>0
        HP(i-1,find(HP(i-1,1:n)<0))=0;
    end
    if sum(HE(i-1,1:n)<0)>0
        HE(i-1,find(HE(i-1,1:n)<0))=0;
    end
    if sum(HS(i-1,1:n)<0)>0
        HS(i-1,find(HS(i-1,1:n)<0))=0;
    end
    if sum(HN(i-1,1:n)<0)>0
        HN(i-1,find(HN(i-1,1:n)<0))=0;
    end
    % ---wuhan---
    if HAWh(i-1)<0
        HIWh(i-1)=0;
    end
    if HIWh(i-1)<0
        HIWh(i-1)=0;
    end
    if HPWh(i-1)<0
        HPWh(i-1)=0;
    end
    if HSWh(i-1)<0
        HSWh(i-1)=0;
    end
    if HPWh(i-1)<0
        HPWh(i-1)=0;
    end
    if HNWh(i-1)<0
        HNWh(i-1)=0;
    end

    %% wuhan
    % ----beta----
    betaWh(find(time2Wh(i,:)==1))=betaBasicWh(find(time2Wh(i,:)==1)).*perControlT(find(time2Wh(i,:)==1));
    betaWhV(find(time2Wh(i,:)==1))=betaBasicWhV(find(time2Wh(i,:)==1)).*perControlT(find(time2Wh(i,:)==1));
    % ----model----
    HSWh(i)=HSWh(i-1)-(r1*betaWh.*HAWh(i-1)+betaWh.*HIWh(i-1)+r2*betaWh.*HPWh(i-1)).*HSWh(i-1)./HNWh(i-1); % suspecitble
    HVWh(i)=HVWh(i-1)-(r1*betaWhV.*HAWh(i-1)+betaWhV.*HIWh(i-1)+r2*betaWhV.*HPWh(i-1)).*HVWh(i-1)./HNWh(i-1); % suspecitble
    HEWh(i)=HEWh(i-1)+(r1*betaWh.*HAWh(i-1)+betaWh.*HIWh(i-1)+r2*betaWh.*HPWh(i-1)).*HSWh(i-1)./HNWh(i-1)+...
        +(r1*betaWhV.*HAWh(i-1)+betaWhV.*HIWh(i-1)+r2*betaWhV.*HPWh(i-1)).*HVWh(i-1)./HNWh(i-1)...
        -rE*HEWh(i-1)-travelWh(i-1)*HEWh(i-1)/HNWh(i-1); %latent
    HPWh(i)=HPWh(i-1)+(1-pa)*rE*HEWh(i-1)-rP*HPWh(i-1)-travelWh(i-1)*HPWh(i-1)/HNWh(i-1); % presymtomatic
    HAWh(i)=HAWh(i-1)+pa*rE*HEWh(i-1)-rA.*HAWh(i-1)-travelWh(i-1)*HAWh(i-1)/HNWh(i-1); % asymtomatic
    HIWh(i)=HIWh(i-1)+rP*HPWh(i-1)-rI.*HIWh(i-1)-rReportWh(i-1).*HIWh(i-1); % symtomatic
    HRWh(i)=HRWh(i-1)+rI.*HIWh(i-1)+rA.*HAWh(i-1); % removed
    IWh(i)=HIWh(i-1)*rReportWh(i-1); % daily confirmed cases
    HNWh(i)=HSWh(i)+HEWh(i)+HAWh(i)+HPWh(i)+HIWh(i)+HRWh(i)+sum(IWh); %all population
    HNWh(i)=popWh;

    %% other cities
    % ----current month and days of month to load cityflow of month----
    monthNow=mod(floor((i-1)/30)+1,12);
    if(monthNow==0)
        monthNow=12;
    end
    Fcitybalcance=cityflow(:,((1+n*(monthNow-1)):n*(monthNow)))/day.day(monthNow); % daily inter city flow of current month

    % ----beta----
    betaN(find(time2(i,:)==1))=betaBasicN(find(time2(i,:)==1)).*perControlT(find(time2(i,:)==1));
    betaNV(find(time2(i,:)==1))=betaBasicNV(find(time2(i,:)==1)).*perControlT(find(time2(i,:)==1));

    % ----Imported Infected cases----
    importE_wh=HEWh(i-1)*flowWh(i-1,:)/popWh.*perTravelBan(i-1,:); %exposed cases from wuhan to other cities
    importA_wh=HAWh(i-1)*flowWh(i-1,:)/popWh.*perTravelBan(i-1,:); %asymptomatic cases from wuhan to other cities
    importP_wh=HPWh(i-1)*flowWh(i-1,:)/popWh.*perTravelBan(i-1,:); %pre-symptomatic cases from wuhan to other cities

    % ----different class to calculate travel people proportion later----
    HEallcity=repmat(HE(i-1,:),n,1); %matrix 304*304,latent asymptomatic population of cities
    HAallcity=repmat(HA(i-1,:),n,1); %matrix 304*304,infectious asymptomatic population of cities
    HPallcity=repmat(HP(i-1,:),n,1); %matrix 304*304,presymptomatic population of cities
    HNallcity=repmat(HN(i-1,:),n,1);   %matrix 304*304,infectious symptomatic population of cities
    HSallcity=repmat(HS(i-1,:),n,1);   %matrix 304*304,suspecitble population of cities
    HVallcity=repmat(HV(i-1,:),n,1);   %matrix 304*304,vaccinated population of cities

    HEallcity(HEallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HAallcity(HAallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HPallcity(HPallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus

    % ----infected population who travel between cities----
    Fcity=Fcitybalcance*perTravelBan(i-1,:);

    inE=sum(Fcity.*HEallcity./HNallcity,2); %latent asymptomatic population who flow into the city
    outE=sum(Fcity.*HEallcity./HNallcity,1); %latent asymptomatic population who flow out of the city
    inA=sum(Fcity.*HAallcity./HNallcity,2); %infectious asymptomatic population who flow into the city
    outA=sum(Fcity.*HAallcity./HNallcity,1); %infectious asymptomatic population who flow out the city
    inP=sum(Fcity.*HPallcity./HNallcity,2); %presymptomatic population who flow into the city
    outP=sum(Fcity.*HPallcity./HNallcity,1); %presymptomatic population who flow out the city
    inS=sum(Fcity.*HSallcity./HNallcity,2);%suspectible population who flow into the city
    outS=sum(Fcity.*HSallcity./HNallcity,1);%suspectible population who flow out the city
    inV=sum(Fcity.*HVallcity./HNallcity,2);%vaccinated population who flow into the city
    outV=sum(Fcity.*HVallcity./HNallcity,1);%vaccinated population who flow out the city
    
    % --new daily cases---
    HAInfectious = HA(i-1,1:n);
    HPInfectious = HP(i-1,1:n);
    HIInfectious = HI(i-1,1:n);
    HAInfectious(HAInfectious<1)=0;
    HPInfectious(HPInfectious<1)=0;
    HIInfectious(HIInfectious<1)=0;
	tmpExporue_S = (r1*betaN.*HAInfectious+betaN.*HIInfectious+r2*betaN.*HPInfectious).*HS(i-1,:)./HN(i-1,:);
    tmpExporue_V = (r1*betaNV.*HAInfectious+betaNV.*HIInfectious+r2*betaNV.*HPInfectious).*HV(i-1,:)./HN(i-1,:);
    
    % ----model----
    HS(i,:)=HS(i-1,:)-tmpExporue_S+inS'-outS; %suspecitble population
    HV(i,:)=HV(i-1,:)-tmpExporue_V+inV'-outV; % vaccinated
    HE(i,:)=HE(i-1,:)+tmpExporue_S+tmpExporue_V+importE_wh-outE+inE'-rE*HE(i-1,:); %latent
    HP(i,:)=HP(i-1,:)+(1-pa)*rE*HE(i-1,:)-rP*HP(i-1,:)+importP_wh-outP+inP';%presymtomatic
    HA(i,:)=HA(i-1,:)+pa*rE*HE(i-1,:)-rA.*HA(i-1,:)+importA_wh-outA+inA';%asymtomatic
    HI(i,:)=HI(i-1,:)+rP*HP(i-1,:)-rI*HI(i-1,:)-rReportN(i-1,:).*HI(i-1,:); %symtomatic
    HR(i,:)=HR(i-1,:)+rI.*HI(i-1,:)+rA.*HA(i-1,:); % removed population
    I(i,:)=rReportN(i-1,:).*HI(i-1,:); %daily confirmed cases
    HN(i,:)=HS(i,:)+HE(i,:)+HA(i,:)+HP(i,:)+HI(i,:)+HR(i,:)+sum(I(:,:),1); %all population
    HN(i,:)=pop;
    roundI(i,:)=round([I(i,:),IWh(i)]);
    if i>25 && sum(roundI(i-2:i,:),"all") < Cutoff
        break
    end
end

%% result data
ydot=[I(:,:),IWh(:)];% daily reported cases