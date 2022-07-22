function ydot = f4(t,theta,xdata)

t = 1:66; %time period 2020/1.1-3.6
T = max(t);
n = 304;  %number of city 

%% param
% ----dynamic parameters----
R0 = 3.2; 
pa = 0.25; % proportion of asymptomatic
rE = 1./2.9; % latent period^(-1)
rA = 1./5; % recovery rate for asymptomatic
rI = 1./5; % recovery rate for symptomatic
rP = 1./2.3; % latent period but for pre-symptomaic^(-1)
r1 = 1./3.84; % relative infectiousness of asymptomatic infections vs. symptomatic infections
r2 = 0.15; % relative infectiousness of pre-symptomatic infections vs. symptomatic infections
betaBasic = R0./(pa.*r1./rA+(1-pa).*(r2./rP+1./rI)); % basic transmission rate for symptomatic infections
betaBasicWh = betaBasic; % basic transmission rate in Wuhan
betaWh = betaBasicWh;
betaBasicN(1:n) = betaBasic; % basic transmission rate in other cities
betaN = betaBasicN;

% ----unknown parameters----
c2(1:n)=(1-theta(1)); % control
rReportN(t,1:n)=theta(2); % report rate in other citeis before 1.25
rReportN(26:end,1:n)=theta(3); % report rate in other citeis after 1.25
rReportWh(1:T)=theta(4); % report rate in Wuhan before 1.25
rReportWh(26:end)=theta(5); %report rate in Wuhan after 1.25
rReportHb = theta(6); % relative report rate in Hubei vs. other provinces
rReportN(t,164:175)=rReportN(t,164:175)*rReportHb;

%% data
% ----data of wuhan----
popWh=xdata.popWh; % population of Wuhan
flowWh=xdata.flowWh; % daily flow from Wuhan to other city
travelWh=sum(flowWh,2);
time2Wh=xdata.time2Wh; % time of 2nd control of Wuhan

% ----data of others----
pop=xdata.pop; % population of other cities
time2=xdata.time2; % time of 2nd control of other cities

% ----inter-city flow----
flowJan=xdata.flowJan/31;   % intra-city travel in JAN.
flowFeb=xdata.flowFeb/28;  % intra-city travel in FEB.
flowMar=xdata.flowMar/31;  % intra-city travel in MAR.
flow = flowJan;

%% intial state
% ----initial state of wuhan----
HEWh(1) = theta(7); % latent asymptomatic 
HAWh(1) =  theta(8); % infectious asymptomatic
HPWh(1) =  theta(9); % presymptomatic
HIWh(1) = theta(10); % infectious symptomatic 
HRWh(1) = 0; % removed infectious indivial
IWh(1) = 0; % daily reported cases
HNWh(1) = popWh; % all population
HSWh(1) = HNWh(1)-HEWh(1)-HAWh(1)-HPWh(1)-HIWh(1)-HRWh(1); %susceptible

% ----initial state of others----
HE(1,1:n) = 0; % latent asymptomatic 
HA(1,1:n) = 0; % infectious asymptomatic
HP(1,1:n) = 0; % presymptomatic
HI(1,1:n) = 0; %i nfectious symptomatic 
HR(1,1:n) = 0; % removed infectious indivial
I(1,1:n)=0; % daily reported cases
HN(1,1:n)= pop;% all population
HS(1,1:n) = pop; % susceptible

%%
for i = 2:66 %1.1-1.23
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
    betaWh(find(time2Wh(i,:)==1))=betaBasicWh(find(time2Wh(i,:)==1)).*c2(find(time2Wh(i,:)==1));

    % ----model----
    HSWh(i)=HSWh(i-1)-(r1*betaWh.*HAWh(i-1)+betaWh.*HIWh(i-1)+r2*betaWh.*HPWh(i-1)).*HSWh(i-1)./HNWh(i-1); % suspecitble
    HEWh(i)=HEWh(i-1)+(r1*betaWh.*HAWh(i-1)+betaWh.*HIWh(i-1)+r2*betaWh.*HPWh(i-1)).*HSWh(i-1)./HNWh(i-1)-rE*HEWh(i-1)-travelWh(i-1)*HEWh(i-1)/HNWh(i-1); %latent
    HPWh(i)=HPWh(i-1)+(1-pa)*rE*HEWh(i-1)-rP*HPWh(i-1)-travelWh(i-1)*HPWh(i-1)/HNWh(i-1); % presymtomatic
    HAWh(i)=HAWh(i-1)+pa*rE*HEWh(i-1)-rA.*HAWh(i-1)-travelWh(i-1)*HAWh(i-1)/HNWh(i-1); % asymtomatic
    HIWh(i)=HIWh(i-1)+rP*HPWh(i-1)-rI.*HIWh(i-1)-rReportWh(i-1).*HIWh(i-1); % symtomatic
    HRWh(i)=HRWh(i-1)+rI.*HIWh(i-1)+rA.*HAWh(i-1); % removed
    IWh(i)=HIWh(i-1)*rReportWh(i-1); % daily confirmed cases
    HNWh(i)=HSWh(i)+HEWh(i)+HAWh(i)+HPWh(i)+HIWh(i)+HRWh(i)+sum(IWh); %all population
    HNWh(i)=popWh;
    
    %% other cities
    % ----beta----
    betaN(find(time2(i,:)==1))=betaBasicN(find(time2(i,:)==1)).*c2(find(time2(i,:)==1));
    
    % ----Imported Infected cases----
    importE_wh=HEWh(i-1)*flowWh(i-1,:)/popWh; %exposed cases from wuhan to other cities
    importA_wh=HAWh(i-1)*flowWh(i-1,:)/popWh; %asymptomatic cases from wuhan to other cities
    importP_wh=HPWh(i-1)*flowWh(i-1,:)/popWh; %pre-symptomatic cases from wuhan to other cities
    
    % ----flow----
     if(i>31 && i<61) 
                flow=flowFeb; %inter-city movement flow in FEB.
     elseif i>60
                flow=flowMar; %inter-city movement flow in MAR.  
     end

    % ----different class to calculate travel people proportion later----
    HLallcity=repmat(HE(i-1,:),n,1); %matrix 304*304,latent asymptomatic population of cities
    HIaallcity=repmat(HA(i-1,:),n,1); %matrix 304*304,infectious asymptomatic population of cities
    HLpallcity=repmat(HP(i-1,:),n,1); %matrix 304*304,presymptomatic population of cities
    HNallcity=repmat(HN(i-1,:),n,1);   %matrix 304*304,infectious symptomatic population of cities
    HSallcity=repmat(HS(i-1,:),n,1);   %matrix 304*304,suspecitble population of cities
    
    HLallcity(HLallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HIaallcity(HIaallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HLpallcity(HLpallcity<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    
    % ----infected population who travel between cities----
    inE=sum(flow.*HLallcity./HNallcity,2); %latent asymptomatic population who flow into the city
    outE=sum(flow.*HLallcity./HNallcity,1); %latent asymptomatic population who flow out of the city
    inA=sum(flow.*HIaallcity./HNallcity,2); %infectious asymptomatic population who flow into the city
    outA=sum(flow.*HIaallcity./HNallcity,1); %infectious asymptomatic population who flow out the city
    inP=sum(flow.*HLpallcity./HNallcity,2); %presymptomatic population who flow into the city
    outP=sum(flow.*HLpallcity./HNallcity,1); %presymptomatic population who flow out the city
    inS=sum(flow.*HSallcity./HNallcity,2);%suspectible population who flow into the city
    outS=sum(flow.*HSallcity./HNallcity,1);%suspectible population who flow out the city
    
    % --new daily cases---
    HAInfectious = HA(i-1,1:n);
    HPInfectious = HP(i-1,1:n);
    HIInfectious = HI(i-1,1:n);
    HAInfectious(HAInfectious<1)=0;
    HPInfectious(HPInfectious<1)=0;
    HIInfectious(HIInfectious<1)=0;
	tmpExporue = (r1*betaN.*HAInfectious+betaN.*HIInfectious+r2*betaN.*HPInfectious).*HS(i-1,:)./HN(i-1,:);
    
    % ----model----
    HS(i,:)=HS(i-1,:)-tmpExporue+inS'-outS; % suspecitble population
    HE(i,:)=HE(i-1,:)+tmpExporue+importE_wh-outE+inE'-rE*HE(i-1,:); % latent
    HP(i,:)=HP(i-1,:)+(1-pa)*rE*HE(i-1,:)-rP*HP(i-1,:)+importP_wh-outP+inP'; % presymtomatic
    HA(i,:)=HA(i-1,:)+pa*rE*HE(i-1,:)-rA.*HA(i-1,:)+importA_wh-outA+inA'; % asymtomatic
    HI(i,:)=HI(i-1,:)+rP*HP(i-1,:)-rI*HI(i-1,:)-rReportN(i-1,:).*HI(i-1,:); % symtomatic
    HR(i,:)=HR(i-1,:)+rI.*HI(i-1,:)+rA.*HA(i-1,:); % removed
    I(i,:)=rReportN(i-1,:).*HI(i-1,:); % daily confirmed cases
    HN(i,:)=HS(i,:)+HE(i,:)+HA(i,:)+HP(i,:)+HI(i,:)+HR(i,:)+sum(I(:,:),1); % all population
    HN(i,:)=pop;
end
ydot=[I(:,:),IWh(:)];

 
