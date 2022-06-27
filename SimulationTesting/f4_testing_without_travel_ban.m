function ydot = f4(t,theta,xdata,iCenter)

T = max(t); % days of simulation
dur = T; % days of epidemic
n = theta.n; % total number of cities

%% read param
% ---reproduction number---
R0=theta.R0;

% ----population level testing----
lagCenter = theta.lag0 ; % time lag in outbreak city
lagOther = theta.lag1; % time lag in other cities
breakCenter = theta.break0 ; % break interval in outbreak city
breakOther = theta.break1; % break interval in other cities
testingCenter = theta.testing0 ; % testing interval in outbreak city
testingOther = theta.testing1; % testing interval in other cities

% ---contact tracing---
k = theta.k; % fraction of contacts that can be successfully traced

% ---basic NPI---
betan = theta.betan; % percentage reduction in transmission rate due to other NPI

% ----vaccine---
ve = theta.ve; % vaccine effectiveness

% ---other---
iCenter = theta.iCenter; % epicenter
cityflow = theta.cityflow; % movement flow

% ---initial state---
perV=theta.perV; % vaccine coverage
E0 = theta.E0; % initial state of outbreak city
A0 = theta.A0;
P0 = theta.P0;
I0 = theta.I0;

%% param
% ----population level testing----
Sens1 = 0.4116; % the sensitivity of PCR when the individual in E status
Sens2 = 0.9438; % the sensitivity of PCR when the individual in A, P, I status

Tau=zeros(T,n);
Tau(1:lagCenter,iCenter)=0;
if breakCenter == 0
    Tau(lagCenter+1:end,iCenter)=testingCenter;
else
    dtesting = 1/testingCenter;
    dround = dtesting+breakCenter;
    nround = floor((T-lagCenter)/dround);
    tmp = rem((T-lagCenter),dround);
    Tau(lagCenter+1:lagCenter+nround*dround,iCenter)=repmat([repmat(testingCenter,dtesting,1);zeros(breakCenter,1)],nround,1);
    if tmp <= dtesting
        Tau(lagCenter+nround*dround+1:end,iCenter)=testingCenter;
    else
         Tau(lagCenter+nround*dround+1:lagCenter+nround*dround+dtesting,iCenter)=testingCenter;
         Tau(lagCenter+1+nround*dround+dtesting:end,iCenter)=0;
    end
end

% ---contact tracing---
t0 = 0; % contact tracing delay
q = 1./14; % Duration of quarantine
L = 14; % Pre-defined contact tracing time window

% ---outbreak flag---
bCase = zeros(T,n); % day of first imported case (binary)
bCase(1,iCenter) = 1;

currentTime(:,1:n) = zeros(1,n); % day of outbreak
currentTime(1,iCenter)=1;
startDay = zeros(1,n);
startDay(1,iCenter) = 1; % start day of outbreak

% ---population size---
pop(:,1:n)=xdata.pop(:,1:n);

% ---basic NPI---
M = 14.6;

% ----dynamic parameters----
pa=0.5; % asymptomatic proportion
rE = 1./2.9; % latent period
rA = 1./5; % recovery period for asymptomatic
rI = 1./5; % recovery period for symptomatic
rP = 1./2.3; % latent period but for pre-symptomaic
tL = 25; % the duration of COVID-19
r1 = 1./3.84; % Relative infectiousness of asymptomatic infections vs. symptomatic infections
r2 = 0.15; % Relative infectiousness of pre-symptomatic infections vs. symptomatic infections
betaI= R0./(pa.*r1./rA+(1-pa).*(r2./rP+1./rI)); % Transmission rate for symptomatic infections

betar=betaI.*betan; %transmission rate for S
betaA = r1.*betar;
betaP = r2.*betar;

betarV = betaI.*betan.*(1-ve); %transmission rate for V
betaAV = r1.*betarV;
betaPV = r2.*betarV;

w = 0; % immunity waned rate

%% intial state
% ----vaccination state----
HV(1,1:n)=perV.*pop; % vaccinated
HQV(1,1:n) = 0; % quarantined vaccinated

% ----initial state----
HE(1,1:n) = 0; % latent
HA(1,1:n) = 0; % asymptomatic 
HP(1,1:n) = 0; % pre-symptomatic
HI(1,1:n) = 0; % symptomatic 
HS(1,1:n) = pop; %susceptible
HN(1,1:n) = pop; %all population

HE(1,iCenter) = E0;
HA(1,iCenter) = A0;
HP(1,iCenter) = P0;
HI(1,iCenter) = I0;
HS(1,1:n) = HN(1,1:n)-HI(1,1:n)-HP(1,1:n)-HA(1,1:n)-HE(1,1:n)-HV(1,1:n);

HC(1,1:n) = 0; %the identified cases through the contact tracing.
HT(1,1:n) = 0; %the identified cases through the population-level-testing
DHC(1,1:n) = 0;
DHT(1,1:n) = 0;

HR(1,1:n) = 0; %removed infectious indivial
HQS(1,1:n) = 0; % quarantined suspectible

dailyCases(1,1:n)=HE(1,:)+HA(1,:)+HP(1,:)+HI(1,:);
dailyCasesInt(1,1:n)=HE(1,:)+HA(1,:)+HP(1,:)+HI(1,:);

%% transition probability
pij = ContactTracing(HA(1,1:n), HP(1,1:n), HI(1,1:n), Tau, currentTime, rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2);
pAE = pij(1,1:n);
pAA = pij(2,1:n);
pAP = pij(3,1:n);
pAI = pij(4,1:n);
pPE = pij(5,1:n);
pPA = pij(6,1:n);
pPP = pij(7,1:n);
pPI = pij(8,1:n);
pIE = pij(9,1:n);
pIA = pij(10,1:n);
pIP = pij(11,1:n);
pII = pij(12,1:n);

Cutoff = 1;

%% model
for i = 2:T 
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
    if sum(HV(i-1,1:n)<0)>0
            HV(i-1,find(HV(i-1,1:n)<0))=0;
    end
    if sum(HN(i-1,1:n)<0)>0
        HN(i-1,find(HN(i-1,1:n)<0))=0;
    end

    % ----average transmission rate----
    meanTau = mean(Tau);
    col=find(startDay~=0);
    if breakOther ~= 0
        meanTau(col) = 1/(1/testingOther+1/breakOther);
    else
        meanTau(col) = 1/(1/testingOther);
    end

    aveB = (betaP.*(rP+Sens2.*meanTau).^(-1)+betar.*(rI+Sens2.*meanTau).^(-1))./...
    ((rP+Sens2.*meanTau).^(-1)+(rI+Sens2.*meanTau).^(-1));
    aveBV = (betaPV.*(rP+Sens2.*meanTau).^(-1)+betarV.*(rI+Sens2.*meanTau).^(-1))./...
    ((rP+Sens2.*meanTau).^(-1)+(rI+Sens2.*meanTau).^(-1));
    currentTime(1,iCenter)=i;

    % ----cityflow----
    flow=cityflow; % daily inter city flow of current month
%     flow = zeros(n,1)
    
    % ----different class to calculate travel people proportion later----
    HEn=repmat(HE(i-1,:),n,1); %matrix 305*305,latent population of cities
    HAn=repmat(HA(i-1,:),n,1); %matrix 305*305,infectious asymptomatic population of cities
    HPn=repmat(HP(i-1,:),n,1); %matrix 305*305,presymptomatic population of cities
    HSn=repmat(HS(i-1,:),n,1);   %matrix 305*305,suspecitble population of cities
    HVn=repmat(HV(i-1,:),n,1);   %matrix 305*305,vaccinated population of cities
    HNn=repmat(HN(i-1,:),n,1);   %matrix 305*305,total population of cities
    
    HEn(HEn<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HAn(HAn<1)=0; %we assume that Less than one person does not have the ability to transmit virus
    HPn(HPn<1)=0; %we assume that Less than one person does not have the ability to transmit virus
   
    % ----infected population who travel between cities----
    inE=sum(flow.*HEn./HNn,2); %latent asymptomatic population who flow into the city
    outE=sum(flow.*HEn./HNn,1); %latent asymptomatic population who flow out of the city
    inA=sum(flow.*HAn./HNn,2); %infectious asymptomatic population who flow into the city
    outA=sum(flow.*HAn./HNn,1); %infectious asymptomatic population who flow out the city
    inP=sum(flow.*HPn./HNn,2); %presymptomatic population who flow into the city
    outP=sum(flow.*HPn./HNn,1); %presymptomatic population who flow out the city
    inS=sum(flow.*HSn./HNn,2);%suspectible population who flow into the city
    outS=sum(flow.*HSn./HNn,1);%suspectible population who flow out the city
    inV=sum(flow.*HVn./HNn,2);%vaccinated population who flow into the city
    outV=sum(flow.*HVn./HNn,1);%vaccinated population who flow out the city

    % --suspectible removed due to CT---
	uS = Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n)+...
    Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n)+...
    Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n)+...
    Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n);
	uS = uS.*HS(i-1,1:n)./HN(i-1,1:n);

    uV = Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n)+...
    Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n)+...
    Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n)+...
    Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n);
	uV = uV.*HV(i-1,1:n)./HN(i-1,1:n);
    
    % --new daily cases---
    HAInfectious = HA(i-1,1:n);
    HPInfectious = HP(i-1,1:n);
    HIInfectious = HI(i-1,1:n);
    HAInfectious(HAInfectious<1)=0;
    HPInfectious(HPInfectious<1)=0;
    HIInfectious(HIInfectious<1)=0;
	tmpExporueS = (betaA.*HAInfectious+betaP.*HPInfectious+betar.*HIInfectious).*HS(i-1,1:n)./HN(i-1,1:n);
    tmpExporueV = (betaAV.*HAInfectious+betaPV.*HPInfectious+betarV.*HIInfectious).*HV(i-1,1:n)./HN(i-1,1:n);
    dailyCases(i,:) = tmpExporueS+tmpExporueV;
    dailyCasesInt(i,:) = round(dailyCases(i,:));

    % ---end of the epidemic outbreak---
    if i>lagCenter && sum(dailyCasesInt(i-2:i,:),"all") < Cutoff
        dur = i-1; %no cases reported in three consecutive days
        break
    end
    
    % ---model---
	HS(i,1:n) = HS(i-1,1:n)-tmpExporueS+w.*HR(i-1,1:n)+q.*HQS(i-1,1:n)-uS+inS'-outS;
    HV(i,1:n) = HV(i-1,1:n)-tmpExporueV+q.*HQV(i-1,1:n)-uV+inV'-outV;
	HE(i,1:n) = HE(i-1,1:n)+tmpExporueV+tmpExporueS-(Sens1.*Tau(i-1,1:n)+rE).*HE(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n).*pIE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n).*pPE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n).*pAE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n).*pIE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n).*pPE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n).*pAE./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)-...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+inE'-outE;
	HA(i,1:n) = HA(i-1,1:n)+pa.*rE.*HE(i-1,1:n)-(Sens2.*Tau(i-1,1:n)+rA).*HA(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n).*pIA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n).*pPA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n).*pAA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n).*pIA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n).*pPA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n).*pAA./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)-...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+inA'-outA;
	HP(i,1:n) = HP(i-1,1:n)+(1-pa).*rE.*HE(i-1,1:n)-(Sens2.*Tau(i-1,1:n)+rP).*HP(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n).*pIP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n).*pPP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n).*pAP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n).*pIP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n).*pPP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n).*pAP./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)-...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+inP'-outP;
	HI(i,1:n) = HI(i-1,1:n)+rP.*HP(i-1,1:n)-(Sens2.*Tau(i-1,1:n)+rI).*HI(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n).*pII./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n).*pPI./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n).*pAI./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n).*pII./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n).*pPI./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n).*pAI./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)-...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)-...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n);
	HR(i,1:n) = HR(i-1,1:n)+rA.*HA(i-1,1:n)+rI.*HI(i-1,1:n)-w.*HR(i-1,1:n);
	HQS(i,1:n) = HQS(i-1,1:n)-q.*HQS(i-1,1:n)+uS;
    HQV(i,1:n) = HQV(i-1,1:n)-q.*HQV(i-1,1:n)+uV;
	HC(i,1:n) = HC(i-1,1:n)+Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n);
	HT(i,1:n) = HT(i-1,1:n)+Tau(i-1,1:n).*(Sens1.*HE(i-1,1:n)+Sens2.*HA(i-1,1:n)+...
        Sens2.*HP(i-1,1:n)+Sens2.*HI(i-1,1:n));
	HN(i,1:n) = HS(i-1,1:n)+HV(i-1,1:n)+HE(i-1,1:n)+HA(i-1,1:n)+HP(i-1,1:n)+HI(i-1,1:n)+HR(i-1,1:n)+HQS(i-1,1:n)+HQV(i-1,1:n)+HT(i-1,1:n)+HC(i-1,1:n); %all population
    HN(i,1:n) = pop;    

    % ---daily reported cases from contact tracing---
	DHC(i,1:n) = Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*aveB.*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaP.*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HS(i-1,1:n).*betaA.*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*aveBV.*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaPV.*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*L.*HV(i-1,1:n).*betaAV.*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HE(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HA(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HP(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*aveB./HN(i-1,1:n)./M-HV(i-1,1:n).*aveBV./HN(i-1,1:n)./M).*HI(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaP./HN(i-1,1:n)./M-HV(i-1,1:n).*betaPV./HN(i-1,1:n)./M).*HP(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens2.*Tau(i-1,1:n).*k.*M.*L.*(1-HS(i-1,1:n).*betaA./HN(i-1,1:n)./M-HV(i-1,1:n).*betaAV./HN(i-1,1:n)./M).*HA(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n)+...
        Sens1.*Tau(i-1,1:n).*k.*M.*L.*HE(i-1,1:n).*HI(i-1,1:n)./HN(i-1,1:n);

    % ---daily reported cases from population level testing---
	DHT(i,1:n) = Tau(i-1,1:n).*(Sens1.*HE(i-1,1:n)+Sens2.*HA(i-1,1:n)+Sens2.*HP(i-1,1:n)+Sens2.*HI(i-1,1:n));

    % ---transition probability---
	pij = ContactTracing(HA(1:i,1:n), HP(1:i,1:n), HI(1:i,1:n), Tau, currentTime,rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2);
	%--[vAE vAA vAP vAI; vPE vPA vPP vPI; vIE vIA vIP vII]--
    pAE = pij(1,1:n);
    pAA = pij(2,1:n);
    pAP = pij(3,1:n);
    pAI = pij(4,1:n);
    pPE = pij(5,1:n);
    pPA = pij(6,1:n);
    pPP = pij(7,1:n);
    pPI = pij(8,1:n);
    pIE = pij(9,1:n);
    pIA = pij(10,1:n);
    pIP = pij(11,1:n);
    pII = pij(12,1:n);

    % ---start of pupulation-testing in other cities---
    if sum(dailyCasesInt(i,1:n)>=Cutoff)>0
        bCase(i,find(dailyCasesInt(i,1:n)>=Cutoff))=1;
    end

    if sum(bCase(i,:) ~= bCase(i-1,:))> 0
        iOther = find((bCase(i,:)-bCase(i-1,:))==1);
        if find(startDay(iOther)==0) ~= 0
            startDay(iOther(find(startDay(iOther)==0))) = i;
        end
        if breakOther == 0
            Tau(lagOther+i:end,iOther)=testingOther;
        else
            dtesting = 1/testingOther;
            dround = dtesting+breakOther;
            nround = floor((T-lagOther-i+1)/dround);
            tmp = rem((T-lagOther-i+1),dround);
            Tau(lagOther+i:lagOther+nround*dround+i-1,iOther)=...
                repmat([repmat(testingOther,dtesting,1);zeros(breakOther,1)],nround,size(iOther,2));
            if tmp <= dtesting
                Tau(lagOther+nround*dround+i:end,iOther)=testingOther;
            else
                 Tau(lagOther+nround*dround+i:lagOther+nround*dround+dtesting+i-1,iOther)=testingOther;
                 Tau(lagOther+nround*dround+dtesting+i:end,iOther)=0;
            end
        end
    end
    currentTime(find(startDay~=0)) = i - startDay(find(startDay~=0)) + 1; 
end
i=dur;
bLockdown = (sum(DHC(:,:),1)./pop >= 0.8); %whether the city is lockdown

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
if sum(HV(i-1,1:n)<0)>0
	HV(i-1,find(HV(i-1,1:n)<0))=0;
end
if sum(HN(i-1,1:n)<0)>0
	HN(i-1,find(HN(i-1,1:n)<0))=0;
end

ydot.case=[dailyCases(:,:)]; % daily new cases
ydot.bLockdown=[bLockdown]; % whether the city is lockdown