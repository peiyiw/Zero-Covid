function [tsCaseAge, tsCaseTotal, tsQuarAge, tsQuarTotal, tsHosAge, tsHosTotal,tsICUAge, tsICUTotal, cumHAge, cumHTotal, cumUAge, cumUTotal, bLockdown] = f4(t,theta,xdata)

%% time span
T = max(t); %days of simulation

%% get data
pop = xdata.pop; %population
ageStructure = xdata.ageStructure; %age structure
cm = (xdata.cm)'; %contact matrix
flow = xdata.flow; % movement flow
% flow = zeros(366,366);

%% set parameter
% ---setting---
nG = 4; %number of age groups;
nC = 366; %total number of cities
iCenter = 168; %epicenter

% ---testing and tracing---
Sens1 = 0.4116; % the sensitivity of PCR when the individual in E status
Sens2 = 0.9438; % the sensitivity of PCR when the individual in A, P, I status
t0 = 0; % contact tracing delay
q = 1./14; % Duration of quarantine
L = 14; % Pre-defined contact tracing time window

% ---dynamic parameters---
sus = [1 1 1 1]; %susceptibility to infection for age group
tL = 25; % the duration of COVID-19
r1 = 1./3.84; % Relative infectiousness of asymptomatic infections vs. symptomatic infections
r2 = 0.15; % Relative infectiousness of pre-symptomatic infections vs. symptomatic infections

%% get parameter
% ---transmission---
ve = theta.ve; %vaccine effectiveness
npi = theta.npi; %percentage reduction in transmission rate due to other NPI

% ---testing and tracing---
k = theta.k; %fraction of contacts that can be successfully traced
lag = theta.lag; %time lag
interval = theta.interval ; %testing interval

% ---dynamic parameters---
R0 = theta.R0;
pa = theta.pa; % asymptomatic proportion
phs = theta.phs; %hospitalization rate for symptomatic Omicron-infection in unvaccinated individuals
phv = theta.phv; %hospitalization rate for symptomatic Omicron-infection in vaccinated individuals
pus = theta.pus; %icu rate for symptomatic Omicron-infection in unvaccinated hospitalized patients
puv = theta.puv; %icu rate for symptomatic Omicron-infection in vaccinated hospitalized patients
rE = theta.rE; % latent period
rA = theta.rA; % recovery period for asymptomatic
rI = theta.rI; % recovery period for symptomatic
rP = theta.rP; % latent period but for pre-symptomaic
rI2H = theta.rI2H; % 1/interval between symptom onset and admission
rI2U = theta.rI2U; % 1/interval between symptom onset and severe illness
rH = theta.rH; % 1/hospital length of stay
rU = theta.rU; % 1/icu length of stay


%% calcluation
eigValue = eigs(((repmat(pa,nG,1).*r1./rA+(1-repmat(pa,nG,1)).*(r2./rP+1./rI)).*repmat(sus,nG,1)'.*cm),1);
beta= R0./eigValue; % Transmission rate for symptomatic infections
[~, seedStructure] = CalcR0(ageStructure(:,1),sus,pa,cm,beta,r1,r2,rE,rP,rA,rI);

% ---initial state---
perV = theta.perV; % vaccine coverage
E0 = theta.E0; % initial state of outbreak city
A0 = theta.A0;
P0 = theta.P0;
I0 = theta.I0;

%% initialization
% ---initialization---
pAE = zeros(1,nG,nC);
pAA = zeros(1,nG,nC);
pAP = zeros(1,nG,nC);
pAI = zeros(1,nG,nC);
pPE = zeros(1,nG,nC);
pPA = zeros(1,nG,nC);
pPP = zeros(1,nG,nC);
pPI = zeros(1,nG,nC);
pIE = zeros(1,nG,nC);
pIA = zeros(1,nG,nC);
pIP = zeros(1,nG,nC);
pII = zeros(1,nG,nC);

% ---transmission---
betaI = beta.*(1-npi); %transmission rate for S
betaA = r1.*betaI;
betaP = r2.*betaI;
betaIV = beta.*(1-npi).*(1-ve); %transmission rate for V
betaAV = r1.*betaIV;
betaPV = r2.*betaIV;

% ---testing---
Tau = zeros(T,nG,nC); %testing rate
Tau(1:lag,1:nG,iCenter) = 0;
if interval ~=0
    Tau(lag+1:end,1:nG,iCenter) = 1./interval;
else
    Tau(lag+1:end,1:nG,iCenter) = 0;
end

% ---intial state---
HV(1,1:nG,1:nC) = perV.*pop.*ageStructure; % vaccinated
HE(1,1:nG,1:nC) = 0; % latent
HA(1,1:nG,1:nC) = 0; % asymptomatic 
HP(1,1:nG,1:nC) = 0; % pre-symptomatic
HI(1,1:nG,1:nC) = 0; % symptomatic
HH(1,1:nG,1:nC) = 0; % hospitalization
HU(1,1:nG,1:nC) = 0; % ICU
HR(1,1:nG,1:nC) = 0; %removed infectious individual
HN(1,1:nG,1:nC) = pop.*ageStructure; %all population

HE(1,1:nG,iCenter) = E0.*seedStructure; 
HA(1,1:nG,iCenter) = A0.*seedStructure;
HP(1,1:nG,iCenter) = P0.*seedStructure;
HI(1,1:nG,iCenter) = I0.*seedStructure;
HS(1,1:nG,1:nC) = HN(1,1:nG,1:nC)-HR(1,1:nG,1:nC)-HI(1,1:nG,1:nC)-HP(1,1:nG,1:nC)-HA(1,1:nG,1:nC)-HE(1,1:nG,1:nC)-HV(1,1:nG,1:nC);

HQV(1,1:nG,1:nC) = 0; % quarantined vaccinated
HQS(1,1:nG,1:nC) = 0; % quarantined suspectible
HC(1,1:nG,1:nC) = 0; %the identified cases through the contact tracing.
HT(1,1:nG,1:nC) = 0; %the identified cases through the population-level-testing
DHC(1,1:nG,1:nC) = 0;
DHT(1,1:nG,1:nC) = 0;

dailyCases(1,1:nG,1:nC) = HE(1,:,:)+HA(1,:,:)+HP(1,:,:)+HI(1,:,:);
dailyCasesInt(1,1:nG,1:nC) = round(HE(1,:,:))+round(HA(1,:,:))+round(HP(1,:,:))+round(HI(1,:,:));
cumQ(1,1:nG,1:nC) = 0; % cumulative quarantine
cumH(1,1:nG,1:nC) =0;
cumU(1,1:nG,1:nC) =0;

sumInE(1,1:nG,1:nC) = 0;
sumOutE(1,1:nG,1:nC) = 0;
sumInA(1,1:nG,1:nC) = 0;
sumOutA(1,1:nG,1:nC) = 0;
sumInP(1,1:nG,1:nC) = 0;
sumOutP(1,1:nG,1:nC) = 0;

% ---switch flag---
bCase = zeros(T,nC); % day with case (binary)
bCase(1,iCenter) = 1;
startDay = zeros(1,nC); %first day with imported case (number)
startDay(1,iCenter) = 1; 
currentTime(1,1:nC) = zeros(1,nC); %day of outbreak (number)
currentTime(1,iCenter) = 1;

%% transition probability
for g = 1:nG
	pij = ContactTracing(reshape(HA(1,g,1:nC),1,nC), ...
        reshape(HP(1,g,1:nC),1,nC), reshape(HI(1,g,1:nC),1,nC), Tau, currentTime,rE,rA,rP,rI,pa(g),t0,tL,Sens1,Sens2);
	%--[vAE vAA vAP vAI; vPE vPA vPP vPI; vIE vIA vIP vII]-----------------
	pAE(1,g,1:nC) = pij(1,1:nC);
	pAA(1,g,1:nC) = pij(2,1:nC);
	pAP(1,g,1:nC) = pij(3,1:nC);
	pAI(1,g,1:nC) = pij(4,1:nC);
	pPE(1,g,1:nC) = pij(5,1:nC);
	pPA(1,g,1:nC) = pij(6,1:nC);
	pPP(1,g,1:nC) = pij(7,1:nC);
	pPI(1,g,1:nC) = pij(8,1:nC);
	pIE(1,g,1:nC) = pij(9,1:nC);
	pIA(1,g,1:nC) = pij(10,1:nC);
	pIP(1,g,1:nC) = pij(11,1:nC);
	pII(1,g,1:nC) = pij(12,1:nC);
end

Cutoff = 1;

%% model
for i = 2:T
    for iC = 1:nC
        HA(i-1,(HA(i-1,:,iC)<0),iC) = 0;
        HI(i-1,(HI(i-1,:,iC)<0),iC) = 0;
        HP(i-1,(HP(i-1,:,iC)<0),iC) = 0;
        HE(i-1,(HE(i-1,:,iC)<0),iC) = 0;
        HS(i-1,(HS(i-1,:,iC)<0),iC) = 0;
        HV(i-1,(HV(i-1,:,iC)<0),iC) = 0;
    end

    % ---average transmission rate---
    if interval ~= 0
        meanTau = 1./interval;
    else
        meanTau = 0;
    end
    aveBetaA = betaA.*repmat(sus,nG,1)'.*cm;
    aveBetaP = betaP.*repmat(sus,nG,1)'.*cm;
    tmpI = (betaP.*(rP+Sens2.*meanTau).^(-1)+betaI.*(rI+Sens2.*meanTau).^(-1))./...
        ((rP+Sens2.*meanTau).^(-1)+(rI+Sens2.*meanTau).^(-1));
    aveBetaI = tmpI.*repmat(sus,nG,1)'.*cm;
    aveBetaAV = betaAV.*repmat(sus,nG,1)'.*cm;
    aveBetaPV = betaPV.*repmat(sus,nG,1)'.*cm;
    tmpIV = (betaPV.*(rP+Sens2.*meanTau).^(-1)+betaIV.*(rI+Sens2.*meanTau).^(-1))./...
        ((rP+Sens2.*meanTau).^(-1)+(rI+Sens2.*meanTau).^(-1));
    aveBetaIV = tmpIV*repmat(sus,nG,1)'.*cm;

    % ---different class to calculate travel people proportion later---
    %matrix 366*366*3
    HEn = repmat(permute(HE(i-1,1:nG,:),[1 3 2]),nC,1,1);
    HAn = repmat(permute(HA(i-1,1:nG,:),[1 3 2]),nC,1,1);
    HPn = repmat(permute(HP(i-1,1:nG,:),[1 3 2]),nC,1,1);
    HSn = repmat(permute(HS(i-1,1:nG,:),[1 3 2]),nC,1,1);
    HVn = repmat(permute(HV(i-1,1:nG,:),[1 3 2]),nC,1,1);
    HNn = repmat(permute(HN(i-1,1:nG,:),[1 3 2]),nC,1,1);
    
    HEn(HEn<1) = 0; %we assume that Less than one person does not have the ability to transmit virus
    HAn(HAn<1) = 0; %we assume that Less than one person does not have the ability to transmit virus
    HPn(HPn<1) = 0; %we assume that Less than one person does not have the ability to transmit virus
   
    
    % ---infected population who travel between cities---
    inE = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HEn./HNn,2); %latent asymptomatic population who flow into the city
    outE = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HEn./HNn,1); %latent asymptomatic population who flow out of the city
    inA = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HAn./HNn,2); %infectious asymptomatic population who flow into the city
    outA = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HAn./HNn,1); %infectious asymptomatic population who flow out the city
    inP = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HPn./HNn,2); %presymptomatic population who flow into the city
    outP = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HPn./HNn,1); %presymptomatic population who flow out the city
    inS = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HSn./HNn,2);%suspectible population who flow into the city
    outS = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HSn./HNn,1);%suspectible population who flow out the city
    inV = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HVn./HNn,2);%vaccinated population who flow into the city
    outV = sum(repmat(flow,1,1,nG).*reshape(ageStructure',1,366,nG).*HVn./HNn,1);%vaccinated population who flow out the city

    inE = permute(pagetranspose(inE),[1 3 2]);
    outE = permute(outE,[1 3 2]);
    inA = permute(pagetranspose(inA),[1 3 2]);
    outA = permute(outA,[1 3 2]);
    inP = permute(pagetranspose(inP),[1 3 2]);
    outP = permute(outP,[1 3 2]);
    inS = permute(pagetranspose(inS),[1 3 2]);
    outS = permute(outS,[1 3 2]);
    inV = permute(pagetranspose(inV),[1 3 2]);
    outV = permute(outV,[1 3 2]);

    sumInE(i,:,:) = sumInE(i-1,:,:)+inE;
    sumOutE(i,:,:) = sumOutE(i-1,:,:)+outE;
    sumInA(i,:,:) = sumInA(i-1,:,:)+inA;
    sumOutA(i,:,:) = sumOutA(i-1,:,:)+outA;
    sumInP(i,:,:) = sumInP(i-1,:,:)+inP;
    sumOutP(i,:,:) = sumOutP(i-1,:,:)+outP;

    % ---tracing---
    precisionA = pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaA./cm)+...
        pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaAV./cm);
	precisionP = pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaP./cm)+...
        pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaPV./cm);
	precisionI = pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaI./cm)+...
        pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*(aveBetaIV./cm);

   
    % --suspectible removed due to CT---
	u = pagemtimes(Sens2.*Tau(i-1,1:nG,1:nC).*k.*cm.*L.*(1-precisionI),pagetranspose(HI(i-1,1:nG,1:nC)))+...
    pagemtimes(Sens2.*Tau(i-1,1:nG,1:nC).*k.*cm.*L.*(1-precisionP),pagetranspose(HP(i-1,1:nG,1:nC)))+...
    pagemtimes(Sens2.*Tau(i-1,1:nG,1:nC).*k.*cm.*L.*(1-precisionA),pagetranspose(HA(i-1,1:nG,1:nC)))+...
    pagemtimes(Sens1.*Tau(i-1,1:nG,1:nC).*k.*cm.*L,pagetranspose(HE(i-1,1:nG,1:nC)));
	uS = pagetranspose(u.*pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)));
    uE = pagetranspose((u.*pagetranspose(HE(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC))));
    uA = pagetranspose((u.*pagetranspose(HA(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC))));
    uP = pagetranspose((u.*pagetranspose(HP(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC))));
    uI = pagetranspose((u.*pagetranspose(HI(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC))));
	uV = pagetranspose(u.*pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)));

    % --new daily cases---
    HAInfectious = HA(i-1,1:nG,1:nC);
    HPInfectious = HP(i-1,1:nG,1:nC);
    HIInfectious = HI(i-1,1:nG,1:nC);
    HAInfectious(HAInfectious<1) = 0;
    HPInfectious(HPInfectious<1) = 0;
    HIInfectious(HIInfectious<1) = 0;
	tmpExporueS = pagemtimes(betaA.*pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HAInfectious))+...
        pagemtimes(betaP.*pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HPInfectious))+...
        pagemtimes(betaI.*pagetranspose(HS(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HIInfectious));
    tmpExporueV = pagemtimes(betaAV.*pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HAInfectious))+...
        pagemtimes(betaPV.*pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HPInfectious))+...
        pagemtimes(betaIV.*pagetranspose(HV(i-1,1:nG,1:nC)./HN(i-1,1:nG,1:nC)).*sus'.*cm,pagetranspose(HIInfectious));
    tmpExporueS = pagetranspose(tmpExporueS);
    tmpExporueV = pagetranspose(tmpExporueV);
    dailyCases(i,:,:) = tmpExporueS+tmpExporueV;
    dailyCasesInt(i,:,:) = round(dailyCases(i,:,:));
  
    % ---hospitalization and icu---
    cumH = cumH+(tmpExporueS.*(1-pa).*phs.*(1-pus)+tmpExporueV.*(1-pa).*phv.*(1-puv));
    cumU = cumU+(tmpExporueS.*(1-pa).*phs.*pus+tmpExporueV.*(1-pa).*phv.*puv);
    HH(i,:,:) = HH(i-1,1:nG,1:nC)+(tmpExporueS.*(1-pa).*phs.*(1-pus)+tmpExporueV.*(1-pa).*phv.*(1-puv)).*(1./(1./rE+1./rP+1./rI2H))-rH.*HH(i-1,1:nG,1:nC); % daily hospitalization
    HU(i,:,:) = HU(i-1,1:nG,1:nC)+(tmpExporueS.*(1-pa).*phs.*pus+tmpExporueV.*(1-pa).*phv.*puv).*(1./(1./rE+1./rP+1./rI2U))-rU.*HU(i-1,1:nG,1:nC); % daily icu
    
    % ---end of the epidemic outbreak---
    if i>lag && sum(dailyCasesInt(i-6:i,:,:),"all") < Cutoff
        break % no cases reported in 7 consecutive days
    end
    
    % ---model---
	HS(i,:,:) = HS(i-1,1:nG,1:nC)-tmpExporueS+q.*HQS(i-1,1:nG,1:nC)-uS+inS-outS;
    HV(i,:,:) = HV(i-1,1:nG,1:nC)-tmpExporueV+q.*HQV(i-1,1:nG,1:nC)-uV+inV-outV;
	HE(i,:,:) = HE(i-1,1:nG,1:nC)+tmpExporueV+tmpExporueS-...
        (Sens1.*Tau(i-1,1:nG,1:nC)+(1-Sens1.*Tau(i-1,1:nG,1:nC)).*rE).*HE(i-1,1:nG,1:nC)-...
         Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC).*pIE)))-...
         Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC).*pPE)))-...
         Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC).*pAE)))-...
         uE+inE-outE;
	HA(i,:,:) = HA(i-1,1:nG,1:nC)+pa.*(1-Sens1.*Tau(i-1,1:nG,1:nC)).*rE.*HE(i-1,1:nG,1:nC)-...
        (Sens2.*Tau(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rA).*HA(i-1,1:nG,1:nC)-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC).*pIA)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC).*pPA)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC).*pAA)))-...
        uA+inA-outA;
	HP(i,:,:) = HP(i-1,1:nG,1:nC)+(1-pa).*(1-Sens1.*Tau(i-1,1:nG,1:nC)).*rE.*HE(i-1,1:nG,1:nC)-...
        (Sens2.*Tau(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rP).*HP(i-1,1:nG,1:nC)-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC).*pIP)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC).*pPP)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC).*pAP)))-...
        uP+inP-outP;
	HI(i,:,:) = HI(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rP.*HP(i-1,1:nG,1:nC)-...
        (Sens2.*Tau(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rI).*HI(i-1,1:nG,1:nC)-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC).*pII)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC).*pPI)))-...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC).*pAI)))-...
        uI;
	HR(i,:,:) = HR(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rA.*HA(i-1,1:nG,1:nC)+(1-Sens2*Tau(i-1,1:nG,1:nC)).*rI.*HI(i-1,1:nG,1:nC);
	HQS(i,:,:) = HQS(i-1,1:nG,1:nC)-q.*HQS(i-1,1:nG,1:nC)+uS;
    HQV(i,:,:) = HQV(i-1,1:nG,1:nC)-q.*HQV(i-1,1:nG,1:nC)+uV;
	HC(i,:,:) = HC(i-1,1:nG,1:nC)+...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC))))+...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC))))+...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC))))+...
        uE+uA+uP+uI;
    HT(i,:,:) = HT(i-1,1:nG,1:nC)+Tau(i-1,1:nG,1:nC).*(Sens1.*HE(i-1,1:nG,1:nC)+Sens2.*HA(i-1,1:nG,1:nC)+...
            Sens2.*HP(i-1,1:nG,1:nC)+Sens2.*HI(i-1,1:nG,1:nC));
    HN(i,:,:) = HS(i-1,1:nG,1:nC)+HV(i-1,1:nG,1:nC)+HE(i-1,1:nG,1:nC)+HA(i-1,1:nG,1:nC)+HP(i-1,1:nG,1:nC)+...
        HI(i-1,1:nG,1:nC)+HR(i-1,1:nG,1:nC)+HQS(i-1,1:nG,1:nC)+HQV(i-1,1:nG,1:nC)+...
        HT(i-1,1:nG,1:nC)+HC(i-1,1:nG,1:nC); %all population
    HN(i,:,:) = pop.*ageStructure;    
    
	DHC(i,:,:) = Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionI),pagetranspose(HI(i-1,1:nG,1:nC))))+...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionP),pagetranspose(HP(i-1,1:nG,1:nC))))+...
        Sens2.*Tau(i-1,1:nG,1:nC).*k.*L.*pagetranspose(pagemtimes((cm.*precisionA),pagetranspose(HA(i-1,1:nG,1:nC))))+...
        uE+uA+uP+uI; %daily reported cases from contact tracing 
	DHT(i,:,:) = Tau(i-1,1:nG,1:nC).*(Sens1.*HE(i-1,1:nG,1:nC)+Sens2.*HA(i-1,1:nG,1:nC)+...
        Sens2.*HP(i-1,1:nG,1:nC)+Sens2.*HI(i-1,1:nG,1:nC));  %daily reported cases from population level testing
    cumQ(i,:,:) = cumQ(i-1,1:nG,1:nC)+DHC(i,1:nG,1:nC)+DHT(i,1:nG,1:nC)+uS+uV; %quarantined population

    % ---start of pupulation-testing in other cities---
    bCase(i,reshape(sum(dailyCasesInt(i,1:nG,1:nC),2)>=Cutoff,1,iC)) = 1; %switch of case
    iOtherFirst = find(((bCase(i,:)-bCase(i-1,:))==1)&(startDay==0)); %cities with new case for the first time
    startDay(iOtherFirst)=i; %firt day with imported case
    if interval ~=0
        Tau(lag+i:end,:,iOtherFirst) = 1./interval; %testing rate in outbreak city
    else
        Tau(lag+i:end,:,iOtherFirst) = 0; %testing rate in outbreak city
    end
    currentTime(startDay~=0) = i - startDay(startDay~=0) + 1; 

    % ---transition probability---
    for g = 1:nG
	    pij = ContactTracing(reshape(HA(1:i,g,1:nC),i,nC), ...
            reshape(HP(1:i,g,1:nC),i,nC), reshape(HI(1:i,g,1:nC),i,nC), Tau, currentTime,rE,rA,rP,rI,pa(g),t0,tL,Sens1,Sens2);
	    %--[vAE vAA vAP vAI; vPE vPA vPP vPI; vIE vIA vIP vII]-----------------
	    pAE(1,g,1:nC) = pij(1,1:nC);
	    pAA(1,g,1:nC) = pij(2,1:nC);
	    pAP(1,g,1:nC) = pij(3,1:nC);
	    pAI(1,g,1:nC) = pij(4,1:nC);
	    pPE(1,g,1:nC) = pij(5,1:nC);
	    pPA(1,g,1:nC) = pij(6,1:nC);
	    pPP(1,g,1:nC) = pij(7,1:nC);
	    pPI(1,g,1:nC) = pij(8,1:nC);
	    pIE(1,g,1:nC) = pij(9,1:nC);
	    pIA(1,g,1:nC) = pij(10,1:nC);
	    pIP(1,g,1:nC) = pij(11,1:nC);
	    pII(1,g,1:nC) = pij(12,1:nC);
    end
end



tsCaseAge = dailyCasesInt(:,:,:); %daily cases of each age group
tsCaseTotal = reshape(sum(dailyCasesInt(:,:,:),2),size(dailyCasesInt,1),nC); % daily new cases of all age group
tsQuarAge =  HQS+HQV+HC+HT; %quarantined individuals of each age group
tsQuarTotal =   reshape(sum(HQS+HQV+HC+HT,2),size(tsQuarAge,1),nC); %quarantined individuals of all age group
tsHosAge = HH; %quarantined individuals of each age group
tsHosTotal = reshape(sum(HH,2),size(HH,1),nC); %quarantined individuals of each age group
tsICUAge = HU; %quarantined individuals of each age group
tsICUTotal = reshape(sum(HU,2),size(HU,1),nC); %quarantined individuals of each age group
cumHAge = cumH;  %hospitalized individuals of each age group
cumHTotal = reshape(sum(cumH,2),size(cumH,1),nC);  %hospitalized individuals of all age group
cumUAge = cumU; %icu individuals of each age group
cumUTotal = reshape(sum(cumU,2),size(cumU,1),nC);  %icu individuals of all age group
bLockdown = (max(tsQuarTotal(:,:),[],1)./pop >= 0.8); %whether the city is lockdown
