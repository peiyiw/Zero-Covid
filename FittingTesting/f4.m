function ydot = f4(t,theta,xdata)

T = 27; % the days of outbreak

Sens1 = 0.4116; % the sensitivity of PCR when the individual in E status
Sens2 = 0.9438; % the sensitivity of PCR when the individual in A, P, I status

%---population level testing rate----
Tau(1:3) = theta(7);
Tau(4:7) = 1.0/4; %the first round. Jan 15, 2021---Jan 18, 2021.
Tau(8) = theta(7);
Tau(9:11) = 1.0/3; %the second round. the start date is Jan 20, 2021.
Tau(12:13) = theta(7);
Tau(14:16) = 1.0/3; % the third round. the start date is Jan 25, 2021.
Tau(17:18) = theta(7);
Tau(19:27) = theta(7); % 



%--population size------
Ptonghua = 2147400; %population of Tonghua

%--dynamic parameters----
pa = 0.284;  %proportion of asymptomatic
rE = 1/2.9; % latent period
rA = 1/5; % recovery period for asymptomatic
rI = 1/5; % recovery period for symptomatic
rP = 1/2.3; % latent period but for pre-symptomaic
tL = 25; % the duration of COVID-19
r1 = 1/3.84; % Relative infectiousness of asymptomatic infections vs. symptomatic infections
r2 = 0.15; % Relative infectiousness of pre-symptomatic infections vs. symptomatic infections

w = 0;

%--NPI parameters--------
t0 = 0; % contact tracing delay
q = 1/14; % Duration of quarantine
L = 14; % Pre-defined contact tracing time window
M = 14.6;
VE=0;%the VE for two-dose vaccination against Delta was 59.0%
pV=0;%the two dose vaccined proportion of XiAn
R0 = 4.6;
betan=theta(1);
%--unknown parameters----
betaI = R0/(pa*r1/rA+(1-pa)*(r2/rP+1/rI));
betar= (1-betan)*betaI; % Transmission rate for symptomatic infections
betar1=betar*(1-VE);%Transmission rate for symptomatic infections in real world
k = theta(2); % fraction of contacts that can be successfully traced
HE(1) = theta(3); % latent
HA(1) = theta(4); % asymptomatic 
HP(1) = theta(5); % pre-symptomatic
HI(1) = theta(6); % symptomatic 


HS(1) = Ptonghua*(1-pV); %susceptible
HV(1) = Ptonghua*pV;
HN(1) = Ptonghua; %all population
HC(1) = 5; %the identified cases through the contact tracing.
HT(1) = 0; %the identified cases through the population-level-testing
DHC(1) = 5;
DHT(1) = 0;

HR(1) = 0; %removed infectious indivial
HQS(1) = 0;
HQV(1) = 0;

betaA = r1*betar;
betaP = r2*betar;
aveB = (betaP*(rP+Sens2*mean(Tau))^(-1)+betar*(rI+Sens2*mean(Tau))^(-1))/((rP+Sens2*mean(Tau))^(-1)+(rI+Sens2*mean(Tau))^(-1));
betaA1 = r1*betar1;
betaP1 = r2*betar1;
aveB1 = (betaP1*(rP+Sens2*mean(Tau))^(-1)+betar1*(rI+Sens2*mean(Tau))^(-1))/((rP+Sens2*mean(Tau))^(-1)+(rI+Sens2*mean(Tau))^(-1));

Cutoff = 1;

pij = ContactTracing(HA(1), HP(1), HI(1), Tau, 1,rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2);
pAE = pij(1,1);
pAA = pij(1,2);
pAP = pij(1,3);
pAI = pij(1,4);
pPE = pij(2,1);
pPA = pij(2,2);
pPP = pij(2,3);
pPI = pij(2,4);
pIE = pij(3,1);
pIA = pij(3,2);
pIP = pij(3,3);
pII = pij(3,4);
for i = 2:T 
if HA(i-1)<Cutoff
		HA(i-1)=0;
	end
	if HI(i-1)<Cutoff
		HI(i-1)=0;
	end
	if HP(i-1)<Cutoff
		HP(i-1)=0;
	end
	if HE(i-1)<Cutoff
		HE(i-1)=0;
    end
    if HS(i-1)<Cutoff
		HS(i-1)=0;
    end
    if HN(i-1)<Cutoff
		HN(i-1)=0;
    end
    if HV(i-1)<Cutoff
		HV(i-1)=0;
	end
	uS_sum = Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)...,
        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)...,
        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)...,
        +Sens1*Tau(i-1)*k*M*L*HE(i-1);
	uS = uS_sum*HS(i-1)/HN(i-1);
    uS1 =uS_sum*HV(i-1)/HN(i-1);

	tmp_Exporue = (betaA.*HA(i-1)+betaP.*HP(i-1)+betar.*HI(i-1))*HS(i-1)/HN(i-1);
    tmp_Exporue1 = (betaA1.*HA(i-1)+betaP1.*HP(i-1)+betar1.*HI(i-1))*HV(i-1)/HN(i-1);
	HS(i) = HS(i-1)-tmp_Exporue+w*HR(i-1)+q*HQS(i-1)-uS;
    HV(i) = HV(i-1)-tmp_Exporue1+w*HR(i-1)+q*HQV(i-1)-uS1;
	HE(i) = HE(i-1)+tmp_Exporue+tmp_Exporue1-(Sens1*Tau(i-1)+rE)*HE(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)*pIE/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)*pPE/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)*pAE/HN(i-1) ...,
	       -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HE(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HE(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HE(i-1)/HN(i-1)...,
           -Sens1*Tau(i-1)*k*M*L*HE(i-1)*HE(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)*pIE/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)*pPE/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)*pAE/HN(i-1);
	
	HA(i) = HA(i-1)+pa*rE*HE(i-1)-(Sens2*Tau(i-1)+rA)*HA(i-1)...,
    -Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)*pIA/HN(i-1)...,
    -Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)*pPA/HN(i-1)...,
    -Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)*pAA/HN(i-1) ...,
	-Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HA(i-1)/HN(i-1)...,
    -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HA(i-1)/HN(i-1)...,
    -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HA(i-1)/HN(i-1)...,
    -Sens1*Tau(i-1)*k*M*L*HE(i-1)*HA(i-1)/HN(i-1)...,
    -Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)*pIA/HN(i-1)...,
    -Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)*pPA/HN(i-1)...,
    -Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)*pAA/HN(i-1);


	HP(i) = HP(i-1)+(1-pa)*rE*HE(i-1)-(Sens2*Tau(i-1)+rP)*HP(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)*pIP/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)*pPP/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)*pAP/HN(i-1) ...,
           -Sens1*Tau(i-1)*k*M*L*HE(i-1)*HP(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)*pIP/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)*pPP/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)*pAP/HN(i-1) ...,
	       -Sens2*Tau(i-1)*k*M*L*(1-HV(i-1)*aveB1/HN(i-1)/M-HS(i-1)*aveB/HN(i-1)/M)*HI(i-1)*HP(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HV(i-1)*betaP1/HN(i-1)/M-HS(i-1)*betaP/HN(i-1)/M)*HP(i-1)*HP(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HV(i-1)*betaA1/HN(i-1)/M-HS(i-1)*betaA/HN(i-1)/M)*HA(i-1)*HP(i-1)/HN(i-1);
	
	
	HI(i) = HI(i-1)+rP*HP(i-1)-(Sens2*Tau(i-1)+rI)*HI(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)*pII/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)*pPI/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)*pAI/HN(i-1) ...,
	       -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HI(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HI(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HI(i-1)/HN(i-1)...,
           -Sens1*Tau(i-1)*k*M*L*HE(i-1)*HI(i-1)/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)*pII/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)*pPI/HN(i-1)...,
           -Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)*pAI/HN(i-1) ;
	
	HR(i) = HR(i-1)+rA*HA(i-1)+rI*HI(i-1)-w*HR(i-1);
    HQS(i) = HQS(i-1)-q*HQS(i-1)+uS;
    HQV(i) = HQV(i-1)-q*HQV(i-1)+uS1;
	HC(i) = HC(i-1)...,
            +Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)/HN(i-1) ...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HE(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HE(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HE(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HE(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HA(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HA(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HA(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HA(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HP(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HP(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HI(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HI(i-1)/HN(i-1);

	
	HT(i) = HT(i-1)+Tau(i-1)*(Sens1*HE(i-1)+Sens2*HA(i-1)+Sens2*HP(i-1)+Sens2*HI(i-1));

	%all population
	HN(i) = HS(i-1)+HE(i-1)+HA(i-1)+HP(i-1)+HI(i-1)+HR(i-1)+HQS(i-1)+HQV(i-1)+HT(i-1)+HC(i-1)+HV(i-1); 
	% the daily cases
	DHC(i) = Sens2*Tau(i-1)*k*L*HS(i-1)*aveB*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HS(i-1)*betaP*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HS(i-1)*betaA*HA(i-1)/HN(i-1) ...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*aveB1*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*betaP1*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*L*HV(i-1)*betaA1*HA(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HE(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HE(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HE(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HE(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HA(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HA(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HA(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HA(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HP(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HP(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HP(i-1)/HN(i-1) ...,
	        +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*aveB/HN(i-1)/M-HV(i-1)*aveB1/HN(i-1)/M)*HI(i-1)*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaP/HN(i-1)/M-HV(i-1)*betaP1/HN(i-1)/M)*HP(i-1)*HI(i-1)/HN(i-1)...,
            +Sens2*Tau(i-1)*k*M*L*(1-HS(i-1)*betaA/HN(i-1)/M-HV(i-1)*betaA1/HN(i-1)/M)*HA(i-1)*HI(i-1)/HN(i-1)...,
            +Sens1*Tau(i-1)*k*M*L*HE(i-1)*HI(i-1)/HN(i-1);
	DHT(i) = Tau(i-1)*(Sens1*HE(i-1)+Sens2*HA(i-1)+Sens2*HP(i-1)+Sens2*HI(i-1));
     if DHC(i)<0
		DHC(i)=0;
	end
	if DHT(i)<0
		DHT(i)=0;
	end

	pij = ContactTracing(HA(1:i), HP(1:i), HI(1:i), Tau, i,rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2);
	%--[vAE vAA vAP vAI; vPE vPA vPP vPI; vIE vIA vIP vII]-----------------
	pAE = pij(1,1);
	pAA = pij(1,2);
	pAP = pij(1,3);
	pAI = pij(1,4);
	pPE = pij(2,1);
	pPA = pij(2,2);
	pPP = pij(2,3);
	pPI = pij(2,4);
	pIE = pij(3,1);
	pIA = pij(3,2);
	pIP = pij(3,3);
	pII = pij(3,4);
end
i=T+1;
if HA(i-1)<Cutoff
		HA(i-1)=0;
	end
	if HI(i-1)<Cutoff
		HI(i-1)=0;
	end
	if HP(i-1)<Cutoff
		HP(i-1)=0;
	end
	if HE(i-1)<Cutoff
		HE(i-1)=0;
    end
    if HS(i-1)<Cutoff
		HS(i-1)=0;
    end
    if HN(i-1)<Cutoff
		HN(i-1)=0;
    end
    if HV(i-1)<Cutoff
		HV(i-1)=0;
    end
ydot=[DHT(:),DHC(:)]; % the daily obervations











