function pij = ContactTracing(A, P, I, tau, currentTime,rE,rA,rP,rI,pa,t0,tL,Sens1,Sens2)

iDay = size(A,1);
n = size(A,2);

currentTime(currentTime==0)=1;

tL = repmat(tL,1,n);
indexLessTL = find(currentTime<=tL);
tL(indexLessTL) = currentTime(indexLessTL)-1;

% if currentTime<=tL
% 	tL = currentTime;
% 	tL = tL-1;
% end

%--the transition rate from compartment K-----
%--Note the tau is changing in reality. To make the calculation easy, the average tau is used here. 
% tau = mean(tau);
tau = max(tau);
rE_star = rE+tau.*Sens1;
rA_star = rA+tau.*Sens2;
rP_star = rP+tau.*Sens2;
rI_star = rI+tau.*Sens2;

%--fraction of individuals who leave J that reack K---
qEA = pa.*rE./rE_star;
qEP = (1-pa).*rE./rE_star;
qPI = rP./rP_star;
qIR = rI./rI_star;
qAR = rA./rA_star;

%--transition probability function---------------
function a = P_AE(qEA,rE_star,rA_star,Time)
	a = qEA.*rE_star.*(exp(-rE_star.*Time)-exp(-rA_star.*Time))./(rA_star-rE_star);
end

function a = P_PE(qEP,rE_star,rP_star,Time)
	a = qEP.*rE_star.*(exp(-rE_star.*Time)-exp(-rP_star.*Time))./(rP_star-rE_star);
end

function a = P_IP(qPI,rI_star,rP_star,Time)
	a = qPI.*rP_star.*(exp(-rP_star.*Time)-exp(-rI_star.*Time))./(rI_star-rP_star);
end

function a = P_IE(qPI,qEP,rI_star,rE_star,rP_star,Time)
	a = qEP.*qPI.*rE_star.*rP_star.*((exp(-rE_star.*Time)-exp(-rP_star.*Time))./...
        (rE_star-rP_star)-(exp(-rE_star.*Time)-exp(-rI_star.*Time))./(rE_star-rI_star))./(rP_star-rI_star);
end

function a = P_EE(rE_star,Time)
	a = exp(-rE_star.*Time);
end

function a = P_AA(rA_star,Time)
	a = exp(-rA_star.*Time);
end

function a = P_PP(rP_star,Time)
	a = exp(-rP_star.*Time);
end

function a = P_II(rI_star,Time)
	a = exp(-rI_star.*Time);
end

%----contact tracing removing probability-----

%---removing from compartment A-----
vAE = zeros(1,n);
vAA = zeros(1,n);
vAP = zeros(1,n);
vAI = zeros(1,n);

for city = 1:n
    if tL(city)>1 && A(iDay,city)>0
	    for i = 1:tL(city)
    	    vAE(city) = vAE(city)+P_EE(rE_star(city),i+t0)*P_AA(rA_star(city),i)*A(iDay-i,city)/A(iDay,city);
		    vAA(city) = vAA(city)+P_AE(qEA(city),rE_star(city),rA_star(city),i+t0)*P_AA(rA_star(city),i)*A(iDay-i,city)/A(iDay,city);
		    vAP(city) = vAP(city)+P_PE(qEP(city),rE_star(city),rP_star(city),i+t0)*P_AA(rA_star(city),i)*A(iDay-i,city)/A(iDay,city);
		    vAI(city) = vAI(city)+P_IE(qPI(city),qEP(city),rI_star(city),rE_star(city),rP_star(city),i+t0)*P_AA(rA_star(city),i)*A(iDay-i,city)/A(iDay,city);
	    end
    else
	    vAE(city) = 1;
	    vAA(city) = 1;
	    vAP(city) = 1;
	    vAI(city) = 1;
    end
end


%---removing from compartment P-----
vPE = zeros(1,n);
vPA = zeros(1,n);
vPP = zeros(1,n);
vPI = zeros(1,n);

for city = 1:n
    if tL(city)>1 && P(iDay,city)>0
	    for i = 1:tL(city)
    	    vPE(city) = vPE(city)+P_EE(rE_star(city),i+t0)*P_PP(rP_star(city),i)*P(iDay-i,city)/P(iDay,city);
		    vPA(city) = vPA(city)+P_AE(qEA(city),rE_star(city),rA_star(city),i+t0)*P_PP(rP_star(city),i)*P(iDay-i,city)/P(iDay,city);
		    vPP(city) = vPP(city)+P_PE(qEP(city),rE_star(city),rP_star(city),i+t0)*P_PP(rP_star(city),i)*P(iDay-i,city)/P(iDay,city);
		    vPI(city) = vPI(city)+P_IE(qPI(city),qEP(city),rI_star(city),rE_star(city),rP_star(city),i+t0)*P_PP(rP_star(city),i)*P(iDay-i,city)/P(iDay,city);
	    end
    else  
	    vPE(city) = 1;
	    vPA(city) = 1;
	    vPP(city) = 1;
	    vPI(city) = 1;
    end
end
       
%---removing from compartment I-----

tmpIE1 = zeros(1,n);
tmpIE2 = zeros(1,n);
tmpIA1 = zeros(1,n);
tmpIA2 = zeros(1,n);
tmpIP1 = zeros(1,n);
tmpIP2 = zeros(1,n);
tmpII1 = zeros(1,n);
tmpII2 = zeros(1,n);

vIE = zeros(1,n);
vIA = zeros(1,n);
vIP = zeros(1,n);
vII = zeros(1,n);

Indicator = P+I;
Zero = Indicator==0;
for city = 1:n
    if tL(city)>1 && I(iDay,city)>0 && sum(Zero(end-currentTime(city)+1:end,city))==0
	    for i = 1:tL(city)
		    tmpIE1(city) = tmpIE1(city)+P_EE(rE_star(city),i+t0)*P_IP(qPI(city),rI_star(city),rP_star(city),i)*P(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=P---
		    tmpIE2(city) = tmpIE2(city)+P_EE(rE_star(city),i+t0)*P_II(rI_star(city),i)*I(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=I---
		    tmpIA1(city) = tmpIA1(city)+P_AE(qEA(city),rE_star(city),rA_star(city),i+t0)*P_IP(qPI(city),rI_star(city),rP_star(city),i)*P(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=P---
		    tmpIA2(city) = tmpIA2(city)+P_AE(qEA(city),rE_star(city),rA_star(city),i+t0)*P_II(rI_star(city),i)*I(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=I---
		    tmpIP1(city) = tmpIP1(city)+P_PE(qEP(city),rE_star(city),rP_star(city),i+t0)*P_IP(qPI(city),rI_star(city),rP_star(city),i)*P(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=P---
		    tmpIP2(city) = tmpIP2(city)+P_PE(qEP(city),rE_star(city),rP_star(city),i+t0)*P_II(rI_star(city),i)*I(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=I---
		    tmpII1(city) = tmpII1(city)+P_IE(qPI(city),qEP(city),rI_star(city),rE_star(city),rP_star(city),i+t0)*P_IP(qPI(city),rI_star(city),rP_star(city),i)*P(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=P---
		    tmpII2(city) = tmpII2(city)+P_IE(qPI(city),qEP(city),rI_star(city),rE_star(city),rP_star(city),i+t0)*P_II(rI_star(city),i)*I(iDay-i,city)^2/I(iDay,city)/(P(iDay-i,city)+I(iDay-i,city)); %--K=I---
	    end
	    vIE(city) = tmpIE1(city)+tmpIE2(city);
	    vIA(city) = tmpIA1(city)+tmpIA2(city);
	    vIP(city) = tmpIP1(city)+tmpIP2(city);
	    vII(city) = tmpII1(city)+tmpII2(city);
    else
	    vIE(city) = 1;
	    vIA(city) = 1;
	    vIP(city) = 1;
	    vII(city) = 1;
    end
end

%---rownames: A, P, I; colnames: E, A, P, I;------
tmp = vAE+vAA+vAP+vAI;
tmp(find(tmp==0))=1;

vAE = vAE./tmp;
vAA = vAA./tmp;
vAP = vAP./tmp;
vAI = vAI./tmp;

tmp = vPE+vPA+vPP+vPI;
tmp(find(tmp==0))=1;

vPE = vPE./tmp;
vPA = vPA./tmp;
vPP = vPP./tmp;
vPI = vPI./tmp;

tmp = vIE+vIA+vIP+vII;
tmp(find(tmp==0))=1;

vIE = vIE./tmp;
vIA = vIA./tmp;
vIP = vIP./tmp;
vII = vII./tmp;

pij = [vAE;vAA;vAP;vAI;vPE;vPA;vPP;vPI;vIE;vIA;vIP;vII];

end