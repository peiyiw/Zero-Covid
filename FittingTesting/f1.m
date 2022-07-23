clear
Case = xlsread('TonghuaDailyInfection_withRollingmeean3D.xlsx');

load Data_TH.mat %data 
data.ydata = [];
data.ydata(:,2) = Case(:,4);
data.ydata(:,3) = Case(:,6);


%%
% The model sum of squares given in the model structure.
model.ssfun = @f2;

%%
% All parameters are constrained to be positive. The initial
% concentrations are also unknown and are treated as extra parameters.
params1 = {    
      %parameter
      %name,mean,min.max,mean,std
      {'betan', 0.5,0,1}  %transmission rate 
      {'k',0.5,0,1} % the The fraction of contacts that can be successfully traced
      {'E0',10,0,50}
      {'A0',10,0,50}
      {'P0',10,0,50}
      {'I0',10,0,50}
      {'Tau0',0.1,0,1}
    };

%%
% We assume having at least some prior information on the
% repeatability of the observation and assign rather non informational
% prior for the residual variances of the observed states. The default
% prior distribution is sigma2 ~ invchisq(S20,N0), the inverse chi
% squared distribution (see for example Gelman et al.). The 3
% components (_A_, _Z_, _P_) all have separate variances.
model.S20 = [4];
model.N0  = [1];

City = 'Tonghua'

%%
% First generate an initial chain.
options.nsimu = 400000;
options.stats = 1;
[results, chain, s2chain,sschain]= mcmcrun(model,data,params1,options);

%regenerate chain to convergence
options.nsimu = 200000;
options.stats = 1;
[results2, chain2, s2chain2] = mcmcrun(model,data,params1,options,results);

save('results2.mat','results2');
save('chain2.mat','chain2');
save('s2chain2.mat','s2chain2');


%%
% Chain plots should reveal that the chain has converged and we can
% % use the results for estimation and predictive inference.
figure
mcmcplot(chain2,[],results2,'denspanel',2);
saveas(gcf,strcat(City,'_Chain2_posterior_probability'),'epsc');

figure
mcmcplot(chain2,[],results2); %,'pairs'
saveas(gcf,strcat(City,'_Chain2_MCMC'),'epsc');




