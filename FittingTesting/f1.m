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
% options.nsimu = 400000;
options.nsimu = 400000;
options.stats = 1;
[results, chain, s2chain,sschain]= mcmcrun(model,data,params1,options);

%regenerate chain to convergence
options.nsimu = 200000;
% options.nsimu = 200000;
options.stats = 1;
[results2, chain2, s2chain2] = mcmcrun(model,data,params1,options,results);

save('results2.mat','results2');
save('chain2.mat','chain2');
save('s2chain2.mat','s2chain2');


%%
% Chain plots should reveal that the chain has converged and we can
% % use the results for estimation and predictive inference.
load('chain2.mat');
load('results2.mat');
load('s2chain2.mat');
City = 'Tonghua'

  
figure
mcmcplot(chain2,[],results2,'denspanel',2);
saveas(gcf,strcat(City,'_Chain2_posterior_probability'),'epsc');

figure
mcmcplot(chain2,[],results2); %,'pairs'
saveas(gcf,strcat(City,'_Chain2_MCMC'),'epsc');

%%
% Function |chainstats| calculates mean ans std from the chain and
% estimates the Monte Carlfigure
% the integrated autocorrelation time and |geweke| is a simple test
% for a null hypothesis that the chain has converged.

results2.sstype = 1; % needed for mcmcpred and sqrt transformation
tmp = chainstats(chain2,results2); %statistic results of parameter estimation

Names = ["betan";"k";"E0";"A0";"P0";"I0";"Tau0"];
T = table(Names,tmp(:,1),tmp(:,2),tmp(:,3),tmp(:,4),tmp(:,5), 'VariableNames', { 'Parameter','Mean', 'std','MC_err','tau','geweke'} );
writetable(T, strcat(City,'_chain2_parameter_stats.csv'));

writematrix(tmp, strcat(City,'_chain_parameter_stats.csv'));



%%
% In order to use the |mcmcpred| function we need
% function |modelfun| with input arguments given as
% |modelfun(xdata,theta)|. We construct this as an anonymous function.
% load('chain2.mat');
% load('results2.mat');
% load('s2chain2.mat');

modelfun = @(d,th) f3(d(:,1),th,d);

%%o
% We sample 1000 parameter realizations from |chain| and |s2chain|
% and calculate the predictive plots.
nsample = 1000;
results2.sstype = 1;
out = mcmcpred(results2,chain2,s2chain2,data.xdata,modelfun,nsample);%data.ydata-->data

time=27;
% DHT
figure
subplot(2,1,1)
fillyy(1:time,out.obslims{1,1}{1,1}(3,:),out.obslims{1,1}{1,1}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,2),'.r');
plot(out.obslims{1,1}{1,1}(2,:));
hold off
ylabel('DHT')
%DHC
subplot(2,1,2)
fillyy(1:time,out.obslims{1,1}{1,2}(3,:),out.obslims{1,1}{1,2}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,3),'.r');
plot(out.obslims{1,1}{1,2}(2,:));
hold off
ylabel('DHC')

saveas(gcf,strcat(City,'_Fitting_wide'),'epsc');



figure

plot(data.ydata(:,2),'.r');
hold on;
plot(out.obslims{1,1}{1,1}(2,:),'r');

plot(data.ydata(:,3),'.b');
plot(out.obslims{1,1}{1,2}(2,:),'b');
hold off
ylabel('DHC')




DownCI = transpose(out.obslims{1,1}{1,1}(1,:));
Mean = transpose(out.obslims{1,1}{1,1}(2,:));
UpCI = transpose(out.obslims{1,1}{1,1}(3,:));

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} );
writetable(T, strcat(City,'_fitting_PopulationTesting_wide.csv'));

DownCI = transpose(out.obslims{1,1}{1,2}(1,:));
Mean = transpose(out.obslims{1,1}{1,2}(2,:));
UpCI = transpose(out.obslims{1,1}{1,2}(3,:));

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} );
writetable(T, strcat(City,'_fitting_ContactTracing_wide.csv'));





figure
subplot(2,1,1)
fillyy(1:time,out.predlims{1,1}{1,1}(3,:),out.predlims{1,1}{1,1}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,2),'.r');
plot(out.predlims{1,1}{1,1}(2,:));
hold off
ylabel('DHT')
%DHC
subplot(2,1,2)
fillyy(1:time,out.predlims{1,1}{1,2}(3,:),out.predlims{1,1}{1,2}(1,:),[0.8 0.8 0.8]);
hold on 
plot(data.ydata(:,3),'.r');
plot(out.predlims{1,1}{1,2}(2,:));
hold off
ylabel('DHC')
saveas(gcf,strcat(City,'_Fitting_narrow'),'epsc');


CIFcn = @(x,p)prctile(x,abs([0,100]-(100-p)/2));

Mean = transpose(mean(chain2,1));
tmp=transpose(CIFcn(chain2,95));
UpCI = tmp(:,2);
DownCI = tmp(:,1);
T = table(Names,DownCI,Mean,UpCI, 'VariableNames', {'Parameter' ,'DownCI', 'Mean','UpCI'} );

writetable(T, strcat(City,'_fitting_parameter.csv'));




DownCI = transpose(out.predlims{1,1}{1,1}(1,:));
Mean = transpose(out.predlims{1,1}{1,1}(2,:));
UpCI = transpose(out.predlims{1,1}{1,1}(3,:));

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} );
writetable(T, strcat(City,'_fitting_PopulationTesting_narrow.csv'));

DownCI = transpose(out.predlims{1,1}{1,2}(1,:));
Mean = transpose(out.predlims{1,1}{1,2}(2,:));
UpCI = transpose(out.predlims{1,1}{1,2}(3,:));

T = table(DownCI,Mean,UpCI, 'VariableNames', { 'DownCI', 'Mean','UpCI'} );
writetable(T, strcat(City,'_fitting_ContactTracing_narrow.csv'));



DHTr2=1-sum((data.ydata(:,2)-transpose(out.obslims{1,1}{1,1}(2,:))).^2)/ sum((data.ydata(:,2)-mean(data.ydata(:,2))).^2)
DHCr2=1-sum((data.ydata(:,3)-transpose(out.obslims{1,1}{1,2}(2,:))).^2)/ sum((data.ydata(:,3)-mean(data.ydata(:,3))).^2)


DHTr2=1-sum((data.ydata(:,2)-transpose(out.predlims{1,1}{1,1}(2,:))).^2)/sum((data.ydata(:,2)-mean(data.ydata(:,2))).^2)
DHCr2=1-sum((data.ydata(:,3)-transpose(out.predlims{1,1}{1,2}(2,:))).^2)/sum((data.ydata(:,3)-mean(data.ydata(:,3))).^2)




