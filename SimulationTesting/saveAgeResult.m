function [tsTotalMean,tsTotalUp,tsTotalDown] = saveAgeResult(y1,y2,y3,y4,interval,lag,iType,cityName,xdata)
    %%--------type-----------
    %1:dailycases
    %2:quarantined
    %3:hospitalization (current)
    %4:hospitalization considering primaty disease(current)
    %5:ICU (current)
    %6:ICU considering primaty disease(current)
    %7:hospitalization (cumulative)
    %8:hospitalization considering primaty disease(cumulative)
    %9:ICU (cumulative)
    %10:ICU considering primaty disease (cumulative)
    %-----------------------------

    typeSet = {'dailycases','quarantined','hospitalization','hospitalization_dis','icu','icu_dis','cum_hospitalization','cum_hospitalization_dis_extra','cum_icu','cum_icu_dis_extra'};
    [~,nDay,nCity] = size(y1); 
    lims = [0.025,0.5,0.975];

    % confidence interval
    CI1 = zeros(3,nDay,nCity);
    CI2 = zeros(3,nDay,nCity);
    CI3 = zeros(3,nDay,nCity);
    CI4 = zeros(3,nDay,nCity);
    for iCity=1:nCity
        CI1(:,:,iCity) = plims(y1(:,:,iCity),lims);
        CI2(:,:,iCity) = plims(y2(:,:,iCity),lims);
        CI3(:,:,iCity) = plims(y3(:,:,iCity),lims);
        CI4(:,:,iCity) = plims(y4(:,:,iCity),lims);
    end

    tsTotalMean = reshape(CI1(2,:,:),nDay,nCity)+reshape(CI2(2,:,:),nDay,nCity)+...
        reshape(CI3(2,:,:),nDay,nCity)+reshape(CI4(2,:,:),nDay,nCity);
    tsTotalUp = reshape(CI1(3,:,:),nDay,nCity)+reshape(CI2(3,:,:),nDay,nCity)+...
        reshape(CI3(3,:,:),nDay,nCity)+reshape(CI4(3,:,:),nDay,nCity);
    tsTotalDown = reshape(CI1(1,:,:),nDay,nCity)+reshape(CI2(1,:,:),nDay,nCity)+...
        reshape(CI3(1,:,:),nDay,nCity)+reshape(CI4(1,:,:),nDay,nCity);

    %% save result
    %% considering primary disease
    if (iType == 8) || (iType ==10)
        for iDis = 1:8
            if iType == 8
                OR = xdata.hosOR(iDis,:,:);
            else
                OR = xdata.icuOR(iDis,:,:);
            end
            CIDis1 = CI1.*(reshape(OR(:,1,:),1,1,nCity)-1);
            CIDis2 = CI2.*(reshape(OR(:,2,:),1,1,nCity)-1);
            CIDis3 = CI3.*(reshape(OR(:,3,:),1,1,nCity)-1);
            CIDis4 = CI4.*(reshape(OR(:,4,:),1,1,nCity)-1);
            tMean1 = array2table(reshape(CIDis1(2,:,:),nDay,nCity),'VariableNames',cityName);
            tMean2 = array2table(reshape(CIDis2(2,:,:),nDay,nCity),'VariableNames',cityName);
            tMean3 = array2table(reshape(CIDis3(2,:,:),nDay,nCity),'VariableNames',cityName);
            tMean4 = array2table(reshape(CIDis4(2,:,:),nDay,nCity),'VariableNames',cityName);
            tMean = array2table(reshape(CIDis1(2,:,:),nDay,nCity)+reshape(CIDis2(2,:,:),nDay,nCity)+...
                reshape(CIDis3(2,:,:),nDay,nCity)+reshape(CIDis4(2,:,:),nDay,nCity),'VariableNames',cityName);


            tUp1 = array2table(reshape(CIDis1(3,:,:),nDay,nCity),'VariableNames',cityName);
            tUp2 = array2table(reshape(CIDis2(3,:,:),nDay,nCity),'VariableNames',cityName);
            tUp3 = array2table(reshape(CIDis3(3,:,:),nDay,nCity),'VariableNames',cityName);
            tUp4 = array2table(reshape(CIDis4(3,:,:),nDay,nCity),'VariableNames',cityName);
            tUp = array2table(reshape(CIDis1(3,:,:),nDay,nCity)+reshape(CIDis2(3,:,:),nDay,nCity)+...
                reshape(CIDis3(3,:,:),nDay,nCity)+reshape(CIDis4(3,:,:),nDay,nCity),'VariableNames',cityName);

            tDown1 = array2table(reshape(CIDis1(1,:,:),nDay,nCity),'VariableNames',cityName);
            tDown2 = array2table(reshape(CIDis2(1,:,:),nDay,nCity),'VariableNames',cityName);
            tDown3 = array2table(reshape(CIDis3(1,:,:),nDay,nCity),'VariableNames',cityName);
            tDown4 = array2table(reshape(CIDis4(1,:,:),nDay,nCity),'VariableNames',cityName);
            tDown = array2table(reshape(CIDis1(1,:,:),nDay,nCity)+reshape(CIDis2(1,:,:),nDay,nCity)+...
                reshape(CIDis3(1,:,:),nDay,nCity)+reshape(CIDis4(1,:,:),nDay,nCity),'VariableNames',cityName);
        
            writetable(tMean1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_mean_",num2str(iDis),".xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tMean2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_mean_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tMean3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_mean_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tMean4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_mean_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");

            writetable(tUp1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_up_",num2str(iDis),".xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tUp2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_up_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tUp3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_up_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tUp4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_up_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");

            writetable(tDown1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_down_",num2str(iDis),".xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tDown2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_down_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tDown3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_down_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tDown4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_down_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");

            writetable(tMean, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_mean_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tUp, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_up_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
            writetable(tDown, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_down_",num2str(iDis),".xlsx"), ...
                    'Sheet',num2str(lag),"WriteMode","overwritesheet");
        end
    else
        tMean1 = array2table(reshape(CI1(2,:,:),nDay,nCity),'VariableNames',cityName);
        tMean2 = array2table(reshape(CI2(2,:,:),nDay,nCity),'VariableNames',cityName);
        tMean3 = array2table(reshape(CI3(2,:,:),nDay,nCity),'VariableNames',cityName);
        tMean4 = array2table(reshape(CI4(2,:,:),nDay,nCity),'VariableNames',cityName);
        tMean = array2table(reshape(CI1(2,:,:),nDay,nCity)+reshape(CI2(2,:,:),nDay,nCity)+...
            reshape(CI3(2,:,:),nDay,nCity)+reshape(CI4(2,:,:),nDay,nCity),'VariableNames',cityName);
        writetable(tMean1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_mean.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tMean2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_mean.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tMean3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_mean.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tMean4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_mean.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tMean, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_mean.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");

        tUp1 = array2table(reshape(CI1(3,:,:),nDay,nCity),'VariableNames',cityName);
        tUp2 = array2table(reshape(CI2(3,:,:),nDay,nCity),'VariableNames',cityName);
        tUp3 = array2table(reshape(CI3(3,:,:),nDay,nCity),'VariableNames',cityName);
        tUp4 = array2table(reshape(CI4(3,:,:),nDay,nCity),'VariableNames',cityName);
        tUp = array2table(reshape(CI1(3,:,:),nDay,nCity)+reshape(CI2(3,:,:),nDay,nCity)+...
            reshape(CI3(3,:,:),nDay,nCity)+reshape(CI4(3,:,:),nDay,nCity),'VariableNames',cityName);
        writetable(tUp1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_up.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tUp2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_up.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tUp3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_up.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tUp4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_up.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tUp, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_up.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");

        tDown1 = array2table(reshape(CI1(1,:,:),nDay,nCity),'VariableNames',cityName);
        tDown2 = array2table(reshape(CI2(1,:,:),nDay,nCity),'VariableNames',cityName);
        tDown3 = array2table(reshape(CI3(1,:,:),nDay,nCity),'VariableNames',cityName);
        tDown4 = array2table(reshape(CI4(1,:,:),nDay,nCity),'VariableNames',cityName);
        tDown = array2table(reshape(CI1(1,:,:),nDay,nCity)+reshape(CI2(1,:,:),nDay,nCity)+...
            reshape(CI3(1,:,:),nDay,nCity)+reshape(CI4(1,:,:),nDay,nCity),'VariableNames',cityName);
        writetable(tDown1, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_down.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tDown2, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_20_40_down.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tDown3, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_40_60_down.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tDown4, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_60_down.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
        writetable(tDown, strcat(".\simu_result\interval_",num2str(interval),"_",typeSet(iType),"_down.xlsx"), ...
                'Sheet',num2str(lag),"WriteMode","overwritesheet");
    end

end