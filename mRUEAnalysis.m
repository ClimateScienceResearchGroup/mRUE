%%% KBS flux manuscript analysis
%%% Multiple Resource Use Efficiency

clear all
close all

% read in all site data NOW IN TABLE 


T1_day = readtable('D:\desktop\reprocessed data\T1_day_data_v2.csv','ReadVariableNames',true);
T2_day = readtable('D:\desktop\reprocessed data\T2_day_data_v2.csv','ReadVariableNames',true);
T3_day = readtable('D:\desktop\reprocessed data\T3_day_data_v2.csv','ReadVariableNames',true);
T4_day = readtable('D:\desktop\reprocessed data\T4_day_data_v2.csv','ReadVariableNames',true);
T5_day = readtable('D:\desktop\reprocessed data\T5_day_data_v2.csv','ReadVariableNames',true);
T6_day = readtable('D:\desktop\reprocessed data\T6_day_data_v2.csv','ReadVariableNames',true);
T7_day = readtable('D:\desktop\reprocessed data\T7_day_data_v2.csv','ReadVariableNames',true);




% calculate water use effiniency NEP/ET
T1_day.WUE=T1_day.GPP_day./T1_day.LE_f_day;
T2_day.WUE=T2_day.GPP_day./T2_day.LE_f_day;
T3_day.WUE=T3_day.GPP_day./T3_day.LE_f_day;
T4_day.WUE=T4_day.GPP_day./T4_day.LE_f_day;
T5_day.WUE=T5_day.GPP_day./T5_day.LE_f_day;
T6_day.WUE=T6_day.GPP_day./T6_day.LE_f_day;
T7_day.WUE=T7_day.GPP_day./T7_day.LE_f_day;




% calculate light use efficiency NEP/(Rn)
T1_day.LUE=T1_day.GPP_day./T1_day.Rg_f_day;
T2_day.LUE=T2_day.GPP_day./T2_day.Rg_f_day;
T3_day.LUE=T3_day.GPP_day./T3_day.Rg_f_day;
T4_day.LUE=T4_day.GPP_day./T4_day.Rg_f_day;
T5_day.LUE=T5_day.GPP_day./T5_day.Rg_f_day;
T6_day.LUE=T6_day.GPP_day./T6_day.Rg_f_day;
T7_day.LUE=T7_day.GPP_day./T7_day.Rg_f_day;


% calculate carbon use efficiency NEP/Reco
T1_day.CUE=T1_day.GPP_day./T1_day.Reco_day;
T2_day.CUE=T2_day.GPP_day./T2_day.Reco_day;
T3_day.CUE=T3_day.GPP_day./T3_day.Reco_day;
T4_day.CUE=T4_day.GPP_day./T4_day.Reco_day;
T5_day.CUE=T5_day.GPP_day./T5_day.Reco_day;
T6_day.CUE=T6_day.GPP_day./T6_day.Reco_day;
T7_day.CUE=T7_day.GPP_day./T7_day.Reco_day;



%Figure 1
%%%% COLOR GUIDE
%Water (40, 109, 238) [0.1569 0.4275 0.9333]
%Light (253, 240, 90) [0.9922 0.9412 0.3529]
%Carbon (238, 70, 81) [0.9333 0.2745 0.3176]

subplot(3,1,1)
hold on

plot(movmean(T1_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T2_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T3_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T4_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T5_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T6_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(T7_day.WUE,7,'omitnan'),'color',[0.1569 0.4275 0.9333])
plot(movmean(nanmean([T1_day.WUE,T2_day.WUE,T3_day.WUE,T4_day.WUE,T5_day.WUE,T6_day.WUE,T7_day.WUE],2),7,'omitnan'),'linewidth',.5,'color',[0 0 0])
box on
axis([0 2922 -.08 .28])
ylabel({'WUE','[?mol m^-^2 s^-^1 / W m^-^2]'})
xticks([0 1*365 2*365 3*365 4*365 5*365 6*365 7*365])
xticklabels({'2009','2010','2011','2012','2013','2014','2015','2016'})
text(2922-160,0.25,'(a)')


subplot(3,1,2)
hold on
plot(movmean(T1_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T2_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T3_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T4_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T5_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T6_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(T7_day.LUE,7,'omitnan'),'color',[0.9922 0.9412 0.3529])
plot(movmean(nanmean([T1_day.LUE,T2_day.LUE,T3_day.LUE,T4_day.LUE,T5_day.LUE,T6_day.LUE,T7_day.LUE],2),7,'omitnan'),'linewidth',.5,'color',[0 0 0])
box on
axis([0 2922 -.04 .14])
ylabel({'LUE','[?mol m^-^2 s^-^1 / W m^-^2]'})
xticks([0 1*365 2*365 3*365 4*365 5*365 6*365 7*365])
xticklabels({'2009','2010','2011','2012','2013','2014','2015','2016'})
text(2922-70,0.12,'(b)')


subplot(3,1,3)
hold on
plot(movmean(T1_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T2_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T3_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T4_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T5_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T6_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(T7_day.CUE,7,'omitnan'),'color',[0.9333 0.2745 0.3176])
plot(movmean(nanmean([T1_day.CUE,T2_day.CUE,T3_day.CUE,T4_day.CUE,T5_day.CUE,T6_day.CUE,T7_day.CUE],2),7,'omitnan'),'linewidth',.5,'color',[0 0 0])
box on
axis([0 2922 -2.5 5.5])
ylabel({'CUE','[?mol m^-^2 s^-^1 / ?mol m^-^2 s^-^1]'})
xticks([0 1*365 2*365 3*365 4*365 5*365 6*365 7*365])
xticklabels({'2009','2010','2011','2012','2013','2014','2015','2016'})
text(2922-70,4.9,'(c)')







% Rescaling mRUE's, -5 to 5
T1_day.LUEs = scaledata(T1_day.LUE,-5,5);
T2_day.LUEs = scaledata(T2_day.LUE,-5,5);
T3_day.LUEs = scaledata(T3_day.LUE,-5,5);
T4_day.LUEs = scaledata(T4_day.LUE,-5,5);
T5_day.LUEs = scaledata(T5_day.LUE,-5,5);
T6_day.LUEs = scaledata(T6_day.LUE,-5,5);
T7_day.LUEs = scaledata(T7_day.LUE,-5,5);

T1_day.WUEs = scaledata(T1_day.WUE,-5,5);
T2_day.WUEs = scaledata(T2_day.WUE,-5,5);
T3_day.WUEs = scaledata(T3_day.WUE,-5,5);
T4_day.WUEs = scaledata(T4_day.WUE,-5,5);
T5_day.WUEs = scaledata(T5_day.WUE,-5,5);
T6_day.WUEs = scaledata(T6_day.WUE,-5,5);
T7_day.WUEs = scaledata(T7_day.WUE,-5,5);

T1_day.CUEs = scaledata(T1_day.CUE,-5,5);
T2_day.CUEs = scaledata(T2_day.CUE,-5,5);
T3_day.CUEs = scaledata(T3_day.CUE,-5,5);
T4_day.CUEs = scaledata(T4_day.CUE,-5,5);
T5_day.CUEs = scaledata(T5_day.CUE,-5,5);
T6_day.CUEs = scaledata(T6_day.CUE,-5,5);
T7_day.CUEs = scaledata(T7_day.CUE,-5,5);



% ONE RUE TO RULE THEM ALL
% output mRUE for R code

Year=cat(2,ones(1,365).*2009,ones(1,365).*2010,ones(1,365).*2011,ones(1,366).*2012,ones(1,365).*2013,ones(1,365).*2014,ones(1,365).*2015,ones(1,366).*2016)';
DOY=cat(2,[1:365],[1:365],[1:365],[1:366],[1:365],[1:365],[1:365],[1:366])';

DataOut = table(...
Year,...
DOY,...
T1_day.LUEs,...
T2_day.LUEs,...
T3_day.LUEs,...
T4_day.LUEs,...
T5_day.LUEs,...
T6_day.LUEs,...
T7_day.LUEs,...
T1_day.WUEs,...
T2_day.WUEs,...
T3_day.WUEs,...
T4_day.WUEs,...
T5_day.WUEs,...
T6_day.WUEs,...
T7_day.WUEs,...
T1_day.CUEs,...
T2_day.CUEs,...
T3_day.CUEs,...
T4_day.CUEs,...
T5_day.CUEs,...
T6_day.CUEs,...
T7_day.CUEs,...
-T1_day.NEE_day,...
-T2_day.NEE_day,...
-T3_day.NEE_day,...
-T4_day.NEE_day,...
-T5_day.NEE_day,...
-T6_day.NEE_day,...
-T7_day.NEE_day,'VariableNames',{'Year' 'DOY'...
'T1_day_LUEs' 'T2_day_LUEs' 'T3_day_LUEs' 'T4_day_LUEs' 'T5_day_LUEs' 'T6_day_LUEs' 'T7_day_LUEs'...
'T1_day_WUEs' 'T2_day_WUEs' 'T3_day_WUEs' 'T4_day_WUEs' 'T5_day_WUEs' 'T6_day_WUEs' 'T7_day_WUEs'...
'T1_day_CUEs' 'T2_day_CUEs' 'T3_day_CUEs' 'T4_day_CUEs' 'T5_day_CUEs' 'T6_day_CUEs' 'T7_day_CUEs'...
'T1_day_NEP' 'T2_day_NEP' 'T3_day_NEP' 'T4_day_NEP' 'T5_day_NEP' 'T6_day_NEP' 'T7_day_NEP',});


writetable(DataOut,'D:\desktop\reprocessed data\mRUE_data_v2.csv')




