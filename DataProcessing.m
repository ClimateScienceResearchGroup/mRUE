%%% KBS flux data processing
%%% Using non-gapfileld data from REddyProc

clear all
close all

%%% Read in T1 site data
T1 = dlmread("D:\desktop\T1\output.txt",'\t',2,0);
T1(T1 == -9999) = NaN;
T1(T1 == -10000) = NaN;

%%% Data Selection


NEE=T1(:,15);
NEE_qc=T1(:,17);
LE_f=T1(:,25);
Rg_f=T1(:,43);
Reco=T1(:,92);
GPP=T1(:,93);


%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T1 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);

writetable(DataOut,'D:\desktop\reprocessed data\T1_data_v2.csv')


%%% calc DAILY T1 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;


%%% Output DAILY T1 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T1_day_data_v2.csv')



%%% End Output T1 Data %%%
clear all
close all


%%% Read in T2 site data
T2 = dlmread("D:\desktop\T2\output.txt",'\t',2,0);
T2(T2 == -9999) = NaN;
T2(T2 == -10000) = NaN;

%%% Data Selection

NEE=T2(:,15);
NEE_qc=T2(:,17);
LE_f=T2(:,25);
Rg_f=T2(:,43);
Reco=T2(:,92);
GPP=T2(:,93);



%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T2 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T2_data_v2.csv')

%%% calc DAILY T2 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;


%%% Output DAILY T2 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T2_day_data_v2.csv')


%%% End Output T2 Data %%%
clear all
close all



%%% Read in T3 site data
T3 = dlmread("D:\desktop\T3\output.txt",'\t',2,0);
T3(T3 == -9999) = NaN;
T3(T3 == -10000) = NaN;

%%% Data Selection


NEE=T3(:,15);
NEE_qc=T3(:,17);
LE_f=T3(:,25);
Rg_f=T3(:,43);
Reco=T3(:,92);
GPP=T3(:,93);


%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T3 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T3_data_v2.csv')


%%% calc DAILY T3 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;


%%% Output DAILY T3 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T3_day_data_v2.csv')

%%% End Output T3 Data %%%
clear all
close all




%%% Read in T4 site data
T4 = dlmread("D:\desktop\T4\output.txt",'\t',2,0);
T4(T4 == -9999) = NaN;
T4(T4 == -10000) = NaN;

%%% Data Selection

NEE=T4(:,15);
NEE_qc=T4(:,17);
LE_f=T4(:,25);
Rg_f=T4(:,43);
Reco=T4(:,92);
GPP=T4(:,93);

%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;

%%% Output T4 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T4_data_v2.csv')


%%% calc DAILY T4 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;



%%% Output DAILY T4 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T4_day_data_v2.csv')


%%% End Output T4 Data %%%
clear all
close all




%%% Read in T5 site data
T5 = dlmread("D:\desktop\T5\output.txt",'\t',2,0);
T5(T5 == -9999) = NaN;
T5(T5 == -10000) = NaN;

%%% Data Selection

NEE=T5(:,15);
NEE_qc=T5(:,17);
LE_f=T5(:,25);
Rg_f=T5(:,43);
Reco=T5(:,92);
GPP=T5(:,93);


%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T5 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T5_data_v2.csv')


%%% calc DAILY T5 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;



%%% Output DAILY T5 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T5_day_data_v2.csv')

%%% End Output T5 Data %%%
clear all
close all




%%% Read in T6 site data
T6 = dlmread("D:\desktop\T6\output.txt",'\t',2,0);
T6(T6 == -9999) = NaN;
T6(T6 == -10000) = NaN;

%%% Data Selection

NEE=T6(:,15);
NEE_qc=T6(:,17);
LE_f=T6(:,25);
Rg_f=T6(:,43);
Reco=T6(:,92);
GPP=T6(:,93);


%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T6 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T6_data_v2.csv')


%%% calc DAILY T6 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;



%%% Output DAILY T6 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T6_day_data_v2.csv')

%%% End Output T6 Data %%%
clear all
close all



%%% Read in T7 site data
T7 = dlmread("D:\desktop\T7\output.txt",'\t',2,0);
T7(T7 == -9999) = NaN;
T7(T7 == -10000) = NaN;

%%% Data Selection

NEE=T7(:,15);
NEE_qc=T7(:,17);
LE_f=T7(:,25);
Rg_f=T7(:,43);
Reco=T7(:,92);
GPP=T7(:,93);


%%%% selecting only data were NEE obs was made (ie, no gap filled data)
gapfilled_loc=NEE_qc>0;
NEE(gapfilled_loc) = NaN;
LE_f(gapfilled_loc) = NaN;
Rg_f(gapfilled_loc) = NaN;
Reco(gapfilled_loc) = NaN;
GPP(gapfilled_loc) = NaN;


%%% Output T7 Site Data %%%

DataOut = table(...
NEE,...
LE_f,...
Rg_f,...
Reco,...
GPP);


writetable(DataOut,'D:\desktop\reprocessed data\T7_data_v2.csv')


%%% calc DAILY T7 Site Data %%%

n=48;
gapfilled_loc_day=arrayfun(@(i) nanmean(gapfilled_loc(i:i+n-1)),1:n:length(gapfilled_loc)-n+1)';
NEE_day=arrayfun(@(i) nanmean(NEE(i:i+n-1)),1:n:length(NEE)-n+1)';
LE_f_day=arrayfun(@(i) nanmean(LE_f(i:i+n-1)),1:n:length(LE_f)-n+1)';
Rg_f_day=arrayfun(@(i) nanmean(Rg_f(i:i+n-1)),1:n:length(Rg_f)-n+1)';
Reco_day=arrayfun(@(i) nanmean(Reco(i:i+n-1)),1:n:length(Reco)-n+1)';
GPP_day=arrayfun(@(i) nanmean(GPP(i:i+n-1)),1:n:length(GPP)-n+1)';


%daily gapfilled threshold
dailypercentgapfilled_loc=gapfilled_loc_day>0.6;
%sum(dailypercentgapfilled_loc)/length(dailypercentgapfilled_loc)*100

NEE_day(dailypercentgapfilled_loc) = NaN;
LE_f_day(dailypercentgapfilled_loc) = NaN;
Rg_f_day(dailypercentgapfilled_loc) = NaN;
Reco_day(dailypercentgapfilled_loc) = NaN;
GPP_day(dailypercentgapfilled_loc) = NaN;



%%% Output DAILY T7 Site Data %%%

DataOut = table(...
NEE_day,...
LE_f_day,...
Rg_f_day,...
Reco_day,...
GPP_day);

writetable(DataOut,'D:\desktop\reprocessed data\T7_day_data_v2.csv')


%%% End Output T7 Data %%%
clear all
close all













