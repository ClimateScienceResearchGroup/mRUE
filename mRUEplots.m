clear all
close all
%%%% COLOR GUIDE
%Water (40, 109, 238) [0.1569 0.4275 0.9333]
%Light (253, 240, 90) [0.9922 0.9412 0.3529]
%Carbon (238, 70, 81) [0.9333 0.2745 0.3176]
xtime=2009:2016;

input = csvread('D:\desktop\reprocessed data\all_loadingfactors_v2.1.dat',1,1);

all_loadingfactors=zeros(3,8,7);

all_loadingfactors(:,:,1) = input(1:3,1:8);
all_loadingfactors(:,:,2) = input(1:3,9:16);
all_loadingfactors(:,:,3) = input(1:3,17:24);
all_loadingfactors(:,:,4) = input(1:3,25:32);
all_loadingfactors(:,:,5) = input(1:3,33:40);
all_loadingfactors(:,:,6) = input(1:3,41:48);
all_loadingfactors(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\all_nan_count_v2.1.dat',1,1);

all_nan_count=zeros(3,8,7);

all_nan_count(:,:,1) = input(1:3,1:8);
all_nan_count(:,:,2) = input(1:3,9:16);
all_nan_count(:,:,3) = input(1:3,17:24);
all_nan_count(:,:,4) = input(1:3,25:32);
all_nan_count(:,:,5) = input(1:3,33:40);
all_nan_count(:,:,6) = input(1:3,41:48);
all_nan_count(:,:,7) = input(1:3,49:56);



%%%% thershold for full year, 70 (20% of 365)

all_loadingfactors(all_nan_count>70)=nan;
mean_all_loadingfactors = nanmean(all_loadingfactors,3);


%%%%% Growing Season data

clear input
input = csvread('D:\desktop\reprocessed data\all_loadingfactors_GS_v2.1.dat',1,1);

all_loadingfactors_GS=zeros(3,8,7);

all_loadingfactors_GS(:,:,1) = input(1:3,1:8);
all_loadingfactors_GS(:,:,2) = input(1:3,9:16);
all_loadingfactors_GS(:,:,3) = input(1:3,17:24);
all_loadingfactors_GS(:,:,4) = input(1:3,25:32);
all_loadingfactors_GS(:,:,5) = input(1:3,33:40);
all_loadingfactors_GS(:,:,6) = input(1:3,41:48);
all_loadingfactors_GS(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\all_nan_count_GS_v2.1.dat',1,1);

all_nan_count_GS=zeros(3,8,7);

all_nan_count_GS(:,:,1) = input(1:3,1:8);
all_nan_count_GS(:,:,2) = input(1:3,9:16);
all_nan_count_GS(:,:,3) = input(1:3,17:24);
all_nan_count_GS(:,:,4) = input(1:3,25:32);
all_nan_count_GS(:,:,5) = input(1:3,33:40);
all_nan_count_GS(:,:,6) = input(1:3,41:48);
all_nan_count_GS(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 295(20% of four months)

all_loadingfactors_GS(all_nan_count_GS>295)=nan;
mean_all_loadingfactors_GS = nanmean(all_loadingfactors_GS,3);



%%%%% 8 years data

clear input
input = csvread('D:\desktop\reprocessed data\all_loadingfactors_allyears_v2.dat',1,1);

all_loadingfactors_allyears=zeros(3,8,7);

all_loadingfactors_allyears(:,:,1) = input(1:3,1:8);
all_loadingfactors_allyears(:,:,2) = input(1:3,9:16);
all_loadingfactors_allyears(:,:,3) = input(1:3,17:24);
all_loadingfactors_allyears(:,:,4) = input(1:3,25:32);
all_loadingfactors_allyears(:,:,5) = input(1:3,33:40);
all_loadingfactors_allyears(:,:,6) = input(1:3,41:48);
all_loadingfactors_allyears(:,:,7) = input(1:3,49:56);





%%%%%% MONTH BY MONTH DATA
%%% Jan

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_jan_v2.1.dat',1,1);

all_loadingfactors_jan=zeros(3,8,7);

all_loadingfactors_jan(:,:,1) = input(1:3,1:8);
all_loadingfactors_jan(:,:,2) = input(1:3,9:16);
all_loadingfactors_jan(:,:,3) = input(1:3,17:24);
all_loadingfactors_jan(:,:,4) = input(1:3,25:32);
all_loadingfactors_jan(:,:,5) = input(1:3,33:40);
all_loadingfactors_jan(:,:,6) = input(1:3,41:48);
all_loadingfactors_jan(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_jan_v2.1.dat',1,1);

all_nan_count_jan=zeros(3,8,7);

all_nan_count_jan(:,:,1) = input(1:3,1:8);
all_nan_count_jan(:,:,2) = input(1:3,9:16);
all_nan_count_jan(:,:,3) = input(1:3,17:24);
all_nan_count_jan(:,:,4) = input(1:3,25:32);
all_nan_count_jan(:,:,5) = input(1:3,33:40);
all_nan_count_jan(:,:,6) = input(1:3,41:48);
all_nan_count_jan(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_jan(all_nan_count_jan>7)=nan;
mean_all_loadingfactors_jan = nanmean(all_loadingfactors_jan,3);




%%% feb

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_feb_v2.dat',1,1);

all_loadingfactors_feb=zeros(3,8,7);

all_loadingfactors_feb(:,:,1) = input(1:3,1:8);
all_loadingfactors_feb(:,:,2) = input(1:3,9:16);
all_loadingfactors_feb(:,:,3) = input(1:3,17:24);
all_loadingfactors_feb(:,:,4) = input(1:3,25:32);
all_loadingfactors_feb(:,:,5) = input(1:3,33:40);
all_loadingfactors_feb(:,:,6) = input(1:3,41:48);
all_loadingfactors_feb(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_feb_v2.dat',1,1);

all_nan_count_feb=zeros(3,8,7);

all_nan_count_feb(:,:,1) = input(1:3,1:8);
all_nan_count_feb(:,:,2) = input(1:3,9:16);
all_nan_count_feb(:,:,3) = input(1:3,17:24);
all_nan_count_feb(:,:,4) = input(1:3,25:32);
all_nan_count_feb(:,:,5) = input(1:3,33:40);
all_nan_count_feb(:,:,6) = input(1:3,41:48);
all_nan_count_feb(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_feb(all_nan_count_feb>7)=nan;
mean_all_loadingfactors_feb = nanmean(all_loadingfactors_feb,3);



%%% mar

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_mar_v2.dat',1,1);

all_loadingfactors_mar=zeros(3,8,7);

all_loadingfactors_mar(:,:,1) = input(1:3,1:8);
all_loadingfactors_mar(:,:,2) = input(1:3,9:16);
all_loadingfactors_mar(:,:,3) = input(1:3,17:24);
all_loadingfactors_mar(:,:,4) = input(1:3,25:32);
all_loadingfactors_mar(:,:,5) = input(1:3,33:40);
all_loadingfactors_mar(:,:,6) = input(1:3,41:48);
all_loadingfactors_mar(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_mar_v2.dat',1,1);

all_nan_count_mar=zeros(3,8,7);

all_nan_count_mar(:,:,1) = input(1:3,1:8);
all_nan_count_mar(:,:,2) = input(1:3,9:16);
all_nan_count_mar(:,:,3) = input(1:3,17:24);
all_nan_count_mar(:,:,4) = input(1:3,25:32);
all_nan_count_mar(:,:,5) = input(1:3,33:40);
all_nan_count_mar(:,:,6) = input(1:3,41:48);
all_nan_count_mar(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_mar(all_nan_count_mar>7)=nan;
mean_all_loadingfactors_mar = nanmean(all_loadingfactors_mar,3);




%%% apr

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_apr_v2.dat',1,1);

all_loadingfactors_apr=zeros(3,8,7);

all_loadingfactors_apr(:,:,1) = input(1:3,1:8);
all_loadingfactors_apr(:,:,2) = input(1:3,9:16);
all_loadingfactors_apr(:,:,3) = input(1:3,17:24);
all_loadingfactors_apr(:,:,4) = input(1:3,25:32);
all_loadingfactors_apr(:,:,5) = input(1:3,33:40);
all_loadingfactors_apr(:,:,6) = input(1:3,41:48);
all_loadingfactors_apr(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_apr_v2.dat',1,1);

all_nan_count_apr=zeros(3,8,7);

all_nan_count_apr(:,:,1) = input(1:3,1:8);
all_nan_count_apr(:,:,2) = input(1:3,9:16);
all_nan_count_apr(:,:,3) = input(1:3,17:24);
all_nan_count_apr(:,:,4) = input(1:3,25:32);
all_nan_count_apr(:,:,5) = input(1:3,33:40);
all_nan_count_apr(:,:,6) = input(1:3,41:48);
all_nan_count_apr(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_apr(all_nan_count_apr>7)=nan;
mean_all_loadingfactors_apr = nanmean(all_loadingfactors_apr,3);



%%% may

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_may_v2.dat',1,1);

all_loadingfactors_may=zeros(3,8,7);

all_loadingfactors_may(:,:,1) = input(1:3,1:8);
all_loadingfactors_may(:,:,2) = input(1:3,9:16);
all_loadingfactors_may(:,:,3) = input(1:3,17:24);
all_loadingfactors_may(:,:,4) = input(1:3,25:32);
all_loadingfactors_may(:,:,5) = input(1:3,33:40);
all_loadingfactors_may(:,:,6) = input(1:3,41:48);
all_loadingfactors_may(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_may_v2.dat',1,1);

all_nan_count_may=zeros(3,8,7);

all_nan_count_may(:,:,1) = input(1:3,1:8);
all_nan_count_may(:,:,2) = input(1:3,9:16);
all_nan_count_may(:,:,3) = input(1:3,17:24);
all_nan_count_may(:,:,4) = input(1:3,25:32);
all_nan_count_may(:,:,5) = input(1:3,33:40);
all_nan_count_may(:,:,6) = input(1:3,41:48);
all_nan_count_may(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_may(all_nan_count_may>7)=nan;
mean_all_loadingfactors_may = nanmean(all_loadingfactors_may,3);






%%% Jun

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_jun_v2.dat',1,1);

all_loadingfactors_jun=zeros(3,8,7);

all_loadingfactors_jun(:,:,1) = input(1:3,1:8);
all_loadingfactors_jun(:,:,2) = input(1:3,9:16);
all_loadingfactors_jun(:,:,3) = input(1:3,17:24);
all_loadingfactors_jun(:,:,4) = input(1:3,25:32);
all_loadingfactors_jun(:,:,5) = input(1:3,33:40);
all_loadingfactors_jun(:,:,6) = input(1:3,41:48);
all_loadingfactors_jun(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_jun_v2.dat',1,1);

all_nan_count_jun=zeros(3,8,7);

all_nan_count_jun(:,:,1) = input(1:3,1:8);
all_nan_count_jun(:,:,2) = input(1:3,9:16);
all_nan_count_jun(:,:,3) = input(1:3,17:24);
all_nan_count_jun(:,:,4) = input(1:3,25:32);
all_nan_count_jun(:,:,5) = input(1:3,33:40);
all_nan_count_jun(:,:,6) = input(1:3,41:48);
all_nan_count_jun(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_jun(all_nan_count_jun>7)=nan;
mean_all_loadingfactors_jun = nanmean(all_loadingfactors_jun,3);



%%% Jul
clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_jul_v2.dat',1,1);

all_loadingfactors_jul=zeros(3,8,7);

all_loadingfactors_jul(:,:,1) = input(1:3,1:8);
all_loadingfactors_jul(:,:,2) = input(1:3,9:16);
all_loadingfactors_jul(:,:,3) = input(1:3,17:24);
all_loadingfactors_jul(:,:,4) = input(1:3,25:32);
all_loadingfactors_jul(:,:,5) = input(1:3,33:40);
all_loadingfactors_jul(:,:,6) = input(1:3,41:48);
all_loadingfactors_jul(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_jul_v2.dat',1,1);

all_nan_count_jul=zeros(3,8,7);

all_nan_count_jul(:,:,1) = input(1:3,1:8);
all_nan_count_jul(:,:,2) = input(1:3,9:16);
all_nan_count_jul(:,:,3) = input(1:3,17:24);
all_nan_count_jul(:,:,4) = input(1:3,25:32);
all_nan_count_jul(:,:,5) = input(1:3,33:40);
all_nan_count_jul(:,:,6) = input(1:3,41:48);
all_nan_count_jul(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_jul(all_nan_count_jul>7)=nan;
mean_all_loadingfactors_jul = nanmean(all_loadingfactors_jul,3);



%%% aug
clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_aug_v2.dat',1,1);

all_loadingfactors_aug=zeros(3,8,7);

all_loadingfactors_aug(:,:,1) = input(1:3,1:8);
all_loadingfactors_aug(:,:,2) = input(1:3,9:16);
all_loadingfactors_aug(:,:,3) = input(1:3,17:24);
all_loadingfactors_aug(:,:,4) = input(1:3,25:32);
all_loadingfactors_aug(:,:,5) = input(1:3,33:40);
all_loadingfactors_aug(:,:,6) = input(1:3,41:48);
all_loadingfactors_aug(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_aug_v2.dat',1,1);

all_nan_count_aug=zeros(3,8,7);

all_nan_count_aug(:,:,1) = input(1:3,1:8);
all_nan_count_aug(:,:,2) = input(1:3,9:16);
all_nan_count_aug(:,:,3) = input(1:3,17:24);
all_nan_count_aug(:,:,4) = input(1:3,25:32);
all_nan_count_aug(:,:,5) = input(1:3,33:40);
all_nan_count_aug(:,:,6) = input(1:3,41:48);
all_nan_count_aug(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_aug(all_nan_count_aug>7)=nan;
mean_all_loadingfactors_aug = nanmean(all_loadingfactors_aug,3);




%%% sep

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_sep_v2.dat',1,1);

all_loadingfactors_sep=zeros(3,8,7);

all_loadingfactors_sep(:,:,1) = input(1:3,1:8);
all_loadingfactors_sep(:,:,2) = input(1:3,9:16);
all_loadingfactors_sep(:,:,3) = input(1:3,17:24);
all_loadingfactors_sep(:,:,4) = input(1:3,25:32);
all_loadingfactors_sep(:,:,5) = input(1:3,33:40);
all_loadingfactors_sep(:,:,6) = input(1:3,41:48);
all_loadingfactors_sep(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_sep_v2.dat',1,1);

all_nan_count_sep=zeros(3,8,7);

all_nan_count_sep(:,:,1) = input(1:3,1:8);
all_nan_count_sep(:,:,2) = input(1:3,9:16);
all_nan_count_sep(:,:,3) = input(1:3,17:24);
all_nan_count_sep(:,:,4) = input(1:3,25:32);
all_nan_count_sep(:,:,5) = input(1:3,33:40);
all_nan_count_sep(:,:,6) = input(1:3,41:48);
all_nan_count_sep(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_sep(all_nan_count_sep>7)=nan;
mean_all_loadingfactors_sep = nanmean(all_loadingfactors_sep,3);



%%% oct

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_oct_v2.1.dat',1,1);

all_loadingfactors_oct=zeros(3,8,7);

all_loadingfactors_oct(:,:,1) = input(1:3,1:8);
all_loadingfactors_oct(:,:,2) = input(1:3,9:16);
all_loadingfactors_oct(:,:,3) = input(1:3,17:24);
all_loadingfactors_oct(:,:,4) = input(1:3,25:32);
all_loadingfactors_oct(:,:,5) = input(1:3,33:40);
all_loadingfactors_oct(:,:,6) = input(1:3,41:48);
all_loadingfactors_oct(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_oct_v2.1.dat',1,1);

all_nan_count_oct=zeros(3,8,7);

all_nan_count_oct(:,:,1) = input(1:3,1:8);
all_nan_count_oct(:,:,2) = input(1:3,9:16);
all_nan_count_oct(:,:,3) = input(1:3,17:24);
all_nan_count_oct(:,:,4) = input(1:3,25:32);
all_nan_count_oct(:,:,5) = input(1:3,33:40);
all_nan_count_oct(:,:,6) = input(1:3,41:48);
all_nan_count_oct(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_oct(all_nan_count_oct>7)=nan;
mean_all_loadingfactors_oct = nanmean(all_loadingfactors_oct,3);




%%% nov

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_nov_v2.dat',1,1);

all_loadingfactors_nov=zeros(3,8,7);

all_loadingfactors_nov(:,:,1) = input(1:3,1:8);
all_loadingfactors_nov(:,:,2) = input(1:3,9:16);
all_loadingfactors_nov(:,:,3) = input(1:3,17:24);
all_loadingfactors_nov(:,:,4) = input(1:3,25:32);
all_loadingfactors_nov(:,:,5) = input(1:3,33:40);
all_loadingfactors_nov(:,:,6) = input(1:3,41:48);
all_loadingfactors_nov(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_nov_v2.dat',1,1);

all_nan_count_nov=zeros(3,8,7);

all_nan_count_nov(:,:,1) = input(1:3,1:8);
all_nan_count_nov(:,:,2) = input(1:3,9:16);
all_nan_count_nov(:,:,3) = input(1:3,17:24);
all_nan_count_nov(:,:,4) = input(1:3,25:32);
all_nan_count_nov(:,:,5) = input(1:3,33:40);
all_nan_count_nov(:,:,6) = input(1:3,41:48);
all_nan_count_nov(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_nov(all_nan_count_nov>7)=nan;
mean_all_loadingfactors_nov = nanmean(all_loadingfactors_nov,3);




%%% dec

clear input
input = csvread('D:\desktop\reprocessed data\full year\all_loadingfactors_dec_v2.1.dat',1,1);

all_loadingfactors_dec=zeros(3,8,7);

all_loadingfactors_dec(:,:,1) = input(1:3,1:8);
all_loadingfactors_dec(:,:,2) = input(1:3,9:16);
all_loadingfactors_dec(:,:,3) = input(1:3,17:24);
all_loadingfactors_dec(:,:,4) = input(1:3,25:32);
all_loadingfactors_dec(:,:,5) = input(1:3,33:40);
all_loadingfactors_dec(:,:,6) = input(1:3,41:48);
all_loadingfactors_dec(:,:,7) = input(1:3,49:56);


clear input
input = csvread('D:\desktop\reprocessed data\full year\all_nan_count_dec_v2.1.dat',1,1);

all_nan_count_dec=zeros(3,8,7);

all_nan_count_dec(:,:,1) = input(1:3,1:8);
all_nan_count_dec(:,:,2) = input(1:3,9:16);
all_nan_count_dec(:,:,3) = input(1:3,17:24);
all_nan_count_dec(:,:,4) = input(1:3,25:32);
all_nan_count_dec(:,:,5) = input(1:3,33:40);
all_nan_count_dec(:,:,6) = input(1:3,41:48);
all_nan_count_dec(:,:,7) = input(1:3,49:56);


%%%% thershold for full year, 7(20% of one months)

all_loadingfactors_dec(all_nan_count_dec>7)=nan;
mean_all_loadingfactors_dec = nanmean(all_loadingfactors_dec,3);








mean_loadingfactor = nanmean(mean_all_loadingfactors,2);
%nanstd(mean_all_loadingfactors(1,:))
%nanstd(mean_all_loadingfactors(2,:))
%nanstd(mean_all_loadingfactors(3,:))


mean_loadingfactor_GS = nanmean(mean_all_loadingfactors_GS,2);

mean_loadingfactor_allyears = nanmean(all_loadingfactors_allyears,2);
mean_loadingfactor_allyears = nanmean(mean_loadingfactor_allyears,3);


mean_loadingfactor_jan = nanmean(mean_all_loadingfactors_jan,2);
mean_loadingfactor_feb = nanmean(mean_all_loadingfactors_feb,2);
mean_loadingfactor_mar = nanmean(mean_all_loadingfactors_mar,2);
mean_loadingfactor_apr = nanmean(mean_all_loadingfactors_apr,2);
mean_loadingfactor_may = nanmean(mean_all_loadingfactors_may,2);
mean_loadingfactor_jun = nanmean(mean_all_loadingfactors_jun,2);
mean_loadingfactor_jul = nanmean(mean_all_loadingfactors_jul,2);
mean_loadingfactor_aug = nanmean(mean_all_loadingfactors_aug,2);
mean_loadingfactor_sep = nanmean(mean_all_loadingfactors_sep,2);
mean_loadingfactor_oct = nanmean(mean_all_loadingfactors_oct,2);
mean_loadingfactor_nov = nanmean(mean_all_loadingfactors_nov,2);
mean_loadingfactor_dec = nanmean(mean_all_loadingfactors_dec,2);



ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')













a=[mean_loadingfactor_jan(1),mean_loadingfactor_feb(1),mean_loadingfactor_mar(1),mean_loadingfactor_apr(1),mean_loadingfactor_may(1),mean_loadingfactor_jun(1),mean_loadingfactor_jul(1),mean_loadingfactor_aug(1),mean_loadingfactor_sep(1),mean_loadingfactor_oct(1),mean_loadingfactor_nov(1),mean_loadingfactor_dec(1)];
b=[mean_loadingfactor_jan(2),mean_loadingfactor_feb(2),mean_loadingfactor_mar(2),mean_loadingfactor_apr(2),mean_loadingfactor_may(2),mean_loadingfactor_jun(2),mean_loadingfactor_jul(2),mean_loadingfactor_aug(2),mean_loadingfactor_sep(2),mean_loadingfactor_oct(2),mean_loadingfactor_nov(2),mean_loadingfactor_dec(2)];
c=abs([mean_loadingfactor_jan(3),mean_loadingfactor_feb(3),mean_loadingfactor_mar(3),mean_loadingfactor_apr(3),mean_loadingfactor_may(3),mean_loadingfactor_jun(3),mean_loadingfactor_jul(3),mean_loadingfactor_aug(3),mean_loadingfactor_sep(3),mean_loadingfactor_oct(3),mean_loadingfactor_nov(3),mean_loadingfactor_dec(3)]);

% TERN PLOT FIGURE FIVE


subplot(1,2,1)

ternplot(a,b,c,'k:o','majors', 5)
hold on
%ternplot(mean_loadingfactor_jan(1),mean_loadingfactor_jan(2),mean_loadingfactor_jan(3),'ko','MarkerSize',4)
% ternplot(mean_loadingfactor_feb(1),mean_loadingfactor_feb(2),mean_loadingfactor_feb(3),'ko','MarkerSize',4.5)
% ternplot(mean_loadingfactor_mar(1),mean_loadingfactor_mar(2),mean_loadingfactor_mar(3),'ko','MarkerSize',5)
% ternplot(mean_loadingfactor_apr(1),mean_loadingfactor_apr(2),mean_loadingfactor_apr(3),'ko','MarkerSize',5.5)
% ternplot(mean_loadingfactor_may(1),mean_loadingfactor_may(2),mean_loadingfactor_may(3),'ko','MarkerSize',6)
% ternplot(mean_loadingfactor_jun(1),mean_loadingfactor_jun(2),mean_loadingfactor_jun(3),'ko','MarkerSize',6.5)
% ternplot(mean_loadingfactor_jul(1),mean_loadingfactor_jul(2),mean_loadingfactor_jul(3),'ko','MarkerSize',7)
% ternplot(mean_loadingfactor_aug(1),mean_loadingfactor_aug(2),abs(mean_loadingfactor_aug(3)),'ko','MarkerSize',7.5)
% ternplot(mean_loadingfactor_sep(1),mean_loadingfactor_sep(2),abs(mean_loadingfactor_sep(3)),'ko','MarkerSize',8)
% ternplot(mean_loadingfactor_oct(1),mean_loadingfactor_oct(2),mean_loadingfactor_oct(3),'ko','MarkerSize',8.5)
% ternplot(mean_loadingfactor_nov(1),mean_loadingfactor_nov(2),mean_loadingfactor_nov(3),'ko','MarkerSize',9)
% ternplot(mean_loadingfactor_dec(1),mean_loadingfactor_dec(2),mean_loadingfactor_dec(3),'ko','MarkerSize',9.5)

text(.52,.33,'(Jan)')
text(.5,.23,'(Apr)')
text(.51,.1,'(Jul)')
text(.35,.42,'(Oct)')
text(.9,.9,'(a)')

ternlabel({'   ';'  ';'WUE Loading Factors [%]'},{'LUE Loading Factors [%]';'  ';'   '},{'CUE Loading Factors [%]';'  ';'   '})


subplot(1,2,2)

ternplot([mean_loadingfactor(1),mean_loadingfactor_GS(1),mean(a)],[mean_loadingfactor(2),mean_loadingfactor_GS(2),mean(b)],[mean_loadingfactor(3),mean_loadingfactor_GS(3),mean(c)],'k:','majors', 5)
hold on
ternplot(mean_loadingfactor(1),mean_loadingfactor(2),mean_loadingfactor(3),'ko','MarkerSize',9)
ternplot(mean_loadingfactor_GS(1),mean_loadingfactor_GS(2),mean_loadingfactor_GS(3),'ko','MarkerSize',6)
ternplot(mean(a),mean(b),mean(c),'ko','MarkerSize',3)
text(.52,.55,'(Annual)')
text(.58,.12,'(Growing Season)')
text(.45,.3,'(Month)')
text(.9,.9,'(b)')

ternlabel({'   ';'  ';'WUE Loading Factors [%]'},{'LUE Loading Factors [%]';'  ';'   '},{'CUE Loading Factors [%]';'  ';'   '})







%%%% CoorCoef work

coor=cat(1,(mean_all_loadingfactors(1,:)),(mean_all_loadingfactors(2,:)),(mean_all_loadingfactors(3,:)))'
coor_GS=cat(1,(mean_all_loadingfactors_GS(1,:)),(mean_all_loadingfactors_GS(2,:)),(mean_all_loadingfactors_GS(3,:)))'
coor_precip=cat(1,mean_all_loadingfactors(1,:)./gull_lake_precip,mean_all_loadingfactors(2,:)./gull_lake_precip,mean_all_loadingfactors(3,:)./gull_lake_precip)'
coor_precip_GS=cat(1,mean_all_loadingfactors_GS(1,:)./gull_lake_precip,mean_all_loadingfactors_GS(2,:)./gull_lake_precip,mean_all_loadingfactors_GS(3,:)./gull_lake_precip)'
coor_norm_diff=cat(1,((mean_all_loadingfactors(1,:)./mean(mean_all_loadingfactors(1,:)))-1)-((mean_all_loadingfactors_GS(1,:)./mean(mean_all_loadingfactors_GS(1,:)))-1),((mean_all_loadingfactors(2,:)./mean(mean_all_loadingfactors(2,:)))-1)-((mean_all_loadingfactors_GS(2,:)./mean(mean_all_loadingfactors_GS(2,:)))-1),((mean_all_loadingfactors(3,:)./mean(mean_all_loadingfactors(3,:)))-1)-((mean_all_loadingfactors_GS(3,:)./mean(mean_all_loadingfactors_GS(3,:)))-1))'



[a,P]=corrcoef(coor,'rows','complete')
mean([a(1,2),a(1,3),a(2,3)])
[c,P]=corrcoef(coor_GS,'rows','complete')
mean([c(1,2),c(1,3),c(2,3)])
[b,P]=corrcoef(coor_precip,'rows','complete')
mean([b(1,2),b(1,3),b(2,3)])
[d,P]=corrcoef(coor_precip_GS,'rows','complete')
mean([d(1,2),d(1,3),d(2,3)])
[e,P]=corrcoef(coor_norm_diff,'rows','complete')
mean([e(1,2),e(1,3),e(2,3)])



%%%% COLOR GUIDE
%Water (40, 109, 238) [0.1569 0.4275 0.9333]
%Light (253, 240, 90) [0.9922 0.9412 0.3529]
%Carbon (238, 70, 81) [0.9333 0.2745 0.3176]
%%%%%%%%%%%%%%%%%%%%%%%% ALL HAIL FIGURE THREE

subplot(3,1,1)
hold on
plot(xtime,mean_all_loadingfactors(1,:),'-', 'color', [0.1569 0.4275 0.9333])

yyaxis right
plot(xtime,mean_all_loadingfactors_GS(1,:),':', 'color', [0.1569 0.4275 0.9333],'HandleVisibility','off')

yyaxis left
plot(xtime,mean_all_loadingfactors(2,:),'-', 'color', [0.9922 0.9412 0.3529])

yyaxis right
plot(xtime,mean_all_loadingfactors_GS(2,:),':','color', [0.9922 0.9412 0.3529],'HandleVisibility','off')

yyaxis left
plot(xtime,mean_all_loadingfactors(3,:),'-', 'color', [0.9333 0.2745 0.3176])

yyaxis right
plot(xtime,mean_all_loadingfactors_GS(3,:),':', 'color', [0.9333 0.2745 0.3176],'HandleVisibility','off')

title('Annual (n=7)')

ylabel('Growing Season Loading Factors')
yyaxis left
ylabel('Loading Factors')
axis([2008.5 2016.5 -.05 0.35])
yticks([0 0.1 0.2])

yyaxis right
axis([2008.5 2016.5 -3 4.5])
yticks([0 1 2 3 4])
%xticklabels({'2009','2010','2011','2012','2013','2014','2015','2016'})
box on
text(2016,1.15,'(a)')

legend('WUE','LUE','CUE')
legend('boxoff')
legend('Location','northwest')




%Monthly Plots

subplot(3,1,2)
hold on

plotx1=cat(1,mean_loadingfactor_jan(1),mean_loadingfactor_feb(1),mean_loadingfactor_mar(1),mean_loadingfactor_apr(1),mean_loadingfactor_may(1),mean_loadingfactor_jun(1),mean_loadingfactor_jul(1),mean_loadingfactor_aug(1),mean_loadingfactor_sep(1),mean_loadingfactor_oct(1),mean_loadingfactor_nov(1),mean_loadingfactor_dec(1));
plot(1:12,plotx1, 'color', [0.1569 0.4275 0.9333])


plotx2=cat(1,mean_loadingfactor_jan(2),mean_loadingfactor_feb(2),mean_loadingfactor_mar(2),mean_loadingfactor_apr(2),mean_loadingfactor_may(2),mean_loadingfactor_jun(2),mean_loadingfactor_jul(2),mean_loadingfactor_aug(2),mean_loadingfactor_sep(2),mean_loadingfactor_oct(2),mean_loadingfactor_nov(2),mean_loadingfactor_dec(2));
plot(1:12,plotx2, 'color', [0.9922 0.9412 0.3529])


plotx3=cat(1,mean_loadingfactor_jan(3),mean_loadingfactor_feb(3),mean_loadingfactor_mar(3),mean_loadingfactor_apr(3),mean_loadingfactor_may(3),mean_loadingfactor_jun(3),mean_loadingfactor_jul(3),mean_loadingfactor_aug(3),mean_loadingfactor_sep(3),mean_loadingfactor_oct(3),mean_loadingfactor_nov(3),mean_loadingfactor_dec(3));
plot(1:12,plotx3, 'color', [0.9333 0.2745 0.3176])


title('Seasonal (n=8)')
ylabel('Loading Factors')
axis([0.5 12.5 -0.3 1.3])
xticks([1 3 5 7 9 11])
xticklabels({'Jan','Mar','May','July','Sep','Nov'})
box on
text(11.75,1.15,'(b)')

subplot(3,1,3)

Y_month = abs(cat(2,plotx1,plotx2,plotx3));
Y_month = bsxfun(@rdivide, Y_month, nansum(Y_month,2));
map = [0.1569, 0.4275, 0.9333
    0.9922, 0.9412, 0.3529
    0.9333, 0.2745, 0.3176];

S=area(Y_month,'LineStyle','none')
alpha(S,.7)
axis([0.5 12.5 0 1])
xticks([1 3 5 7 9 11])
xticklabels({'Jan','Mar','May','July','Sep','Nov'})
colormap(map)
ylabel('Loading Factors Weight [%]')
text(11.75,1.1,'(c)')














%%%%%%%%%% FIGURE FOUR

%Water (40, 109, 238) [0.1569 0.4275 0.9333]
%Light (253, 240, 90) [0.9922 0.9412 0.3529]
%Carbon (238, 70, 81) [0.9333 0.2745 0.3176]

Soybean = all_loadingfactors_GS(:,1,1:6);
Soybean_W = [reshape(Soybean(1,:,:),[1,6]),nan(1,8)];
Soybean_L = [reshape(Soybean(2,:,:),[1,6]),nan(1,8)];
Soybean_C = [reshape(Soybean(3,:,:),[1,6]),nan(1,8)];

Switchgrass = [all_loadingfactors_GS(:,2:8,1),all_loadingfactors_GS(:,2:8,5)];
Switchgrass_W = reshape(Switchgrass(1,:),[1,14]);
Switchgrass_L = reshape(Switchgrass(2,:),[1,14]);
Switchgrass_C = reshape(Switchgrass(3,:),[1,14]);

Prairie = [all_loadingfactors_GS(:,2:8,2),all_loadingfactors_GS(:,2:8,4)];
Prairie_W = reshape(Prairie(1,:),[1,14]);
Prairie_L = reshape(Prairie(2,:),[1,14]);
Prairie_C = reshape(Prairie(3,:),[1,14]);

Corn = [all_loadingfactors_GS(:,2:8,3),all_loadingfactors_GS(:,2:8,6)];
Corn_W = reshape(Corn(1,:),[1,14]);
Corn_L = reshape(Corn(1,:),[1,14]);
Corn_C = reshape(Corn(1,:),[1,14]);

Grassland = all_loadingfactors_GS(:,2:8,7);
Grassland_W = [Grassland(1,:,:),nan(1,7)];
Grassland_L = [Grassland(2,:,:),nan(1,7)];
Grassland_C = [Grassland(3,:,:),nan(1,7)];



subplot(3,2,1)
boxplot([Soybean_W;Switchgrass_W;Prairie_W;Corn_W;Grassland_W]','symbol','x','color',[0.1569 0.4275 0.9333],'Labels',{'Soybean','Switchgrass','Prairie','Corn','Grassland'})
axis([0.5 5.5 -.5 8.5])
ylabel('Growing Season Loading Factors')
text(0.75,7.5,'(a)')

subplot(3,2,3)
boxplot([Soybean_L;Switchgrass_L;Prairie_L;Corn_L;Grassland_L]','symbol','x','color',[0.9922 0.9412 0.3529],'Labels',{'Soybean','Switchgrass','Prairie','Corn','Grassland'})
axis([0.5 5.5 -.5 8.5])
ylabel('Growing Season Loading Factors')
text(0.75,7.5,'(b)')

subplot(3,2,5)
boxplot([Soybean_C;Switchgrass_C;Prairie_C;Corn_C;Grassland_C]','symbol','x','color',[0.9333 0.2745 0.3176],'Labels',{'Soybean','Switchgrass','Prairie','Corn','Grassland'})
axis([0.5 5.5 -.5 8.5])
ylabel('Growing Season Loading Factors')
text(0.75,7.5,'(c)')



subplot(3,2,[2,4,6])

stem((2009:2016)-.1,((mean_all_loadingfactors(1,:)./mean(mean_all_loadingfactors(1,:)))-1)-((mean_all_loadingfactors_GS(1,:)./mean(mean_all_loadingfactors_GS(1,:)))-1),'color',[0.1569 0.4275 0.9333])
hold on
stem((2009:2016)+.0,((mean_all_loadingfactors(2,:)./mean(mean_all_loadingfactors(2,:)))-1)-((mean_all_loadingfactors_GS(2,:)./mean(mean_all_loadingfactors_GS(2,:)))-1),'color',[0.9922 0.9412 0.3529])
stem((2009:2016)+.1,((mean_all_loadingfactors(3,:)./mean(mean_all_loadingfactors(3,:)))-1)-((mean_all_loadingfactors_GS(3,:)./mean(mean_all_loadingfactors_GS(3,:)))-1),'color',[0.9333 0.2745 0.3176])

ylabel({'Normalized Difference Between','Annual and Growing Season Loading Factors'})
axis([2008.1 2016.9 -1.9 4.9])
view(90,90)
%XAxisLocation('top')
legend('WUE','LUE','CUE')
legend('boxoff')
legend('Location','southeast')
text(2008.5,4.2,'(d)')

















