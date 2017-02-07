% Plot the results of the correlation analyses
%
% By Daniel Jeck 2016

load('corr_rand') % correlation results of randomly shuffled pictures (e.g. pic 1 fixation with pic 2 interest)
load('corr_true'); % correlation results of matched pictures (e.g. pic 1 fixation and interest)
load('R_samp_err_fixtap'); % correlation of fixations resampled with the number of taps (sample error only cause of <1 correlation)
load('R_samp_err_intfix'); % correlation of interest resampled with the number of fixation
load('R_samp_err_inttap'); % correlation of interest resampled with the number of taps
load('R_samp_err_tapsal');
load('R_intNtaps_fix2');
load('R_samp_err_fixsal');
load('R_samp_err_intsal');

Ravg = mean(R_true,3);
figure(1);

%% Fixations vs interest

subplot(2,3,1)

alpha = 0.5;
xlist = [-0.5:0.05:1];

[pRfixint_rand] = phistf(Rfixint_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
[pRfixint] = phistf(Rfixint(:),xlist,'FaceColor','r','FaceAlpha',alpha);
[pRfixtap_samp_err] = phistf(mean(R_samp_err_intfix,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRfixint_rand,'b')
hold all;
plot(xlist,pRfixint,'r')
plot(xlist,pRfixtap_samp_err,'k')
hold off

% xlabel('Correlation Value');
ylabel('Relative Frequency');

% compute p-values
[~, p_Rfixint] = ztest(mean(Rfixint),mean(Rfixint_rand(:)), ...
    sqrt(var(Rfixint(:))/length(Rfixint) + var(Rfixint_rand(:))/length(Rfixint_rand(:))),0.95,'right'); %two sampled z-test by computing sigma as sqrt(var1/n1+var2/n2)
[~, p_Rfixint_samp_err] = ztest(mean(Rfixint),mean(R_samp_err_intfix(:)), ...
    sqrt(var(Rfixint(:))/length(Rfixint) + var(mean(R_samp_err_intfix))/length(mean(R_samp_err_intfix))),0.95,'both');
p_Rfixint_samp_err2 = (sum(mean(Rfixint)>=mean(R_samp_err_intfix,2))+1)/(size(R_samp_err_intfix,1)+1);

disp(['pfixint = ' num2str(p_Rfixint) ' / ' num2str(p_Rfixint_samp_err) '/' num2str(p_Rfixint_samp_err2)]);
title('Fixation vs Interest');

y = 0.5;

% subplot(2,2,1)
hold all;

herrorbar(mean(Rfixint),y,std(Rfixint)/sqrt(length(Rfixint)),'or');
herrorbar(mean(Rfixint_rand(:)),y,std(Rfixint_rand(:))/sqrt(length(Rfixint)),'ob');
herrorbar(mean(R_samp_err_intfix(:)),y,std(R_samp_err_intfix(:))/sqrt(length(Rfixint)),'ok');
ylim([0 y+.1])
hold off



%% Fixations vs taps
figure(1);
% subplot(2,4,2);
% plot(sort(Rfixtap_rand(:)))
% hold on;
% plot(0.95*length(Rfixtap_rand(:)),Rfixtap,'r*');
% plot(0.95*length(Rfixtap_rand(:)),Ravg(1,3),'g*');
% errorbar(0.95*length(Rfixtap_rand(:)),Ravg(1,3),...
%     1.96*std(Rfixtap)/sqrt(length(Rfixtap)),'g');
% plot(0.95*length(Rfixtap_rand(:)),mean(R_samp_err_fixtap(:)),'k*');
% plot(0.95*length(Rfixtap_rand(:)),mean(R_intNtaps_fix(:)),'b*');
% hold off

subplot(2,3,4)

pRfixtap_rand = phistf(Rfixtap_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
pRfixtap = phistf(Rfixtap(:),xlist,'FaceColor','r','FaceAlpha',alpha);
pRfixtap_samp_err = phistf(mean(R_samp_err_fixtap,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRfixtap_rand,'b')
hold all;
plot(xlist,pRfixtap,'r')
plot(xlist,pRfixtap_samp_err,'k')
hold off

[~, p_Rfixtap] = ztest(mean(Rfixtap),mean(Rfixtap_rand(:)),...
    sqrt(var(Rfixtap(:))/length(Rfixtap) + var(Rfixtap_rand(:))/length(Rfixtap_rand(:))),0.95,'right');
[~, p_Rfixtap_samp_err] = ztest(mean(Rfixtap),mean(R_samp_err_fixtap(:)),...
    sqrt(var(Rfixtap(:))/length(Rfixtap) + var(mean(R_samp_err_fixtap))/length(mean(R_samp_err_fixtap))),0.95,'both');
p_Rfixtap_samp_err2 = (sum(mean(Rfixtap)>=mean(R_samp_err_fixtap,2))+1)/(size(R_samp_err_fixtap,1)+1);

title('Fixation vs Taps');
xlabel('Correlation Value');
ylabel('Relative Frequency')

disp(['pfixtap = ' num2str(p_Rfixtap) ' / ' num2str(p_Rfixtap_samp_err) '/' num2str(p_Rfixtap_samp_err2)]);

y = 0.5;
% subplot(2,3,4)
hold all;

herrorbar(mean(Rfixtap),y,std(Rfixtap)/sqrt(length(Rfixtap)),'or');
herrorbar(mean(Rfixtap_rand(:)),y,std(Rfixtap_rand(:))/sqrt(length(Rfixtap)),'ob');
herrorbar(mean(R_samp_err_fixtap(:)),y,std(R_samp_err_fixtap(:))/sqrt(length(Rfixtap)),'ok');
ylim([0 y+.1])

hold off



%% Interest vs taps
figure(1);
% subplot(2,4,3);
% plot(sort(Rinttap_rand(:)))
% hold on;
% plot(0.95*length(Rinttap_rand(:)),Rinttap,'r*');
% plot(0.95*length(Rinttap_rand(:)),Ravg(2,3),'g*');
% errorbar(0.95*length(Rinttap_rand(:)),Ravg(2,3),...
%     1.96*std(Rinttap)/sqrt(length(Rinttap)),'g');
% plot(0.95*length(Rinttap_rand(:)),mean(R_samp_err_inttap(:)),'k*');
% 
% hold off

subplot(2,3,5)

pRinttap_rand = phistf(Rinttap_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
pRinttap = phistf(Rinttap(:),xlist,'FaceColor','r','FaceAlpha',alpha);
pRinttap_samp_err = phistf(mean(R_samp_err_inttap,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRinttap_rand,'b')
hold all;
plot(xlist,pRinttap,'r')
plot(xlist,pRinttap_samp_err,'k')
hold off


[~, p_Rinttap] = ztest(mean(Rinttap),mean(Rinttap_rand(:)),...
    sqrt(var(Rinttap(:))/length(Rinttap) + var(Rinttap_rand(:))/length(Rinttap_rand(:))),0.95,'right');
[~, p_Rinttap_samp_err] = ztest(mean(Rinttap),mean(R_samp_err_inttap(:)),...
    sqrt(var(Rinttap(:))/length(Rinttap) + var(mean(R_samp_err_inttap))/length(mean(R_samp_err_inttap))),0.95,'both');
p_Rinttap_samp_err2 = (sum(mean(Rinttap)>=mean(R_samp_err_inttap,2))+1)/(size(R_samp_err_inttap,1)+1);



title('Interest vs Taps');
xlabel('Correlation Value');
% ylabel('Relative Frequency');
disp(['pinttap = ' num2str(p_Rinttap) ' / ' num2str(p_Rinttap_samp_err) ' / ' num2str(p_Rinttap_samp_err2)]);

y = 0.5;
% subplot(2,3,5)
hold all;

herrorbar(mean(Rinttap),y,std(Rinttap)/sqrt(length(Rinttap)),'or');
herrorbar(mean(Rinttap_rand(:)),y,std(Rinttap_rand(:))/sqrt(length(Rinttap)),'ob');
herrorbar(mean(R_samp_err_inttap(:)),y,std(R_samp_err_inttap(:))/sqrt(length(Rinttap)),'ok');
ylim([0 y+.1])

hold off


%% Taps vs Sal

figure(1);
% subplot(2,4,4);
% plot(sort(Rtapsal_rand(:)))
% 
% hold on;
% plot(0.95*length(Rtapsal_rand(:)),Rtapsal,'r*');
% plot(0.95*length(Rtapsal_rand(:)),Ravg(3,4),'g*');
% errorbar(0.95*length(Rtapsal_rand(:)),Ravg(3,4),...
%     1.96*std(Rtapsal)/sqrt(length(Rtapsal)),'g');
% plot(0.95*length(Rtapsal_rand(:)),mean(R_samp_err_tapsal(:)),'k*');
% 
% hold off

subplot(2,3,6)

pRtapsal_rand = phistf(Rtapsal_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
pRtapsal = phistf(Rtapsal(:),xlist,'FaceColor','r','FaceAlpha',alpha);
pRtapsal_samp_err = phistf(mean(R_samp_err_tapsal,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRtapsal_rand,'b')
hold all;
plot(xlist,pRtapsal,'r')
plot(xlist,pRtapsal_samp_err,'k')
hold off

% legend('Null hypothesis','Measured correllation','Upper bound');

[~, p_Rtapsal] = ztest(mean(Rtapsal),mean(Rtapsal_rand(:)),...
    sqrt(var(Rtapsal(:))/length(Rtapsal) + var(Rtapsal_rand(:))/length(Rtapsal_rand(:))),0.95,'right');
[~, p_Rtapsal_samp_err] = ztest(mean(Rtapsal),mean(R_samp_err_tapsal(:)),...
    sqrt(var(Rtapsal(:))/length(Rtapsal) + var(mean(R_samp_err_tapsal))/length(mean(R_samp_err_tapsal))),0.95,'both');
p_Rtapsal_samp_err2 = (sum(mean(Rtapsal)>=mean(R_samp_err_tapsal,2))+1)/(size(R_samp_err_tapsal,1)+1);

title('Comp Saliency vs Taps');
xlabel('Correlation Value');

disp(['ptapsal = ' num2str(p_Rtapsal) ' / ' num2str(p_Rtapsal_samp_err) ' / ' num2str(p_Rtapsal_samp_err2)]);

legend('Null hypothesis','Measured correllation','Sample Error hypothesis');

y=0.5;
% subplot(2,3,6)
hold all;

herrorbar(mean(Rtapsal),y,std(Rtapsal)/sqrt(length(Rtapsal)),'or');
herrorbar(mean(Rtapsal_rand(:)),y,std(Rtapsal_rand(:))/sqrt(length(Rtapsal)),'ob');
herrorbar(mean(R_samp_err_tapsal(:)),y,std(R_samp_err_tapsal(:))/sqrt(length(Rtapsal)),'ok');
ylim([0 y+.1])

hold off


%% Fixation vs Salience
%

figure(1);
subplot(2,3,2);
[pRfixsal_rand] = phistf(Rfixsal_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
[pRfixsal] = phistf(Rfixsal(:),xlist,'FaceColor','r','FaceAlpha',alpha);
[pRfixsal_samp_err] = phistf(mean(R_samp_err_fixsal,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRfixsal_rand,'b')
hold all;
plot(xlist,pRfixsal,'r')
plot(xlist,pRfixsal_samp_err,'k')
hold off

% xlabel('Correlation Value');
% ylabel('Relative Frequency');

[~, p_Rfixsal] = ztest(mean(Rfixsal),mean(Rfixsal_rand(:)),...
    sqrt(var(Rfixsal(:))/length(Rfixsal) + var(Rfixsal_rand(:))/length(Rfixsal_rand(:))),0.95,'right');
[~, p_Rfixsal_samp_err] = ztest(mean(Rfixsal),mean(R_samp_err_fixsal(:)),...
    sqrt(var(Rfixsal(:))/length(Rfixsal) + var(mean(R_samp_err_fixsal))/length(mean(R_samp_err_fixsal))),0.95,'both');
p_Rfixsal_samp_err2 = (sum(mean(Rfixsal)>=mean(R_samp_err_fixsal,2))+1)/(size(R_samp_err_fixsal,1)+1);


y=0.5;

disp(['pfixsal = ' num2str(p_Rfixsal) ' / ' num2str(p_Rfixsal_samp_err) ' / ' num2str(p_Rfixsal_samp_err2)]);
title('Fixation vs Comp Saliency');
hold on;
herrorbar(mean(Rfixsal),y,std(Rfixsal)/sqrt(length(Rfixsal)),'or');
herrorbar(mean(Rfixsal_rand(:)),y,std(Rfixsal_rand(:))/sqrt(length(Rfixsal)),'ob');
herrorbar(mean(R_samp_err_fixsal(:)),y,std(R_samp_err_fixsal(:))/sqrt(length(Rfixsal)),'ok');
ylim([0 y+.1])

hold off;

%% Int vs Salience

figure(1);
subplot(2,3,3);
[pRintsal_rand] = phistf(Rintsal_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
xlim([-0.2 1]);
hold on;
[pRintsal] = phistf(Rintsal(:),xlist,'FaceColor','r','FaceAlpha',alpha);
[pRintsal_samp_err] = phistf(mean(R_samp_err_intsal,1),xlist,'FaceColor','k','FaceAlpha',alpha);
hold off;

plot(xlist,pRintsal_rand,'b')
hold all;
plot(xlist,pRintsal,'r')
plot(xlist,pRintsal_samp_err,'k')
hold off

% xlabel('Correlation Value');
% ylabel('Relative Frequency');

[~, p_Rintsal] = ztest(mean(Rintsal),mean(Rintsal_rand(:)),...
    sqrt(var(Rintsal(:))/length(Rintsal) + var(Rintsal_rand(:))/length(Rintsal_rand(:))),0.95,'right');
[~, p_Rintsal_samp_err] = ztest(mean(Rintsal),mean(R_samp_err_intsal(:)),...
    sqrt(var(Rintsal(:))/length(Rintsal) + var(mean(R_samp_err_intsal))/length(mean(R_samp_err_intsal))),0.95,'both');
p_Rintsal_samp_err2 = (sum(mean(Rintsal)>=mean(R_samp_err_intsal,2))+1)/(size(R_samp_err_intsal,1)+1);

y=0.5;

disp(['pintsal = ' num2str(p_Rintsal) ' / ' num2str(p_Rintsal_samp_err) '/' num2str(p_Rintsal_samp_err2)]);
title('Interest vs Comp Saliency');
hold on;
herrorbar(mean(Rintsal),y,std(Rintsal)/sqrt(length(Rintsal)),'or');
herrorbar(mean(Rintsal_rand(:)),y,std(Rintsal_rand(:))/sqrt(length(Rintsal)),'ob');
herrorbar(mean(R_samp_err_intsal(:)),y,std(R_samp_err_intsal(:))/sqrt(length(Rintsal)),'ok');
ylim([0 y+.1])

hold off;
subplot(2,3,5)
%% Fix vs. Russell Salience
% 
% pRfixruss_rand = phistf(Rfixruss_rand(:),xlist,'FaceColor','b','FaceAlpha',alpha);
% % xlim([-0.2 1]);
% hold on;
% pRfixruss = phistf(Rfixruss(:),xlist,'FaceColor','r','FaceAlpha',alpha);
% pRfixruss_samp_err = phistf(mean(R_samp_err_fixruss,2),xlist,'FaceColor','k','FaceAlpha',alpha);
% hold off;
% 
% plot(xlist,pRfixruss_rand,'b')
% hold all;
% plot(xlist,pRfixruss,'r')
% plot(xlist,pRfixruss_samp_err,'k')
% hold off
% 
% legend('Null hypothesis','Measured correllations','Sample Error hypothesis');
% 
% [~, p_Rfixruss] = ztest(mean(Rfixruss),mean(Rfixruss_rand(:)),std(Rfixruss_rand(:))/sqrt(length(Rfixruss)),0.95,'right');
% [~, p_Rfixruss_samp_err] = ztest(mean(Rfixruss),mean(R_samp_err_fixruss(:)),std(mean(R_samp_err_fixruss,2)),0.95,'both');
% 
% 
% title('Russ Sal vs Fix');
% xlabel('Correlation Value');
% 
% disp(['pfixruss = ' num2str(p_Rfixruss) ' / ' num2str(p_Rfixruss_samp_err)]);


% y=0.4;
% % subplot(2,3,5)
% hold all;
% 
% herrorbar(mean(Rfixruss),y,std(Rfixruss)/sqrt(length(Rfixruss)),'or');
% herrorbar(mean(Rfixruss_rand(:)),y,std(Rfixruss_rand(:))/sqrt(length(Rfixruss)),'ob');
% herrorbar(mean(R_samp_err_fixruss(:)),y,std(R_samp_err_fixruss(:))/sqrt(length(Rfixruss)),'ok');
% 
% hold off

%%
boldify

%%
figure(2)
pic = 64-30;

% for pic = 1:48
pthresh = 0.025; 

y=0.5;
hist(R_samp_err_fixtap(:,pic));

Pupper = phistf(R_samp_err_fixtap(:,pic),xlist);
Plower = phistf(Rfixtap_rand(:),xlist);

upper_err = std(R_samp_err_fixtap(:,pic));
lower_err = std(Rfixtap_rand(:));

% upper_sort = sort(R_upper_fixtap(:,pic));
% upper_L = upper_sort(ceil(pthresh*length(upper_sort)))-mean(upper_sort);
% upper_sort = sort(R_upper_fixtap(:,pic),'descend');
% upper_U = upper_sort(ceil(pthresh*length(upper_sort)))-mean(upper_sort);

% lower_sort = sort(Rfixtap_rand(:));
% lower_L = lower_sort(ceil(pthresh*length(lower_sort)))-mean(lower_sort);
% lower_sort = sort(Rfixtap_rand(:),'descend');
% lower_U = lower_sort(ceil(pthresh*length(lower_sort)))-mean(lower_sort);



plot(xlist,Pupper,'k');
hold on;
plot(xlist,Plower,'b');
herrorbar(mean(Rfixtap_rand(:)),y,lower_err,'ob')
herrorbar(mean(R_samp_err_fixtap(:,pic)),y,upper_err,'ok')
stem(Rfixtap(pic),y,'or');
hold off;
ylim([0 0.6])
xlim([-0.2 1]);

xlabel('Correlation Value');
ylabel('Relative Frequency');

pnull = sum(Rfixtap_rand(:)>Rfixtap(pic))/length(Rfixtap_rand(:))

p_samp_err = sum(R_samp_err_fixtap(:,pic)<Rfixtap(pic))/length(R_samp_err_fixtap(:,pic))
% title(num2str(pic));

boldify
% pause;
% end