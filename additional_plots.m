
%%%
% plot SD excursions and bin descretization

clear all
close all

fp = 'D:\My Files\Work\BGU\datasets\Panas\';
[eegFile, fp] = uigetfile([fp '*.set'], 'Select set file');
EEG = pop_loadset([fp eegFile]);


chan = 1;
eegData = zscore(reshape(EEG.data,size(EEG.data,1),[]),0,2);
N = size(eegData,2);
t = 0:1/EEG.srate:(N-1)/EEG.srate;
bin = 10;
maxval = max(eegData(chan,:));
minval = min(eegData(chan,:));
binData = nan(1,N);

figure; 
ax(1) = subplot(2,1,1);plot(t,eegData(chan,:),'b');xlabel('sec');ylabel('Z-score'); title('Fz - all trials'); hold on; 
line([0,(N-1)/EEG.srate],[3,3], 'LineStyle',':', 'Color','k', 'LineWidth',1);
line([0,(N-1)/EEG.srate],[-3,-3], 'LineStyle',':', 'Color','k', 'LineWidth',1);
for i = bin:bin:N
    line([i/EEG.srate,i/EEG.srate],[minval,maxval], 'LineStyle','--', 'Color','m');
end
box off;set(gca,'FontSize',16);
hold off;

ax(2) = subplot(2,1,2);hold on;
for i = bin:bin:N
    line([i/EEG.srate,i/EEG.srate],[-0.5,1.5], 'LineStyle','--', 'Color','m');

    if any(abs(eegData(chan,i-bin+1:i)) >= 3)
        binData(i-bin+1:i) = 1;
    end
end
plot(t,binData,'k', 'LineWidth',1);xlabel('sec');title('Time bins of length 39ms discretization');
hold off;
box off;set(gca,'FontSize',16);

linkaxes(ax,'x');


%%%%%%%%%%%%%%%%%%%%%%%
%plot families

figure;
rasterMat = zeros(60,5);
rasterMat([9,12,13,14,16,20,28,32,46,47,48,49,53,55,59],2) = 1;
rasterMat([10,12,17,25,30,33,46,47,49,56],3) = 1;
subplot(1,5,1);spy(rasterMat,'k.',10);xlabel('');ylabel('channels');set(gca,'FontSize',16);

rasterMat = zeros(60,5);
rasterMat([10,11,15,26,30,32,34,46,47,48,50,55],2) = 1;
rasterMat([10,12,17,25,30,33,46,47,49,56],3) = 1;
subplot(1,5,2);spy(rasterMat,'k.',10);xlabel('');set(gca,'FontSize',16);

rasterMat = zeros(60,5);
rasterMat([2,5,18,19,21,25,37,39,40,58,59],1) = 1;
rasterMat([10,11,14,15,25,30,31,32,33,34,47,48,49,50,56],2) = 1;
rasterMat([10,11,15,26,30,32,34,46,47,48,50,55],3) = 1;
rasterMat([10,12,17,25,30,33,46,47,49,56],4) = 1;
subplot(1,5,3);spy(rasterMat,'k.',10);xlabel('');set(gca,'FontSize',16);

rasterMat = zeros(60,5);
rasterMat([11,14,30,31,32,33,34,47],2) = 1;
rasterMat([10,11,15,26,30,32,34,46,47,48,50,55],3) = 1;
rasterMat([10,12,17,25,30,33,46,47,49,56],4) = 1;
subplot(1,5,4);spy(rasterMat,'k.',10);xlabel('');set(gca,'FontSize',16);

rasterMat = zeros(60,5);
rasterMat([5,11,14,26,31,32,34,46,47,48,49,55],3) = 1;
subplot(1,5,5);spy(rasterMat,'k.',10);xlabel('');set(gca,'FontSize',16);
