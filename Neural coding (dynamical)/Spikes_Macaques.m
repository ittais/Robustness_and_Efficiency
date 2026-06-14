function [muAll, stdAll, muMat, stdMat]=Spikes(validExpIDs,Nchannels, EN)
%This function converts raw spike timestamps from one recording session into dynamic neuron-neuron correlation networks, 
% computes spectral and graph-theoretic measures (strength, entropy, Fiedler value, spectral gap, energy) 
% for dACC, amygdala, and whole-brain networks across time, and returns regional and whole-network summary 
% statistics.
%This function is run with varying temporal and spatial resolutions prior
%to figures plots.

%% ============================================================
% SECTION 1 — LOAD INPUT DATA
% Loads spike timestamp metadata and event timing data.
%% ============================================================

load('neuron_cell.mat');
load('event_cell.mat'); %load electrophysiological data and labeling


%% ============================================================
% SECTION 2 — SESSION METADATA AND CURRENT EXPERIMENT SELECTION
% Extracts recording dates, regions, labels, and selects the session
% specified by validExpIDs(EN).
%% ============================================================
dates=cell2mat(cell_neurons(2:end,4));
numdates=length(unique(dates)); %69 dates (~3 hours each)
[U,Dates1,Dates2]=unique(dates);

monkeys=cell2table(cell_neurons(2:end,2));
Regions=cell2table(cell_neurons(2:end,8));

location=cell2table(cell_neurons(2:end,8)); %region (dACC/Amygdala/SI)
[locations,a1,a2]=unique(location);

labelsall=(cell_neurons(2:end,8));
Labelsall=string(labelsall);

% i=4; %date
%[dataAllTMatb,Labelsb,removal] = FILTERandSMOOTH1(i,cell_neurons,Dates1,labelsall);

%t_On = 0; %time stimulus turns on [ms]
%t_Move = 500; %time stimulus begins moving [ms]
%t_Off = 2500; %time stimulus turns off [ms]
%NumAngles = size(spikes,1) - NumControls %number of angles tested, equally spaced;
% %last 2 sets of recordings are controls
i=validExpIDs(EN); %1; %validExpIDs(EN); %8; %4; %8- Amy-dACC BDI (i=4 Amy-SI BCI) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 1; %spacing between sampled time points [ms]
NumTimePoints = 11*10^6; %30hr (=10.8M ms) %size(spikes,2) %number of time points; time was sampled every 1 ms
startTrial=Dates1(i); %first trial in a specific date
endTrial=Dates1(i+1)-1; %last trial in the same specific date
labels=labelsall(startTrial:endTrial);
Labels=string(labels);
%NumTrials = size(spikes,3) %number of trials performed at each angle



%% ============================================================
% SECTION 3 — RAW SPIKE TIMES TO BINARY SPIKE MATRIX
% Converts each channel's spike timestamps into a binary vector over time.
%% ============================================================
t_All = 1:NumTimePoints;
t_vect = 1:(endTrial-startTrial+1); %1:NumTimePoints; %t_On:dt:(NumTimePoints-1)*dt; %time vector for each trial

figure;%(1) RAW DATA RASTER PLOT
dataAllTMat=zeros(endTrial-startTrial+1,NumTimePoints);
%location=table2array(location);
triali=1;
for trial=startTrial:endTrial %length(cell_neurons)
 dataAll=cell_neurons{trial+1,21}; 
 dataAllT=zeros(1,NumTimePoints); dataAllT(dataAll)=1;
 % % plot(t_All,(triali-1)*dataAllT,'+'); %trial*spikes(2,:,trial),'+') %raster plot
 dataAllTMat(triali,:)=dataAllT; %spikes 
 triali=triali+1;
 %data1=cell_neurons{trial,21}(startPoint:endPoint);
 %data=zeros(size(t_vect)); % data(data1)=1;
 %plot(t_vect,data,'+'); %trial*spikes(2,:,trial),'+') %raster plot
 hold on
end
hold off; xlabel('time (ms)'); xlim([0 9*10^6]);
yticks(0:(endTrial-startTrial)); yticklabels(cell_neurons(1+(startTrial:endTrial),8)');
ylabel('Channel'); title('Raster plot- spike sorted'); set(gcf,'Color','white');


%% ============================================================
% SECTION 4 — ORDER CHANNELS BY REGION
% Splits channels into dACC and amygdala and concatenates them in region order.
%% ============================================================
Ldacc=find(Labels=='dacc');
Lamy=find(Labels=='amy');
% Lsi=find(Labels=='si');

%ORDER CHANNELS:
dataAllTMatdacc=dataAllTMat(Ldacc,:); 
dataAllTMatamy=dataAllTMat(Lamy,:);
% dataAllTMatsi=dataAllTMat(Lsi,:);

dataAllTMatb=cat(1,dataAllTMatdacc,dataAllTMatamy);%,dataAllTMatsi);
Labelsb=cat(2,repmat("dacc",1,(length(Ldacc))),repmat("amy",1,(length(Lamy))));%,repmat("si",1,(length(Lsi))));

%SBP (spiking band pass) on original spike data (ordered by region):

% 
% figure;%(1) ORDERED DATA RASTER PLOT
% for i=1:size(dataAllTMatb,1)  %length(cell_neurons)
%  % % plot(t_All,(i-1)*dataAllTMatb(i,:),'+'); %trial*spikes(2,:,trial),'+') %raster plot
%  hold on
% end
% hold off; xlabel('time (ms)'); xlim([0 9*10^6]);
% yticks(0:(endTrial-startTrial)); yticklabels(Labelsb);
% ylabel('Channel'); title('Raster plot- ordered by region'); set(gcf,'Color','white');

%FIRING RATES (per channel, across all recording-> for finding zero/noise):

%% ============================================================
% SECTION 5 — FIRING-RATE QUALITY CONTROL
% Estimates firing rates across the recording and marks low-firing channels
% for removal.
%% ============================================================
T=100;
jump=size(dataAllTMatb,2)/T;
firingrateT=zeros(size(dataAllTMatb,1),T); %firing rate every 1/10 of the recording
for i=1:size(dataAllTMatb,1)
    for j=1:T
        firingrateT(i,j)=sum(dataAllTMatb(i,(1+(j-1)*jump):(j*jump)),2)./(jump*10^-3);
    end
end

[rows,columns,v]=find(firingrateT(:,1:85)<0.5);
length(unique(rows))

% figure;
% for i=1:size(dataAllTMatb,1)
%     subplot(size(dataAllTMatb,1),1,i);
%     bar(1:T,firingrateT(i,:));
%     ylabel(num2str(i),'FontSize',6); hold on;
% end
% set(gcf,'Color','white'); title('Firing rates across time scales');


firingrate=sum(dataAllTMatb,2)./(length(dataAllTMatb)*10^-3); %across entire recording
firingrateT=sum(dataAllTMatb(:,3*10^6:7*10^6),2)./(length(dataAllTMatb(:,3*10^6:7*10^6))*10^-3); %within specific timeframe

% figure;
% bar(1:length(firingrate),firingrate); title('Average firing rates (per channel)');
% set(gca, 'XTick', 1:length(firingrate),'Fontsize',10); xtickangle(0);
% xlabel('Channel'); ylabel('Frequency [Hz]'); set(gcf,'Color','white');
% %figure; histogram(firingrate,length(firingrate));

thresh=0.5; %prctile(firingrate,30); %percentile (30%)
% hold on; yline(thresh);
 %mark channels for removal
thresh2=prctile(firingrate,95); %percentile (30%)
% hold on; yline(0.5); %yline(thresh2);

removal=zeros(1,length(firingrateT));
%removal(firingrate>thresh2)=1;
removal(firingrateT<thresh)=1; %1 for removal, o for no removal
dataAllTMatR=dataAllTMat; dataAllTMatR(removal==1,:)=[];
LabelsbR=Labelsb; LabelsbR(removal==1)=[]; %REMOVE DATA AND LABELS BELOW/ABOVE THRESHOLD

%recalculate averages etc.:

%% ============================================================
% SECTION 6 — REGION INDICES AFTER CHANNEL REMOVAL
% Recomputes dACC and amygdala channel indices after low-firing channels
% have been removed.
%% ============================================================
LdaccR=find(LabelsbR=='dacc');
LamyR=find(LabelsbR=='amy');
% LsiR=find(LabelsbR=='si');
dataAllTMatdaccR=dataAllTMatR(LdaccR,:);
dataAllTMatamyR=dataAllTMatR(LamyR,:);
% dataAllTMatsiR=dataAllTMatR(LsiR,:);


%i=4;

%% ============================================================
% SECTION 7 — SPIKE TRAINS TO FIRING-RATE TIME SERIES
% Applies a moving spike-count window and converts counts to firing rates.
%% ============================================================
startTime=1; %int32(event_cell{3,7}(1)); %start of first successful BCI
endTime=10*10^6; %int32(event_cell{3,8}(end)); %end of first successful BCI
window=500; %size of sliding window
delta=floor(window*0.1); %jump size between windows

% eventIntensity=zeros(size(dataAllTMat,1),floor((endTime-startTime)/window)); %#spikes per window in event
% for i=1:size(eventIntensity,1)
%     for j=1:size(eventIntensity,2)
%         eventIntensity(i,j)=sum(dataAllTMat(i,(startTime+1+window*(j-1)):(startTime+window*j)));
%     end
% end
eventIntensity=zeros(size(dataAllTMat));%,1),size(dataAllTMat,2)/delta); %dataAllTMat,1),floor((endTime-startTime)/(delta))); %#spikes per window in event
for i=1:size(eventIntensity,1)
    %for j=1:size(eventIntensity,2)
        %eventIntensity(i,j)=sum(dataAllTMat(i,(startTime+delta*(j-1)):(startTime+delta*(j-1)+window)));
        eventIntensity(i,:)=movsum(dataAllTMat(i,:),window);
    %end
end
%eventIntensityA=eventIntensity;
eventIntensity=eventIntensity(:,1:delta:end);
Ldacc=find(Labels=='dacc');
Lamy=find(Labels=='amy');
% Lsi=find(Labels=='si');

%PCA:
eventIntensityremoved=eventIntensity;
%eventIntensityremoved(sizesremove(1:size(eventIntensity,1)),:)=[]; %remove N<3
%Labelsremoved=Labels; %SESSIONS ARE REMOVED, NOT SPECIFIC CHANNELS (run through i)
Labelsremoved=eventIntensity; %(sizesremove(1:length(Labels),:))=[];
eventIntensitydacc=eventIntensityremoved(LdaccR,:); %R for removal of zero/noise channels
eventIntensityamy=eventIntensityremoved(LamyR,:);
% eventIntensitysi=eventIntensityremoved(LsiR,:);

eventIntensityb=cat(1,eventIntensitydacc,eventIntensityamy);%,eventIntensitysi);
eventIntensityb=eventIntensityb./(window*10^-3); %firing rate [1/sec]
eventIntensitydacc=eventIntensitydacc./(window*10^-3);
eventIntensityamy=eventIntensityamy./(window*10^-3);
% eventIntensitysi=eventIntensitysi./(window*10^-3);
Labelsb=cat(2,repmat("dacc",1,(length(LdaccR))),repmat("amy",1,(length(LamyR))));%,repmat("si",1,(length(LsiR))));
% 
%%eventIntensityb,Labelsb- data organized, filtered, smoothed


%% ============================================================
% SECTION 8 — FIRING-RATE QC PLOTS
% Plots firing-rate traces for all channels and separately by region.
%% ============================================================
%(1) %ORDERED DATA RASTER PLOT
%dataAllTMat=zeros(endTrial-startTrial+1,NumTimePoints);
%location=table2array(location);
triali=1;
figure;
for triali=1:size(eventIntensityb,1) %length(cell_neurons)
    subplot(size(eventIntensityb,1),1,triali);
    %subplot(max([length(Ldacc),length(Lamy),Length(Lsi)]),3,triali);
    plot((1:length(eventIntensityb)).*delta,eventIntensityb(triali,:)); %trial*spikes(2,:,trial),'+') %raster plot
    xlim([0 9*10^6]);
    %ylim([0 max(max(eventIntensityb))]);
    ylabel(num2str(triali),'FontSize',6); %title(Labelsb(triali)); 
    hold on;
end
set(gcf,'Color','white');
%xlabel('t(ms)'); ylabel('Firing rate (Hz)')
hold on; xline(event_cell{224,7}); %{98,7}); %start time BCIs
hold on; xline(event_cell{224,8}); %{98,8}); %end time BCIs

C=max([length(Ldacc),length(Lamy)]);%,length(Lsi)]);
figure; 
for triali=1:length(LdaccR)
        subplot(C,1,triali);
        plot((1:length(eventIntensityb)).*delta,eventIntensityb(triali,:)); %trial*spikes(2,:,trial),'+') %raster plot
        xlim([0 9*10^6]); %ylim([0 max(max(eventIntensityb))]);
        %title(Labelsb(triali)); 
        hold on;
end
sgtitle('dACC'); set(gcf,'Color','white');

figure;
for triali=1:length(LamyR)
        subplot(C,1,triali);
        plot((1:length(eventIntensityb)).*delta,eventIntensityb(triali,:)); %trial*spikes(2,:,trial),'+') %raster plot
        xlim([0 9*10^6]); %ylim([0 max(max(eventIntensityb))]);
        %title(Labelsb(triali)); 
        hold on;
end
sgtitle('Amy'); set(gcf,'Color','white');

% figure;
% for triali=1:length(LsiR)
%         subplot(C,1,triali);
%         plot((1:length(eventIntensityb)).*delta,eventIntensityb(triali,:)); %trial*spikes(2,:,trial),'+') %raster plot
%         xlim([0 9*10^6]); % ylim([0 max(max(eventIntensityb))]);
%         %title(Labelsb(triali)); 
%         hold on;
% end
% sgtitle('si'); set(gcf,'Color','white');


%% ============================================================
% SECTION 9 — AVERAGE FIRING RATES PER REGION
% Computes and plots mean dACC and amygdala firing-rate traces.
%% ============================================================
%average firing rates per region:
avgeventIntensitydacc=mean(eventIntensitydacc,1);
avgeventIntensityamy=mean(eventIntensityamy,1);
% avgeventIntensitysi=mean(eventIntensitysi,1);
avgall=cat(1,avgeventIntensitydacc,avgeventIntensityamy);%,avgeventIntensitysi);

figure;
subplot(2,1,1);
plot((1:length(eventIntensityb)).*delta,avgeventIntensitydacc); 
xlim([0 9*10^6]); title('dACC'); hold on;
subplot(2,1,2);
plot((1:length(eventIntensityb)).*delta,avgeventIntensityamy); 
xlim([0 9*10^6]); title('Amy'); hold on;
% subplot(3,1,3);
% plot((1:length(eventIntensityb)).*delta,avgeventIntensitysi); 
% xlim([0 9*10^6]); title('si'); 
sgtitle('Average firing rates per region'); set(gcf,'Color','white');


%% ============================================================
% SECTION 10 — REGIONAL AVERAGE CORRELATION MATRIX
% Correlates the average regional firing-rate traces.
%% ============================================================
Ravg=corrcoef(avgall'); Ravg(Ravg<0)=0;
figure; imagesc(Ravg); axis square; %colormap(hot);
set(gca, 'XTick', 1:3, 'XTickLabel', {'dacc','amy','si'},'Fontsize',10); xtickangle(90);
set(gca, 'YTick', 1:3, 'YTickLabel', {'dacc','amy','si'},'Fontsize',10); grid on; colorbar
set(gcf,'Color','white');
figure; histogram(Ravg,9); title('R values'); set(gcf,'Color','white');



%% ============================================================
% SECTION 11 — STATIC WHOLE-NETWORK FUNCTIONAL CONNECTIVITY
% Builds a static correlation network from a selected time range.
%% ============================================================
%GRAPH FOR ENTIRE NETWORK:
%calculate correlation/weighted connectivity matrix:

startT=1*10^6/delta; %5*10^6/delta; %3*10^6/delta; %delta units
endT=5*10^6/delta; %11 %9*10^6/delta; %7*10^6/delta;

Re=corrcoef(eventIntensityb(:,startT:endT)'); %Re=corrcoef(eventIntensityb'); 
Re(isnan(Re))=0;
thresh=prctile(Re,90); Re(Re<thresh)=0;
% Re(Re<0)=0; %Re=abs(Re); %
figure; imagesc(Re); axis square; %colormap(hot);
set(gca, 'XTick', 1:length(Re), 'XTickLabel', Labelsb,'Fontsize',10); xtickangle(90);
set(gca, 'YTick', 1:length(Re), 'YTickLabel', Labelsb ,'Fontsize',10); grid on; colorbar;
set(gcf,'Color','white');
figure; histogram(Re,100); title('R values'); set(gcf,'Color','white');


%% ============================================================
% SECTION 12 — STATIC NETWORK GRAPH MEASURES
% Computes weighted/binary efficiency and clustering using BCT functions.
%% ============================================================
%BRAIN CONNECTIVITY TOOLBOX ANALYSIS (ALTERNATIVE TO NETWORK SCALING):
%R=Rallb{4}; R(R<0)=0; %weighted
R=Re; R(R<0)=0; %weighted
Rbin=R; Rbin(Rbin~=0)=1; %binarized
Effw=efficiency_wei(R); Effw2=round(Effw,2); %global efficiency
Effb=efficiency_bin(Rbin); Effb2=round(Effb,2);
CCw=mean(clustering_coef_wu(R)); CCw2=round(CCw,2); %average clustering coefficient
CCb=mean(clustering_coef_bu(Rbin)); CCb2=round(CCb,2);

LABELS=LabelsbR;
figure; subplot(1,2,1); imagesc(R); axis square
set(gca, 'XTick', 1:length(R), 'XTickLabel', LABELS,'Fontsize',10); xtickangle(90);
set(gca, 'YTick', 1:length(R), 'YTickLabel', LABELS ,'Fontsize',10); grid on; colorbar; 
title(['Wieghted (corr.)  ','Eff.:',num2str(Effw2),'  C=',num2str(CCw2)],'FontSize', 18);
subplot(1,2,2); imagesc(Rbin); axis square
set(gca, 'XTick', 1:length(R), 'XTickLabel', LABELS,'Fontsize',10); xtickangle(90);
set(gca, 'YTick', 1:length(R), 'YTickLabel', LABELS ,'Fontsize',10); grid on; colorbar; 
title(['Binary  ','Eff.:',num2str(Effb2),'  C=',num2str(CCb2)],'FontSize', 18);
set(gcf,'Color','white');


%% ============================================================
% SECTION 13 — STATIC REGIONAL SUBGRAPHS
% Extracts regional subgraphs and computes regional graph measures.
%% ============================================================
%regional subgraphs:
Rdacc=R(LABELS=='dacc',LABELS=='dacc');
% Ramysi=R(LABELS=='amy'|LABELS=='si',LABELS=='amy'|LABELS=='si');
Effwdacc=efficiency_wei(Rdacc); CCwdacc=mean(clustering_coef_wu(Rdacc));
% Effwamysi=efficiency_wei(Ramysi); CCwamysi=mean(clustering_coef_wu(Ramysi));

LABELS=LabelsbR;
figure; subplot(1,1,1); imagesc(Rdacc); axis square
set(gca, 'XTick', 1:length(Rdacc), 'XTickLabel', LABELS(LABELS=='dacc'),'Fontsize',10); xtickangle(90);
set(gca, 'YTick', 1:length(Rdacc), 'YTickLabel', LABELS(LABELS=='dacc') ,'Fontsize',10); grid on; colorbar; 
title(['dACC subgraph  ','Eff.:',num2str(Effwdacc),'  C=',num2str(CCwdacc)],'FontSize', 18);
% subplot(1,2,2); imagesc(Ramysi); axis square
% set(gca, 'XTick', 1:length(Ramysi), 'XTickLabel', LABELS(LABELS=='amy'|LABELS=='si'),'Fontsize',10); xtickangle(90);
% set(gca, 'YTick', 1:length(Ramysi), 'YTickLabel', LABELS(LABELS=='amy'|LABELS=='si') ,'Fontsize',10); grid on; colorbar; 
% title(['Amy+SI subgraph  ','Eff.:',num2str(Effwamysi),'  C=',num2str(CCwamysi)],'FontSize', 18);
% set(gcf,'Color','white');



% %SPECIFIC TIME EVENTS:
% %calculate global efficiency and decrease from original:
% starters=event_cell{98,7}; %{224,7}; %
% enders=event_cell{98,8}; %{224,8}; %




%% ============================================================
% SECTION 14 — DYNAMIC NETWORK CONSTRUCTION
% Splits the recording into N temporal windows and computes correlation
% matrices for the full network, dACC, and amygdala.
%% ============================================================
%changes in mat-2-mat correlation (t and t-1, by flattening matrices to vectors and caulcating correlations), 
% total weight (sum of weights in matrix *1/2), 
% and total number of links (sum of binary matrix *1/2)
%across N=30 time scales in each recording (~6 min):

N=2000; %11 %number of time cycles (across 3hr recording)
interval=round((endT-startT)/N); %round(length(eventIntensity)/N); %delta units
intervalms=interval*delta;
startTime=1; %
Racross={}; Racrossdacc={}; Racrossamy={}; %Racrosssi={};%Racrosst={}; %cell of correlations across time intervals (6min)
totallin=zeros(1,N); totallindacc=zeros(1,N); totallinamy=zeros(1,N); %totallinsi=zeros(1,N);

totalw=zeros(1,N);

CCall=zeros(1,N); %clustering coefficient across time:
A=cat(1,eventIntensitydacc,eventIntensityamy);

for t=1:N
    if (startTime+interval*t)>length(eventIntensity)
        % corrt=corrcoef(eventIntensity(:,(startTime+interval*(t-1)):length(eventIntensity))');
        % % corrt=corrcoef(eventIntensityb(:,(startT+interval*(t-1)):length(eventIntensity))');
        corrt=corrcoef(A(:,(startT+interval*(t-1)):length(eventIntensity))');
        corrtdacc=corrcoef(eventIntensitydacc(:,(startT+interval*(t-1)):length(eventIntensitydacc))');
        corrtamy=corrcoef(eventIntensityamy(:,(startT+interval*(t-1)):length(eventIntensityamy))');
        % corrtsi=corrcoef(eventIntensitysi(:,(startT+interval*(t-1)):length(eventIntensitysi))');
    else
        % corrt=corrcoef(eventIntensity(:,(startTime+interval*(t-1)):(startTime+interval*t))');
        % corrt=corrcoef(eventIntensityb(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrt=corrcoef(A(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrtdacc=corrcoef(eventIntensitydacc(:,(startT+interval*(t-1)):(startT+interval*t))');
        corrtamy=corrcoef(eventIntensityamy(:,(startT+interval*(t-1)):(startT+interval*t))');
        % corrtsi=corrcoef(eventIntensitysi(:,(startT+interval*(t-1)):(startT+interval*t))');
    end
    corrt(isnan(corrt))=0; %replace nans with zeros
    % % % corrt(corrt<0)=0; %replace negative values with zeros %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    corrtdacc(isnan(corrtdacc))=0; %replace nans with zeros
    corrtamy(isnan(corrtamy))=0; %replace nans with zeros
    % corrtsi(isnan(corrtsi))=0; %replace nans with zeros

    CCall(t)=mean(clustering_coef_wu(corrt)); %mean clustering coefficient

    totalw(t)=0.5*sum(sum(corrt)); %total weight
    corrtb=corrt; corrtb(corrtb~=0)=1; %binarize
    totallin(t)=0.5*sum(sum(corrtb)); %total number of links

    totalwdacc(t)=0.5*sum(sum(corrtdacc)); %
    corrtbdacc=corrt; corrtbdacc(corrtbdacc~=0)=1; %
    totallindacc(t)=0.5*sum(sum(corrtbdacc)); %dacc
    totalwamy(t)=0.5*sum(sum(corrtamy)); %
    corrtbamy=corrtamy; corrtbamy(corrtbamy~=0)=1; %
    totallinamy(t)=0.5*sum(sum(corrtbamy)); %amy
    % totalwsi(t)=0.5*sum(sum(corrtsi)); %
    % corrtbsi=corrtsi; corrtbsi(corrtbsi~=0)=1; %
    % totallinsi(t)=0.5*sum(sum(corrtbsi)); %si

    Racross{t}=corrt; %all correlation matrices
    Racrossdacc{t}=corrtdacc; %
    Racrossamy{t}=corrtamy; %
    % Racrosssi{t}=corrtsi; %
end



%% ============================================================
% SECTION 15 — THRESHOLD DYNAMIC NETWORKS
% Applies a top-10-percent threshold to each dynamic correlation matrix.
%% ============================================================
% %set negative/below threshhold values to zero: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:N
    Racross{t}(Racross{t}<prctile(Racross{t},90,"all"))=0; %0; %all correlation matrices
    Racrossdacc{t}(Racrossdacc{t}<prctile(Racrossdacc{t},90,"all"))=0; %0; %
    Racrossamy{t}(Racrossamy{t}<prctile(Racrossamy{t},90,"all"))=0; %0; %
    % Racrosssi{t}(Racrosssi{t}<prctile(Racrosssi{t},90,"all"))=0; %0; %
end


%% ============================================================
% SECTION 16 — DYNAMIC SPECTRAL / GRAPH MEASURES: DACC VS AMY
% Computes Fisher-Z, strength, entropy, Fiedler value, spectral gap,
% and Laplacian energy for dACC and amygdala over time.
%% ============================================================
% === SETTINGS ===
CorrMats = {Racrossdacc, Racrossamy};
regionNames = {'DACC', 'AMY'};
colors = lines(2);
density = 0.1;
nRepeats = 10;              % Number of subsampling repetitions
smoothingWindow = 10;       % For plotting
nRegions = numel(CorrMats);
nTime = numel(CorrMats{1});
rng(1);                     % For reproducibility

% === Determine Minimum Channel Count ===
minCh = Nchannels;
minCh0 = inf;
for r = 1:nRegions
    for t = 1:nTime
        Rt = CorrMats{r}{t};
        if isempty(Rt), continue; end
        minCh0 = min(minCh0, size(Rt,1));
    end
end

% === Initialize Outputs ===
nMeasures = 6;
metricNames = {'Mean Fisher Z', 'Strength (Top 10%)', 'Spectral Entropy', ...
               'Fiedler Value', 'Spectral Gap', 'Laplacian Energy'};
ylabels = {'Z', 'Strength', 'Entropy', 'λ₂', 'Gap', 'Energy'};
results = cell(nRegions, nMeasures);

for r = 1:nRegions
    for m = 1:nMeasures
        results{r,m} = nan(nTime, nRepeats);
    end
end

% === Main Loop ===
for rep = 1:nRepeats
    for r = 1:nRegions
        for t = 1:nTime
            Rt = CorrMats{r}{t};
            if isempty(Rt), continue; end
            nCh = size(Rt,1);
            if nCh < minCh
                warning('Experiment %d has only %d channels (less than minCh = %d)', validExpIDs(EN), nCh, minCh);
            end
            idx = randperm(nCh, minCh);
            R = Rt(idx, idx);
            
            % Zero negative correlations
            % R(R < 0) = 0;
            
            % Ensure symmetry
            R = (R + R') / 2;
            
            % Zero diagonal
            R(logical(eye(minCh))) = 0;

            % 1. Mean Fisher Z
            rv = R(triu(true(minCh),1));
            results{r,1}(t,rep) = mean(atanh(rv));

            % 2. Fixed-Density Strength
            [~, sortedIdx] = sort(rv, 'descend');  % Already abs-positive
            nEdges = max(1, floor(density * numel(rv)));
            results{r,2}(t,rep) = mean(rv(sortedIdx(1:nEdges)));

            % 3. Spectral Entropy
            eigVals = eig(R);   % R is now symmetric and non-negative
            eigVals = abs(eigVals);                 % Safety
            % % % % % % % % % % % % % % % % % % eigVals = eigVals / sum(eigVals + eps); % Normalize
            H = -sum(eigVals .* log(eigVals + eps));
            results{r,3}(t,rep) = H / log(minCh);

            % 4–6. Normalized Laplacian
            A = R;  % Already symmetric, non-negative, diagonal = 0
            degVec = sum(A,2);
            DinvSqrt = diag(1 ./ sqrt(degVec + eps));
            Lnorm = eye(minCh) - DinvSqrt * A * DinvSqrt;
            eigL = sort(real(eig(Lnorm)));  % Take real part in case of numerical noise

            if numel(eigL) >= 2
                results{r,4}(t,rep) = eigL(2);                         % Fiedler
                results{r,5}(t,rep) = eigL(end) - eigL(2);             % Spectral Gap
            end
            results{r,6}(t,rep) = sum((eigL - mean(eigL)).^2);         % Laplacian Energy
        end
    end
end


%% ============================================================
% SECTION 17 — SMOOTH AND SUM DACC / AMY MEASURES
% Averages across subsampling repeats, smooths time series, and sums values.
%% ============================================================
% === Smooth and Sum ===
smoothData = cell(nRegions, nMeasures);
sumData = zeros(nRegions, nMeasures);
for r = 1:nRegions
    for m = 1:nMeasures
        avgVals = mean(results{r,m}, 2, 'omitnan');
        smoothData{r,m} = movmean(avgVals, smoothingWindow, 'omitnan');
        sumData(r,m) = nansum(avgVals);
    end
end


%% ============================================================
% SECTION 18 — PLOT DACC / AMY SPECTRAL TIME SERIES
% Plots smoothed dynamic measures for dACC and amygdala.
%% ============================================================
% === PLOT: SMOOTHED TIME SERIES ===
figure('Name', 'Smoothed Spectral Measures (Subsampled)', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    for r = 1:nRegions
        plot(1:nTime, smoothData{r,m}, 'LineWidth', 1.5, 'Color', colors(r,:));
    end
    title(metricNames{m});
    xlabel('Time'); ylabel(ylabels{m});
    grid on;
    if m == 1
        legend(regionNames, 'Location', 'best');
    end
end

% === PLOT: BAR CHART OF TOTAL VALUES ===
figure('Name', 'Summed Spectral Measures (Unnormalized)', 'Position', [200 200 800 400]);
bar(real(sumData'), 'grouped');  % Use real part to suppress warnings
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Summed Value Across Time');
legend(regionNames, 'Location', 'best');
title('Total Spectral Measure Sum Per Region (Subsampled)');
grid on;



%% ============================================================
% SECTION 19 — TEMPORAL STABILITY OF DYNAMIC NETWORKS
% Computes correlations between consecutive connectivity matrices.
%% ============================================================
% %%%%%%%%
 % correlation coefficients along the off-diagonal
corrmat=zeros(1,N-1); corrmatdacc=zeros(1,N-1); corrmatamy=zeros(1,N-1); %corrmatsi=zeros(1,N-1);

for t=1:N-1
    Rt1=Racross{t}; %rand(5,5); 
    Rt2=Racross{t+1}; %rand(5,5); 
    mat=corrcoef(Rt1,Rt2); 
    corrmat(t)=mat(1,length(mat));

    Rt1dacc=Racrossdacc{t}; %dacc 
    Rt2dacc=Racrossdacc{t+1}; % 
    matdacc=corrcoef(Rt1dacc,Rt2dacc); 
    corrmatdacc(t)=matdacc(1,length(matdacc));
    Rt1amy=Racrossamy{t}; %amy 
    Rt2amy=Racrossamy{t+1}; %
    matamy=corrcoef(Rt1amy,Rt2amy); 
    corrmatamy(t)=matamy(1,length(matamy));
    % Rt1si=Racrosssi{t}; %rand(5,5); 
    % Rt2si=Racrosssi{t+1}; %rand(5,5); 
    % matsi=corrcoef(Rt1si,Rt2si); 
    % corrmatsi(t)=matsi(1,length(matsi));
end

corrmat(isnan(corrmat))=0; corrmatdacc(isnan(corrmatdacc))=0; corrmatamy(isnan(corrmatamy))=0; %corrmatsi(isnan(corrmatsi))=0;

% figure; plot(1:N-1-2,corrmat(1:length(corrmat)-2)); ylabel('Correlation: t and (t-1)'); xlabel('Time intervals');
figure; subplot(4,1,1); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmat(1:length(corrmat))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
hold on; xline(event_cell{224,7}(1)); %{98,7}(1)); % % %start time BCIs
hold on; xline(event_cell{224,8}(end)); %{98,8}(end)); % % %end time BCIs
set(gcf,'Color','white'); title('Entire network');
subplot(4,1,2); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmatdacc(1:length(corrmatdacc))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
hold on; xline(event_cell{224,7}(1)); %{98,7}(1)); % %start time BCIs
hold on; xline(event_cell{224,8}(end)); % {98,8}(end)); %{%end time BCIs
set(gcf,'Color','white'); title('dACC network');
subplot(4,1,3); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmatamy(1:length(corrmatamy))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
hold on; xline(event_cell{224,7}(1)); %{98,7}(1)); % %start time BCIs
hold on; xline(event_cell{224,8}(end)); %{98,8}(end)); %%end time BCIs
set(gcf,'Color','white'); title('Amy network');
% subplot(4,1,4); plot(delta*(startT+2*interval):delta*interval:delta*endT,corrmatsi(1:length(corrmatsi))); ylabel('Correlation: t and (t-1)'); xlabel('Time');
% hold on; xline(event_cell{224,7}(1)); %{98,7}(1)); % %start time BCIs
% hold on; xline(event_cell{224,8}(end)); %{98,8}(end)); % %end time BCIs
% set(gcf,'Color','white'); title('SI network');

% %%%%%%%%
% correlation coefficients along the off-diagonal
corrmat=zeros(1,N-1); corrmatdacc=zeros(1,N-1); corrmatamy=zeros(1,N-1); %corrmatsi=zeros(1,N-1);

for t=1:N-1
    Rt1=Racross{t}; Rt2=Racross{t+1}; 
    mat=corrcoef(Rt1,Rt2); 
    corrmat(t)=mat(1,end);

    Rt1dacc=Racrossdacc{t}; Rt2dacc=Racrossdacc{t+1}; 
    matdacc=corrcoef(Rt1dacc,Rt2dacc); 
    corrmatdacc(t)=matdacc(1,end);

    Rt1amy=Racrossamy{t}; Rt2amy=Racrossamy{t+1}; 
    matamy=corrcoef(Rt1amy,Rt2amy); 
    corrmatamy(t)=matamy(1,end);

    % Rt1si=Racrosssi{t}; Rt2si=Racrosssi{t+1}; 
    % matsi=corrcoef(Rt1si,Rt2si); 
    % corrmatsi(t)=matsi(1,end);
end

corrmat(isnan(corrmat))=0; corrmatdacc(isnan(corrmatdacc))=0; 
corrmatamy(isnan(corrmatamy))=0; %corrmatsi(isnan(corrmatsi))=0;


%% ============================================================
% SECTION 20 — SMOOTH TEMPORAL STABILITY CURVES
% Applies moving-average smoothing to matrix-to-matrix stability traces.
%% ============================================================
% === SMOOTHING ===
window = 10;  % adjust this size if you want more or less smoothing
corrmatSmooth     = movmean(corrmat, window);
corrmatdaccSmooth = movmean(corrmatdacc, window);
corrmatamySmooth  = movmean(corrmatamy, window);
% corrmatsiSmooth   = movmean(corrmatsi, window);

% === PLOTTING ===
timeVec = delta * (startT + 2*interval : interval : endT);

figure;

subplot(4,1,1); 
plot(timeVec, corrmatSmooth, 'k', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('Entire network'); hold on;
xline(event_cell{224,7}(1), 'r--'); 
xline(event_cell{224,8}(end), 'r--');

subplot(4,1,2); 
plot(timeVec, corrmatdaccSmooth, 'b', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('dACC network'); hold on;
xline(event_cell{224,7}(1), 'r--'); 
xline(event_cell{224,8}(end), 'r--');

subplot(4,1,3); 
plot(timeVec, corrmatamySmooth, 'm', 'LineWidth', 1.5); 
ylabel('Corr(t, t+1)'); xlabel('Time');
title('Amy network'); hold on;
xline(event_cell{224,7}(1), 'r--'); 
xline(event_cell{224,8}(end), 'r--');

% subplot(4,1,4); 
% plot(timeVec, corrmatsiSmooth, 'g', 'LineWidth', 1.5); 
% ylabel('Corr(t, t+1)'); xlabel('Time');
% title('SI network'); hold on;
% xline(event_cell{224,7}(1), 'r--'); 
% xline(event_cell{224,8}(end), 'r--');

set(gcf, 'Color', 'white');



%% ============================================================
% SECTION 21 — DYNAMIC SPECTRAL / GRAPH MEASURES: DACC, AMY, ALL
% Repeats the spectral/graph analysis for dACC, amygdala, and full network.
%% ============================================================
% === SETTINGS ===
CorrMats = {Racrossdacc, Racrossamy, Racross};  % Add Racross as third region
regionNames = {'DACC', 'AMY', 'ALL'};
colors = lines(3);
density = 0.1;
nRepeats = 10;
smoothingWindow = 10;
nRegions = numel(CorrMats);
rng(1);  % Reproducibility

% === Minimum Channel Count ===
minCh = Nchannels; %inf;
for r = 1:nRegions
    for t = 1:numel(CorrMats{r})
        Rt = CorrMats{r}{t};
        if isempty(Rt), continue; end
        minCh = Nchannels; %min(minCh, size(Rt,1)); %2*Nchannels; %
    end
end

% === Init ===
nMeasures = 6;
metricNames = {'Mean Fisher Z', 'Strength (Top 10%)', 'Spectral Entropy', ...
               'Fiedler Value', 'Spectral Gap', 'Laplacian Energy'};
ylabels = {'Z', 'Strength', 'Entropy', 'λ₂', 'Gap', 'Energy'};
allResults = cell(nRegions,1);
for r = 1:nRegions
    allResults{r} = cell(1,nMeasures);
    for m = 1:nMeasures
        allResults{r}{m} = nan(numel(CorrMats{r}), nRepeats);
    end
end

% === Main Loop ===
for rep = 1:nRepeats
    for r = 1:nRegions
        for t = 1:numel(CorrMats{r})
            Rt = CorrMats{r}{t};
            if isempty(Rt), continue; end
            nCh = size(Rt,1);
            if nCh < minCh
                warning('Experiment %d has only %d channels (less than minCh = %d)', validExpIDs(EN), nCh, minCh);
            end
            idx = randperm(nCh, minCh);
            R = Rt(idx, idx);
            R = (R + R') / 2;               % Symmetrize
            R(logical(eye(minCh))) = 0;     % Zero diagonal

            % 1. Mean Fisher Z
            rv = R(triu(true(minCh),1));
            allResults{r}{1}(t,rep) = mean(atanh(rv));

            % 2. Fixed-Density Strength
            [~, sortedIdx] = sort(rv, 'descend');
            nEdges = max(1, floor(density * numel(rv)));
            allResults{r}{2}(t,rep) = mean(rv(sortedIdx(1:nEdges)));

            % 3. Spectral Entropy
            eigVals = abs(eig(R));
            H = -sum(eigVals .* log(eigVals + eps));
            allResults{r}{3}(t,rep) = H / log(minCh);

            % 4–6. Normalized Laplacian
            A = R;
            degVec = sum(A,2);
            DinvSqrt = diag(1 ./ sqrt(degVec + eps));
            Lnorm = eye(minCh) - DinvSqrt * A * DinvSqrt;
            eigL = sort(real(eig(Lnorm)));

            if numel(eigL) >= 2
                allResults{r}{4}(t,rep) = eigL(2);                  % Fiedler
                allResults{r}{5}(t,rep) = eigL(end) - eigL(2);      % Spectral Gap
            end
            allResults{r}{6}(t,rep) = sum((eigL - mean(eigL)).^2);  % Energy
        end
    end
end


%% ============================================================
% SECTION 22 — EXTRACT FINAL SUMMARY STATISTICS
% Computes final means and standard deviations returned by the function.
%% ============================================================
% === Smoothed & Summary Stats ===
extractSmoothedAndStats = @(data) deal( ...
    movmean(mean(data, 2, 'omitnan'), 10, 'omitnan'), ...
    mean(data, 2, 'omitnan'), ...
    std(data, 0, 2, 'omitnan'));

sm1 = cell(1, nMeasures); mu1 = cell(1, nMeasures); sigma1 = cell(1, nMeasures);
sm2 = cell(1, nMeasures); mu2 = cell(1, nMeasures); sigma2 = cell(1, nMeasures);
sm3 = cell(1, nMeasures); mu3 = cell(1, nMeasures); sigma3 = cell(1, nMeasures);
for m = 1:nMeasures
    [sm1{m}, mu1{m}, sigma1{m}] = extractSmoothedAndStats(allResults{1}{m});
    [sm2{m}, mu2{m}, sigma2{m}] = extractSmoothedAndStats(allResults{2}{m});
    [sm3{m}, mu3{m}, sigma3{m}] = extractSmoothedAndStats(allResults{3}{m});
end


%% ============================================================
% SECTION 23 — FINAL SUMMARY FIGURES
% Plots regional and whole-network summaries.
%% ============================================================
% === FIGURE 1: DACC and AMY ===
figure('Name', 'Spectral Measures: DACC vs AMY', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    fillBetween(1:numel(sm1{m}), mu1{m}, sigma1{m}, colors(1,:), 0.2);
    fillBetween(1:numel(sm2{m}), mu2{m}, sigma2{m}, colors(2,:), 0.2);
    plot(sm1{m}, 'Color', colors(1,:), 'LineWidth', 1.5);
    plot(sm2{m}, 'Color', colors(2,:), 'LineWidth', 1.5);
    title(metricNames{m});
    ylabel(ylabels{m}); xlabel('Time');
    grid on;
    if m==1, legend(regionNames(1:2), 'Location', 'best'); end
end

% === BAR PLOT: DACC + AMY ===
figure('Name', 'Summary: DACC vs AMY', 'Position', [200 200 800 400]);
muMat = [cellfun(@nanmean, mu1); cellfun(@nanmean, mu2)];
stdMat = [cellfun(@nanstd, mu1); cellfun(@nanstd, mu2)];
b = bar(muMat'); hold on;
for r = 1:2, b(r).FaceColor = colors(r,:); end
ngroups = nMeasures; nbars = 2;
groupwidth = min(0.8, nbars/(nbars + 1.5));
x = (1:ngroups)';
for r = 1:2
    errorbar(x - groupwidth/2 + (r-1)*groupwidth/(nbars-1), ...
        muMat(r,:), stdMat(r,:), '.k', 'CapSize', 4);
end
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Mean ± STD'); title('Total Spectral Measure Summary');
legend(regionNames(1:2), 'Location', 'best'); grid on;

% === FIGURE 2: ALL (Racross) ===
figure('Name', 'Spectral Measures: ALL (Racross)', 'Position', [100 100 1200 800]);
for m = 1:nMeasures
    subplot(3,2,m); hold on;
    fillBetween(1:numel(sm3{m}), mu3{m}, sigma3{m}, colors(3,:), 0.2);
    plot(sm3{m}, 'Color', colors(3,:), 'LineWidth', 1.5);
    title(metricNames{m});
    ylabel(ylabels{m}); xlabel('Time');
    grid on;
    if m==1, legend(regionNames(3), 'Location', 'best'); end
end

% === BAR PLOT: ALL ===
figure('Name', 'Summary: ALL (Racross)', 'Position', [200 200 800 400]);
muAll = cellfun(@nanmean, mu3);
stdAll = cellfun(@nanstd, mu3);
bar(muAll, 'FaceColor', colors(3,:)); hold on;
errorbar(1:nMeasures, muAll, stdAll, '.k', 'CapSize', 4);
set(gca, 'XTickLabel', metricNames, 'XTickLabelRotation', 45);
ylabel('Mean ± STD'); title('ALL Region Summary'); grid on;


%% ============================================================
% SECTION 24 — RETURN OUTPUTS
% Function outputs are muAll, stdAll, muMat, and stdMat.
%% ============================================================
%return: muAll, stdAll, muMat, stMat
end