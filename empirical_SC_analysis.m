% BAtch analysis of empirical SC

clear;clc;close all
data_dir='/young/';%data file directory

cd(data_dir)
list=dir('*.mat');

s=0; t=0;
TL=[];
for i=1:length(list)
    dataname=list(i).name;
    age(i,1)=str2num(dataname(end-5:end-4));
    load(list(i).name);
    TC=round(sc*1e6);
    sc=sc./max(max(sc));
    TL=track_mean;
    Is= TL>0 & TL<=65; sr=TC.*Is; % SR connections
    Il= TL>=140;lr=TC.*Il; % LR connections
    Im=TL>65 & TL<140;  mr=TC.*Im; % MR connections

    mtl(i)=max(max(TL));
    totalc(i)=sum(sum(TC)); % Total count of exsisting connections
    src(i,1)=(sum(TC(Is)))/2; % Total SR connections
    mrc(i,1)=(sum(TC(Im)))/2; % Total MR connections
    lrc(i,1)=(sum(TC(Il)))/2; % Total LR connections
    nor_src(i)=src(i)/totalc(i); % normalized SR connection
    nor_mrc(i)=mrc(i)/totalc(i); % normalized MR connection
    nor_lrc(i)=lrc(i)/totalc(i); % normalized LR connection
    sr_noRoi(i,1)=length(find(Is))/2;
    lr_noRoi(i,1)=length(find(Il))/2;
    mr_noRoi(i,1)=length(find(Im))/2;
end
