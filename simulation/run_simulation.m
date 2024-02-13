% Run simulation to generate metastability and synchrony map
% Determination of the optimal values of k and tau
clear

data_dir='/young/';%data file directory
datafile_name='CC110045-24.mat';% name of data
dataname=datafile_name(1:11);
savedir='/young/';%save directory of simulation file
t_max=80;%maximum time of simulation
dt=0.0001;%time step
sampling=1; %sampling size
sig_n=0.000; %noise amplitude
sig_f=0; %frequncy
frequency_mean = 60;%mean frequency
%% k and tau paramater
k0=.1;%initial value of k
kfin= 50;%end point of k
k_loop=80; %number of iteration

tau0=.1; %initial value of tau
taufin= 50; %end point of k
tau_loop=80; %number of iteration

parameter=[k0,kfin,k_loop,tau0,taufin,tau_loop];
%% Choice of subgraph
option='YA';
%'option'  is a switching variable.
% if option='SR' ==> Model simulation for subgragh consisting SR
% connections only.
% if option='SR_MR' ==> Model simulation for subgragh consisting SR and MR
% connections.
% if option='MR_LR' ==> Model simulation for subgragh consisting MR and LR
% connections.
% if option='LR' ==> Model simulation for subgragh consisting LR
% connections only.
% if option=' YA'/'OA' ==> Model simulation for entire network
% if option=' YA'/'OA' and depending on the dataset  (whether it is a dataset for
% a young subject or an old subject), simulate the model for the intact
% network of the corresponding subject.

%% Model simulation

sim_meta_map(data_dir,datafile_name,savedir,t_max,dt,sampling,sig_n,sig_f,frequency_mean,parameter,option);
%simulated file will save in the savedir

%% Find optimal k and tau from metastability map

load([savedir,dataname,'_',option,'.mat']); %load the simulation file of the subject

[x,y]=find(metaB_map==max(max(S.metaB_map))); %find out the postion of maximum matastability from metastability map 
opt_k=k_range(x); %optimal k
opt_tau=tau_range(y); %optimal tau
