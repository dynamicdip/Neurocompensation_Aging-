
clear

data_dir='/home/pc/Documents/research/Matlab_shubham_codes_files/camcan/code2send/young/data/';%data file directory
datafile_name='CC110045-24.mat';% name of data
savedir='/home/pc/Documents/research/Matlab_shubham_codes_files/camcan/code2send/young/';%save directory of simulation file
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
option='Y'; %'option'  is a switching variable. 
% if option='SR' ==> Model simulation for subgragh consisting SR
% connections only.
% if option='SR_MR' ==> Model simulation for subgragh consisting SR and MR
% connections.
% if option='MR_LR' ==> Model simulation for subgragh consisting MR and LR
% connections.
% if option='LR' ==> Model simulation for subgragh consisting LR
% connections only.
% if option=' ' ==> Model simulation for entire network
% if option=' ' and depending on the dataset  (whether it is a dataset for
% a young subject or an old subject), simulate the model for the intact 
% network of the corresponding subject.

%% Model simulation

sim_meta_map(data_dir,datafile_name,savedir,t_max,dt,sampling,sig_n,sig_f,frequency_mean,parameter,option);
%simulated file will save in the savedir

