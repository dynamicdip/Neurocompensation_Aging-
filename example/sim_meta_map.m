% Synchrony and Metastability Map simulation for parameters $k$ and $tau$

function sim_meta_map(data_dir,datafile_name,savedir,t_max,dt,sampling,sig_n,sig_f,frequency_mean,parameter,option)
delete(gcp);
load(fullfile(data_dir,datafile_name));
dataname=datafile_name(1:11);

weights = sc;
weights=weights/max(max(weights)) ;
N=size(weights,1);
tract_lengths = track_mean;
FC(1:N+1:N*N) = 0 ;

k0=parameter(1);kfin=parameter(2);k_loop=parameter(3);
tau0=parameter(4);taufin=parameter(5);tau_loop=parameter(6);

delta_k=(kfin-k0)/k_loop;
delta_tau=(taufin-tau0)/tau_loop;

switch option
    case 'SR'
        I= 0<tract_lengths & tract_lengths<65;    SC=weights.*I;  tl=tract_lengths.*I; fc=FC.*I;
    case 'SR_MR'
        I= 0<tract_lengths & tract_lengths<140;    SC=weights.*I;  tl=tract_lengths.*I; fc=FC.*I;
    case 'MR_LR'
        I= 65<tract_lengths;    SC=weights.*I;  tl=tract_lengths.*I; fc=FC.*I;
    case 'LR'
        I= 140<tract_lengths & tract_lengths==140;    SC=weights.*I;  tl=tract_lengths.*I; fc=FC.*I;
    otherwise
        SC=weights; tl=tract_lengths; fc=FC;
       
end
%%
f_dist=zeros(N,1);
corr_map=zeros(k_loop,tau_loop); dist_map=zeros(k_loop,tau_loop);
sync_map=zeros(k_loop,tau_loop); meta_map=zeros(k_loop,tau_loop);
syncB_map=zeros(k_loop,tau_loop); metaB_map=zeros(k_loop,tau_loop);
k=k0;
parfor kk=1:k_loop
    k=k0+(kk-1)*delta_k;
    tau=tau0;
    corr_list=zeros(1,tau_loop);dist_list=zeros(1,tau_loop);
    sync_list=zeros(1,tau_loop);meta_list=zeros(1,tau_loop);
    syncB_list=zeros(1,tau_loop);metaB_list=zeros(1,tau_loop);
    for jj=1:tau_loop
        disp([k,jj])
        [ths] = Joana_Kuramoto_GlobalCoupling(SC,tl,frequency_mean, ...
            sig_f, f_dist,k,tau,t_max,dt,sampling, sig_n);
        ths_new = transpose(ths);
        rn=sin(ths_new);
        ord1=abs(sum(exp(1i*ths)))/N;
        sync=mean(ord1); meta=std(ord1);
        nn=length(ths);	 dtt=1e-3;
        T = nn*dtt;
        B = BOLD(T,rn(:,1));
        BOLDact = zeros(length(B),N);
        BOLDact(:,1) = B;
        for nnew=2:N
            B = BOLD(T,rn(:,nnew));
            BOLDact(:,nnew) = B;
        end
        ds=100;bds_FI=BOLDact(500:ds:end,:);
        Cb=corrcoef(bds_FI);
        Cb(1:N+1:N*N) = 0 ;
        corrcoef(Cb,fc);
        BOLD_act=(BOLDact-mean(BOLDact))./std(BOLDact);
        phase_BOLD=angle(hilbert(BOLD_act));
        ord1B=abs(sum(exp(1i* transpose(phase_BOLD) )))/N;
        syncB=mean(ord1B);
        metaB=std(ord1B);
        corr_dummy=corrcoef(Cb, fc);
        dist_dummy=sqrt(sum(sum(abs(Cb-fc).^2)));
        corr_list(1,jj)= corr_dummy(1,2);
        dist_list(1,jj)= dist_dummy;
        sync_list(1,jj)= sync; meta_list(1,jj)= meta;
        syncB_list(1,jj)= syncB; metaB_list(1,jj)= metaB;
        tau=tau+delta_tau;
    end
    corr_map(kk,:)=corr_list;    dist_map(kk,:)=dist_list;
    sync_map(kk,:)=sync_list;    meta_map(kk,:)=meta_list;
    syncB_map(kk,:)=syncB_list;    metaB_map(kk,:)=metaB_list;
end
k_range=linspace(k0,k0+(k_loop-1)*delta_k,k_loop);
tau_range=linspace(tau0,tau0+(tau_loop-1)*delta_tau,tau_loop);
S.corr_map=corr_map;    S.dist_map=dist_map;
S.sync_map=sync_map;    S.meta_map=meta_map;
S.syncB_map=syncB_map;  S.metaB_map=metaB_map;
S.k_range=k_range;      S.tau_range=tau_range;
S.fre=frequency_mean;
save([savedir,dataname,'_sim.mat'], 'S');

end