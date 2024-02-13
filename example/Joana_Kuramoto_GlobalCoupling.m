function [ths] = Joana_Kuramoto_GlobalCoupling(C,D,frequency_mean,f_std,f_dist,k,tau,t_max,dt,sampling, sig_n)

% Network simulation of Kuramoto oscillators with time delays 
% 
% C                    = Matrix of coupling weights (NxN) between pairs of
%                        regions (can be directed and weighted) 
% D                    = Matrix of distances (in mm) (NxN)
% k                    = mean coupling strength
% tau                  = mean delay
% frequency_mean       = Neural populations' average intrinsic frequency
% (Hz) [e.g. 40 Hz]
% f_std                = Standard deviation of intrinsic frequencies across regions
%                        Can be 0 if all equal oscillators.
% f_dist               = Vector containing the distribution of intrinsic frequencies.
%                             if all equal: f_dist=zeros(N,1)
%                             if normally distributed: f_dist=randn(N,1) (saved before running, otherwise every run would result in different intrinsic frequencies)
% t_max                = Total time of simulated activity (seconds) (e.g. 300 s to compare with resting state, for finding parameters)
% dt                   = integration step (smaller than smaller delays) ( in seconds)
%                        (eg. 0.0001 s = 0.1 ms)
% sampling             = sampling for saving simulated activity (eg. 10)
%                        => 10*dt  = 1ms (the actual resolution being saved)
% sig_n                = standard deviation of noise in the phases (can be zero, in radians). Added on each 
%                         node (gaussion, extraction from a distribution in each iteration)
%
%
%    Joana Cabral, 2011
%    
%    Cabral J, Hugues E, Sporns O, Deco G.
%    Role of local network oscillations in resting-state functional connectivity.
%    Neuroimage. 2011 Jul 1;57(1):130-9. Epub 2011 Apr 12.
%
%    Cabral J, Kringelbach ML, Deco G
%    Functional Graph Alterations in Schizophrenia: A Result from a Global Anatomic Decoupling?
%    Pharmacopsychiatry 45 (1), 57
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_t_total   = ceil(t_max/dt);          % total number of time steps e.g 300/0.0001 = 3000000
n_ts  = ceil(n_t_total/sampling);      % same with downsampling ex. 300000
n_tcycle  = 10000;                     % number of time steps in each cycle (for RAM memory reasons)
n_tcs = n_tcycle/sampling;             % same with downsampling ex. 1000
n_cycles   = ceil(n_t_total/n_tcycle); % number of cycles ex. 300

N          = size(C,1); % number of oscillatores (e.g. 90)
C          = k*(C/mean(C(:)));     % Normalize so that mean(C(:))=k; Thus a frontal increase will automatically cause a descrease in the rest of the brain. Be aware of this in gambling study. Might normalize and then increase frontral connection
C          = dt*C;                 % Scale the coupling strengths per time step
I          = C>0;            %
d_m        = mean(D(I));           % mean distance of existing links
v          = d_m/tau;              % conduction velocity (m/s) (speed of transmission) (without delays everything will synchronize over time)
%stepsDelay = round(D/(v*dt*1e3));  % number of time steps for delays
stepsDelay = round(D/(v*dt*1e3));  % number of time steps for delays
sig_noise  = sig_n*sqrt(dt);       % Scale noise per time step

f_diff     = f_dist*f_std;           % define intrinsinc node frequencies.
omega_step = 2*pi*frequency_mean*dt; % 0.0251 radians if f=40Hz and dt = 0.0001. (phase step pr. time step)
omega_diff = 2*pi*f_diff*dt;
omegas     = omega_step+omega_diff;  % Natural phase increment at each time step


n_td = fix(max(stepsDelay(:)))+10; % number of time steps for maximal delays (array for storing old phase values e.g. 90 x n_td)
n_tp = n_td+n_tcycle;              % number of time steps in one cycle with the time for delays (total time)

th   = zeros(N,n_tp,'single');     % initialize phase timeseries for one cycle (e.g. 90 x th)
ths  = zeros(N,n_ts,'single');     % initialize phase timeseries to save (1/10th of thre previous) OK as long as time steps are small enough to avoid aliasing

% Initialization (init without coupling (which would require values that has not yet been calculated)
% Robust in the long run  to different initilisation

    rng('default')
    th(:,1) = 2*pi*rand(N,1);  %<<<<<<<<<<<<INITIAL CONDITION
    for n=1:N % loop through regions
        th(n,1:n_td) = th(n,1)+(0:omegas(n):(n_td-1)*omegas(n));
        th(n,1:n_td) = mod(th(n,1:n_td),2*pi);
    end

% Equations integration

tic

for c = 1:n_cycles
    %disp(['Cycle = ' num2str(c) ' of ' num2str(n_cycles)])
    th(:,n_td+1:n_tp) = 0; % clean all future values (matrix is reused)

    if c < n_cycles      % total number os steps in this cycle
        n_tpc = n_tp;    % normal nr of steps 10000
    else
        n_tpc = n_t_total-(n_cycles-1)*n_tcycle+n_td; % nr of steps to complete total time
    end
    
    for t = n_td:n_tpc-1
        dth = omegas + sig_noise*randn(N,1); % outside summing 
        for n = 1:N
            for p = 1:N
                if C(n,p)>0 % if coupled
                   %dth(n) = dth(n) + C(n,p)*sin(th(p,t-stepsDelay(n,p))-th(n,t)); % main equation
                   dth(n) = dth(n) + C(n,p)*sin(th(p,t-stepsDelay(n,p))-th(n,t)); % main equation
                end
            end
        end
        th(:,t+1) = mod(th(:,t)+dth,2*pi); % update phases 
    end
    ni = (c-1)*n_tcs;
    ns = ceil((n_tpc-n_td)/sampling);
    ths(:,ni+1:ni+ns) = th(:,n_td+1:sampling:n_tpc);    %<<<<<<< Uncomment to show values during calculation
    
    th(:,1:n_td)      = th(:,n_tp-n_td+1:n_tp);
    T_e = toc;
   % disp(['   (Selapsed time = ',num2str(T_e),' s)'])
end

ths(:,1:20/(dt*sampling))=[];  % remove initial 60 seconds of simulations to exclude transient dynamics.

