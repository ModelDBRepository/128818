%%% script to run all stages of fitting Izhikevich MSN neurons to data from
%%% Moyer et al's multi-compartment MSN neuron

%% All X0 values taken from hand-tunings in handtune_Iz_MSneuron.m
clear all

% simulation parameters
T = 5000;      % length of simulation in ms
dt = 0.1;      % time step
f_start = 1000/dt;  % start of averaging period
f_end = T/dt;       % end of averaging period
f_time = (f_end - f_start) * 1e-3 * dt; % length of averaging period in seconds
Istore = 270;

%% dopamine levels - saturated dopamine?
D1 = 0.8;
D2 = 0.8;

% Izi neuron parameters: C, a, b, c, d, k, vt, vr, vpeak
% these parameters of Izi MSN neuron are fixed
% pars = [k,a,b,c,vr,vpeak];
izipars = [1,0.01,-20,-55,-80,40];

% synaptic parameters - ratios of ampa:x conductance
ampa_nmda = 2;
ampa_gaba = 1.4;
Egaba = -60;
Enmda = 0;
Eampa = 0;
Mg = 1;         % magnesium ion concentration in mM
ts_ampa  = 6;    % time constant of the glutamate AMPA PSPs
ts_nmda  = 160;    % time constant of the glutamate NMDA PSPs
ts_gaba  = 4;   % time constant of the GABA PSPs

synpars = [ampa_nmda, ampa_gaba, Egaba, Enmda, Eampa, Mg, ts_ampa, ts_nmda, ts_gaba];

%--------------------------------------------------------------------------
%% step 1: tune parameters of basic Izi MSN to fit f-I curve and
%% time-to-first spike of unmodified Moyer neuron
%% search on x = [C,vt,d];

disp('Stage 1')

allpars = [izipars T dt f_start f_end f_time Istore];
X0 = [15,-30,90];    
opts = optimset('fminsearch');
%opts = optimset(opts,'MaxFunEvals',2000);
[X,FVAL,exitflag,output] = fminsearch(@(x) stageone(x,allpars),X0,opts);

%% test resulting model - check time-to-first spike
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);
% found MS parameters: X = [C,vt,d]
C = X(1); vt = X(2); d =X(3); 

v = vr*ones(1,T/dt); u=0*v;
for i = 1:(T/dt-1)
    %--- unmodulated
    v(i+1) = v(i) + dt*(k*(v(i)-vr)*(v(i)-vt)-u(i) + Istore)/C;
    u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));
    % spikes?   
    if v(i+1)>=vpeak;
        v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+d;
    end
end
temp = find(v >= vpeak);
ttfs = temp(1) * dt

%--------------------------------------------------------------------------
%% step 2: take tuned parameters and use as basis for model for D1 and D2
%% intrinsic models - separately tune D1 and D2 intrinsic model parameters
%% to fit corresponding data from Moyer et al
disp('Stage 2a')

allintpars = [allpars X];   % tuned parameters passed to these searches
X0 = [0.03,0.3]; % [mKIR, mLCA]
[XD1,FVAL,exitflag,output] = fminsearch(@(x) stagetwo(x,[allintpars 1 D1]),X0,opts);

disp('Stage 2b')
X0 = 0.04; % alpha
[XD2,FVAL,exitflag,output] = fminsearch(@(x) stagetwo(x,[allintpars 2 D2]),X0,opts);

%--------------------------------------------------------------------------
%% step 3: take tuned parameters from step 1 and use as basis for tuning
%% synaptic conductances to fit f-f curve of unmodified Moyer neuron
%% First: find functional form of f-f curve to fit to! [We do this because
%% f-f data generated from stochastic input, so would be meaningless to fit
%% exactly - but all f-f curves approx linear]
disp('Stage 3')

% regress against data...
% f-f curves
synapticHz = [850.0800	900.1440	950.2080	1000.2720	1050.3360	1100.4000	1150.4640	1200.5280	1250.5920	1300.6560	1350.7200];
unmodff = [0 0.7778	1.4444	2.5556	3.0000	5.0000	5.3333	7.5556	7.8889	9.5556	10.2222];
[B,BINT,R,RINT,STATS] = regress(unmodff',[ones(numel(synapticHz),1) synapticHz']);

% tune parameters - take tuned Izi MSN and add synaptic input
% note - we only explictly tune G_ampa here - G_nmda and G_gaba are both
% fixed ratios of this value, based on ratios given in Moyer et al
X0 = 7;     % size of G_ampa in nS??     
[Xsyn,FVAL,exitflag,output] = fminsearch(@(x) stagethree(x,[allpars X B' synpars]),X0,opts);

%--------------------------------------------------------------------------
%% step 4: take synaptic+basic model tuned parameters from step 3, add D1
%% and D2 intrinsic models from step 2, and separately tune D1 and D2 
%% synaptic effect models to fit f-f cuves from corresponding data from
%% Moyer et al.
%% First: find functional form of f-f curves to fit to!
d1allff = [0.4444	1.2222	3.2222	4.2222	6.3333	8.7778	9.3333	13.1111	13.6667	16.4444	17.3333];
d2allff = [ 0 0 0.4444	0.5556	1.0000	2.0000	2.4444	4.0000	4.5556	6.4444	6.6667];

[Bd1all,BINTd1all,Rd1all,RINTd1all,STATSd1all] = regress(d1allff',[ones(numel(synapticHz),1) synapticHz']);
[Bd2all,BINTd2all,Rd2all,RINTd2all,STATSd2all] = regress(d2allff',[ones(numel(synapticHz),1) synapticHz']);

% takes in: [basic pars, stage 1 pars, fit curve, stage 2 pars, DA type, DA level, stage 3 pars]
% tune parameters - D1 model
disp('Stage 4a')
X0 = 6;  
[Xd1all,FVAL,exitflag,output] = fminsearch(@(x) stagefour(x,[allpars X Bd1all' synpars Xsyn 1 D1 XD1]),X0,opts);

% tune parameters - D2 model
disp('Stage 4b')
X0 = 0.2;
[Xd2all,FVAL,exitflag,output] = fminsearch(@(x) stagefour(x,[allpars X Bd2all' synpars Xsyn 2 D2 XD2]),X0,opts);


%% also: check error in fit for D1 and D2 intrinsic only models to
%% corresponding Moyer et al data!! Nothing to tune, relies entirely on
%% success of previous tunings....
%% First: find functional form of f-f curves to fit to!
d1intff = [0 0.2222	0.5556	1.0000	2.0000	3.6667	4.3333	7.2222	8.1111	11.1111	11.5556];
d2intff = [0.7778	1.5556	2.6667	4.0000	5.5556	7.0000	7.2222	9.7778	10.2222	12.3333	12.5556];

[Bd1int,BINTd1int,Rd1int,RINTd1int,STATSd1int] = regress(d1intff',[ones(numel(synapticHz),1) synapticHz']);
[Bd2int,BINTd2int,Rd2int,RINTd2int,STATSd2int] = regress(d2intff',[ones(numel(synapticHz),1) synapticHz']);

disp('Stage 5 - test')
fD1int = stagefive([allpars X Bd1int' synpars Xsyn 1 D1 XD1]);
fD2int = stagefive([allpars X Bd2int' synpars Xsyn 2 D2 XD2]);

%%% put these parameters back into  script to see fits!!

save fit_model_results_NEWtuning_ISIs
