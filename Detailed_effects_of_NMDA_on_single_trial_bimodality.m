%%%% script to illustrate exactly what happens during bimodal "jumps"
%%%% induced by large NMDA conductance fluctuations
%%%% Mark Humphries 23/9/2009

clear all

% found values
load fit_model_results_NEWtuning

rand('state',1); randn('state',1);
% -------------------------------------------------------------------------
% Input parameters
% spike-train parameters: start from 500 events/s
% step to range of spike events/s
N_nmda = 84; alpha_nmda = 0;
N_ampa = 84; alpha_ampa = 0;
N_gaba = 84; alpha_gaba = 0; 

%%% the parameter values for the example histogram/spike-trains (though
%%% this spike-train will be different of course]
mNMDA = 100;
r_nmda = 4; 
r_ampa = 4;
r_gaba = 4;

%% or on a huge scale
% mNMDA = 1;
% N_nmda = 10000;
% r_nmda = 1;

% dopamine levels
D1 = 0;
D2 = 0;

% -------------------------------------------------------------------------
% all PSP parameters in saved file
Egaba = -60; 
Enmda = 0;
Eampa = 0;

% these should stay in the same ratio
PSPampa = Xsyn; %% loaded from file...
PSPnmda = PSPampa / ampa_nmda; PSPgaba = PSPampa ./ ampa_gaba;

% MS neuron parameters in saved file
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

% found MS parameters: X = [C,vt,d]
C = X(1); vt =X(2); d = X(3);

% extra DA model parameters in saved file
KIR = XD1(1);    % KIR modifier 
LCA = XD1(2);    % LCA modifier
vrD1 = vr*(1+D1*KIR);
dD1 = d*(1-D1*LCA);

% D2 - intrinsic
alpha = XD2;
kD2 = k*(1-alpha*D2);

% synaptic
cD1 = Xd1all;
cD2 = Xd2all;

% simulation parameters
T = 1000; % max duration of simulation (milliseconds)
dt = 0.1; % time step 

% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);


% scale NMDA conductance
PSPnmda = mNMDA * PSPampa / ampa_nmda;

Ggaba = zeros(1,n); Gampa = zeros(1,n);Gnmda = zeros(1,n);
Igaba = zeros(1,n); Iampa = zeros(1,n); Inmda = zeros(1,n);
BD1all_nmda = zeros(1,n);
vD1all = vr*ones(1,n); uD1all=0*v;

% generate the spike trains
Sampa = spkgen([0:dt:T], N_ampa, r_ampa, alpha_ampa);
Snmda = spkgen([0:dt:T], N_nmda, r_nmda, alpha_nmda);
Sgaba = spkgen([0:dt:T], N_gaba, r_gaba, alpha_gaba);       
S = sum(Sampa + Snmda + Sgaba);   % total spike-events


% do simulation
for i = 1:n-1
    % NMDA gate 
    BD1all_nmda(i+1)  = 1 ./ (1 + (Mg/3.57) * exp(-vD1all(i)*0.062));    % from Moyer et al 
    
    if Snmda(i) >= 1
        i
    end
    % delete events
%     if i == 9709 % 9706
%         Sampa(i) = 0;
%     end
%     if i >= 9650
%         Snmda(i) = 0;
%     end
% 
    Gampa(i+1) = Gampa(i) + (PSPampa .* Sampa(i)./ts_ampa);
    Gampa(i+1) = Gampa(i+1) * SynExp_ampa;

    Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Snmda(i)./ts_nmda);
    Gnmda(i+1) = Gnmda(i+1) * SynExp_nmda;

    Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i)./ ts_gaba); % add the MS PSPs
    Ggaba(i+1) = Ggaba(i+1) * SynExp_gaba;

    % D1 intrinsic + synaptic: 
    
    Iampa(i+1) = (Gampa(i+1) .* (Eampa - vD1all(i)));
    Inmda(i+1) = (1+cD1*D1)*BD1all_nmda(i+1)*(Gnmda(i+1) .* (Enmda - vD1all(i)));
    Igaba(i+1) = (Ggaba(i+1) .* (Egaba - vD1all(i)));
    
    % standard model
    vD1all(i+1) = vD1all(i) + dt*(k*(vD1all(i)-vrD1)*(vD1all(i)-vt)-uD1all(i) + ...
        Iampa(i+1) + Inmda(i+1) + Igaba(i+1))/C;
    uD1all(i+1) = uD1all(i) + dt*a*(b*(vD1all(i)-vrD1)-uD1all(i));
    % spikes?   
    if vD1all(i+1)>=vpeak
        vD1all(i)=vpeak; vD1all(i+1)=c; 
        uD1all(i+1)=uD1all(i+1)+dD1;
    end

end
        

%--------------------------------------------------------------------------
% plot voltage trace 
figure(1); clf
plot(t,vD1all)

spkTs = t(vD1all == vpeak);
%%% pick a first spike to look at... e.g. 975, 1187
tspk = find(spkTs <= 975,1,'last');
ix = spkTs(tspk) / dt;

ix = 9741;

%%% take 100 ms prior to this of v, and all S
seg = 15 / dt; 
vseg = vD1all(ix-seg:ix+1);
Bseg = BD1all_nmda(ix-seg:ix);
Sgaba_seg = Sgaba(ix-seg:ix); gabaix = find(Sgaba_seg >= 1);
Snmda_seg = Snmda(ix-seg:ix); nmdaix = find(Snmda_seg >= 1);
Sampa_seg = Sampa(ix-seg:ix); ampaix = find(Sampa_seg >= 1);
Igaba_seg = Igaba(ix-seg:ix); 
Inmda_seg = Inmda(ix-seg:ix);
Iampa_seg = Iampa(ix-seg:ix);

figure(2); clf
subplot(411), plot(vseg)
subplot(412), plot(Inmda_seg), hold on, plot(nmdaix,ones(numel(nmdaix),1)*300,'r.')
subplot(413), plot(Bseg)

subplot(414), plot(Iampa_seg+Igaba_seg), hold on, plot(ampaix,ones(numel(ampaix),1),'r.'); 
plot(gabaix,ones(numel(gabaix),1),'b.');

% save all 
fname = ['Detailed_single_trial_bimodality_test_D1_' num2str(D1) '.mat'];
save(fname)

