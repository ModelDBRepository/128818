%%%% script to see if increasing NMDA conductance creates
%%%% bimodality within a single trial (cf Durstewitz & Gabriel, 2007) 


clear all

% found values
load fit_model_results_NEWtuning

rand('state',1); randn('state',1);
% -------------------------------------------------------------------------
% Input parameters
% spike-train parameters:
N_nmda = 84; alpha_nmda = 0;
N_ampa = 84; alpha_ampa = 0;
N_gaba = 84; alpha_gaba = 0; 

r_nmda = 1:1:4; 
r_ampa = 4;
r_gaba = 4;

% dopamine levels
D1 = 0;
D2 = 0;

% NMDA multiplier
mNMDA = [50:50:150];
mGABA = 1;   

% or pick single values to look at spike-trains
% mNMDA = 200;
% r_nmda = 4; 
% r_ampa = 4;
% r_gaba = 4;

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
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step

% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);

% storage
nHz = numel(r_nmda);
vnospks = cell(nHz,numel(mNMDA));
S = zeros(nHz,numel(mNMDA));    % mean spike-event count for each combo
nspks = zeros(nHz,numel(mNMDA)); 
stbl = zeros(nHz,numel(mNMDA)); 
subsamp = floor(1:1/dt:n);
% run sims...

for j = 1:numel(mNMDA)
    j
    
    PSPnmda = mNMDA(j) * PSPampa / ampa_nmda;
    
    PSPgaba = mGABA * PSPgaba / ampa_gaba;
       
    tic
    for loop = 1:nHz
        loop
        Ggaba = zeros(1,n);
        Gampa = zeros(1,n);
        Gnmda = zeros(1,n);
        vD1all = vr*ones(1,n); uD1all=0*v;

        % generate the spike trains
        Sampa = spkgen([0:dt:T], N_ampa, r_ampa, alpha_ampa);
        Snmda = spkgen([0:dt:T], N_nmda, r_nmda(loop), alpha_nmda);
        Sgaba = spkgen([0:dt:T], N_gaba, r_gaba, alpha_gaba);       
        S(loop,j) = sum(Sampa + Snmda + Sgaba);   % total spike-events


        % do simulation
        for i = 1:n-1
            Gampa(i+1) = Gampa(i) + (PSPampa .* Sampa(i)./ts_ampa);
            Gampa(i+1) = Gampa(i+1) * SynExp_ampa;

            Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Snmda(i)./ts_nmda);
            Gnmda(i+1) = Gnmda(i+1) * SynExp_nmda;

            Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i)./ ts_gaba); % add the MS PSPs
            Ggaba(i+1) = Ggaba(i+1) * SynExp_gaba;

            % D1 intrinsic + synaptic
            BD1all_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD1all(i)*0.062));    % from Moyer et al 
            
            % standard model
            vD1all(i+1) = vD1all(i) + dt*(k*(vD1all(i)-vrD1)*(vD1all(i)-vt)-uD1all(i) + ...
                    (Gampa(i+1) .* (Eampa - vD1all(i)))+ (1+cD1*D1)*BD1all_nmda*(Gnmda(i+1) .* (Enmda - vD1all(i))) + (Ggaba(i+1) .* (Egaba - vD1all(i))))/C;
           
            uD1all(i+1) = uD1all(i) + dt*a*(b*(vD1all(i)-vrD1)-uD1all(i));
            % spikes?   
            if vD1all(i+1)>=vpeak
                vD1all(i)=vpeak; vD1all(i+1)=c; 
                uD1all(i+1)=uD1all(i+1)+dD1;
            end

        end
        
        % record all
        vrecord{loop,j} = vD1all;

        % record spikes
        nspks(loop,j) = sum(vD1all == vpeak);
        
    end
    toc
    
end


%--------------------------------------------------------------------------
% plot results 
figure(1); clf
plot(vD1all)

figure(3); clf
edges = -90:1:-20;
mids = edges(1:end-1) + diff(edges)/2;

counter = 0;
for j = 1:numel(mNMDA)
    for loop = 1:nHz
        counter = counter + 1;
        vsnip = vrecord{loop,j}(t>=1000 & t<=t(end));
        vnospks = vsnip(vsnip < -20);
        N = histc(vnospks,edges)';
        subplot(numel(mNMDA),nHz,counter), bar(edges,N,'histc')
        title(['V for ' num2str(mNMDA(j)) 'xNMDA, r=' num2str(r_nmda(loop)) ' NMDA Hz']);
    end    
end

fname = ['NMDA_NEW_single_trial_bimodality_test_D1_' num2str(D1) '.mat'];

if nHz > 1
    % do not save if just doing single run...
    save(fname)
end
