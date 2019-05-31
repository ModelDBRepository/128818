%%%% script to show how D1 synaptic pars affect output...


clear all

% load tuned model
load fit_model_results_NEWtuning

% -------------------------------------------------------------------------
% Input parameters
% spike-train parameters
N_ctx = 84; alpha_ctx = 0;  r_ctx = [4:0.5:8]; 
N_gaba = 84; alpha_gaba = 0; r_gaba = r_ctx;
nHz = numel(r_ctx);
% r_gaba = r_gaba * 0;    % no GABA input

% dopamine levels
D1 = 0:0.2:1;

% -------------------------------------------------------------------------
% all PSP parameters in saved file
Egaba = -60;
Enmda = 0;
Eampa = 0;

% these should stay in the same ratio
PSPampa = Xsyn; %% loaded from file...
PSPnmda = PSPampa / ampa_nmda; PSPgaba = PSPampa ./ ampa_gaba;

% fixed parameters
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

% found MS parameters: X = [C,vt,d]
C = X(1); vt =X(2); d = X(3);

% extra DA model parameters in saved file
KIR = XD1(1);    % KIR modifier 
LCA = XD1(2);    % LCA modifier
vrD1 = vr*(1+D1*KIR);
dD1 = d*(1-D1*LCA);
        
cD1 = Xd1all;   % coefficient for increase in NMDA PSP 


% simulation parameters
nbatches = 5;
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step - was 0.01????

% f-f curves
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);

% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
f_start = 1000/dt;
f_end = T/dt;
f_time = (f_end - f_start) * 1e-3 * dt;
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);

% storage
ffD1int = zeros(nHz,numel(D1));
ffD1_DA = zeros(nHz,numel(D1));
ffISID1_DA = zeros(nHz,numel(D1));
brate = zeros(nHz,numel(D1));

% run sims...
for bloop = 1:nbatches
    for j = 1:numel(D1)
        j
        for loop = 1:nHz
            Ggaba = zeros(1,n);
            Gampa = zeros(1,n);
            Gnmda = zeros(1,n);
            vD1int = vr*ones(1,n); uD1int=0*v;
            vD1all = vr*ones(1,n); uD1all=0*v;

            % generate the spike trains
            Sctx = spkgen([0:dt:T], N_ctx, r_ctx(loop), alpha_ctx);
            Sgaba = spkgen([0:dt:T], N_gaba, r_gaba(loop), alpha_gaba);
            S = Sctx + Sgaba;

            for i = 1:n-1
                Gampa(i+1) = Gampa(i) + (PSPampa .* Sctx(i)./ts_ampa);
                % Gampa(i+1) = Gampa(i) + (PSPampa .* Sctx(i));
                Gampa(i+1) = Gampa(i+1) * SynExp_ampa;

                Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Sctx(i)./ts_nmda);
                % Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Sctx(i));
                Gnmda(i+1) = Gnmda(i+1) * SynExp_nmda;

                Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i)./ ts_gaba); % add the MS PSPs
                % Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i));    % without synaptic scaling
                Ggaba(i+1) = Ggaba(i+1) * SynExp_gaba;

                %%% D1 effects
                % D1 intrinsic only
                BD1int_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD1int(i)*0.062));    % from Moyer et al 
                vD1int(i+1) = vD1int(i) + dt*(k*(vD1int(i)-vrD1(j))*(vD1int(i)-vt)-uD1int(i) + ...
                    (Gampa(i+1) .* (Eampa - vD1int(i)))+ BD1int_nmda*(Gnmda(i+1) .* (Enmda - vD1int(i))) + (Ggaba(i+1) .* (Egaba - vD1int(i))))/C;
                uD1int(i+1) = uD1int(i) + dt*a*(b*(vD1int(i)-vrD1(j))-uD1int(i));
                % spikes?   
                if vD1int(i+1)>=vpeak
                    vD1int(i)=vpeak; vD1int(i+1)=c; 
                    uD1int(i+1)=uD1int(i+1)+dD1(j);
                end

                % D1 intrinsic + synaptic: affects NMDA only according to Moyer et al...
                BD1all_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD1all(i)*0.062));    % from Moyer et al 
                vD1all(i+1) = vD1all(i) + dt*(k*(vD1all(i)-vrD1(j))*(vD1all(i)-vt)-uD1all(i) + ...
                    (Gampa(i+1) .* (Eampa - vD1all(i)))+ (1+cD1*D1(j))*BD1all_nmda*(Gnmda(i+1) .* (Enmda - vD1all(i))) + (Ggaba(i+1) .* (Egaba - vD1all(i))))/C;

                uD1all(i+1) = uD1all(i) + dt*a*(b*(vD1all(i)-vrD1(j))-uD1all(i));
                % spikes?   
                if vD1all(i+1)>=vpeak
                    vD1all(i)=vpeak; vD1all(i+1)=c; 
                    uD1all(i+1)=uD1all(i+1)+dD1(j);
                end

            end
            brate(loop,j) = brate(loop,j) + (sum(S) / (T*1e-3)) / nbatches;
            ffD1int(loop,j) = ffD1int(loop,j) + (sum(vD1int(f_start:f_end) >= vpeak) ./ f_time) / nbatches;
            ffD1_DA(loop,j) = ffD1_DA(loop,j) + (sum(vD1all(f_start:f_end) >= vpeak) ./ f_time) / nbatches;
            temp = find(vD1all == vpeak); isis = diff(temp)*dt;
            if isis rate =  1000./mean(isis); else rate = 0; end
            ffISID1_DA (loop,j) = ffISID1_DA (loop,j) + rate / nbatches;
            
            %keyboard
        end
    end
end

%--------------------------------------------------------------------------
% plot results - f-f curves
figure(8); clf; hold on
plot(brate,ffD1_DA); hold on
plot(brate,ffISID1_DA,':'); 

strLgnd = cell(numel(D1),1);
for i = 1:numel(D1)
    strLgnd{i} = num2str(D1(i));
end
legend(strLgnd,'Location','Best')

figure(9); clf; 
plot(brate,ffD1int); 
strLgnd = cell(numel(D1),1);
for i = 1:numel(D1)
    strLgnd{i} = num2str(D1(i));
end
legend(strLgnd,'Location','Best')

% save ffD1_effect