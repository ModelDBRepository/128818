%%%% script to show how D2 synaptic pars affect output...


clear all

% load tuned model
load fit_model_results_NEWtuning

% -------------------------------------------------------------------------
% Input parameters
% spike-train parameters
N_ctx = 84; alpha_ctx = 0;  r_ctx = [4:0.5:8]; 
N_gaba = 84; alpha_gaba = 0; r_gaba = r_ctx;
nHz = numel(r_ctx);

% dopamine levels
D2 = 0:0.2:1;

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
alpha = XD2(1);  
cD2 = Xd2all;

% simulation parameters
nbatches = 5;
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step

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
ffD2_DA = zeros(nHz,numel(D2));
ffISID2_DA = zeros(nHz,numel(D2));
brate = zeros(nHz,numel(D2));

% run sims...
for bloop = 1:nbatches
    for j = 1:numel(D2)
        j
        kD2 = k * (1-alpha*D2(j));
        for loop = 1:nHz
            Ggaba = zeros(1,n);
            Gampa = zeros(1,n);
            Gnmda = zeros(1,n);
            vD2all = vr*ones(1,n); uD2all=0*v;

            % generate the spike trains
            Sctx = spkgen([0:dt:T], N_ctx, r_ctx(loop), alpha_ctx);
            Sgaba = spkgen([0:dt:T], N_gaba, r_gaba(loop), alpha_gaba);
            S = Sctx + Sgaba;

            for i = 1:n-1
                Gampa(i+1) = Gampa(i) + (PSPampa .* Sctx(i)./ts_ampa);
                Gampa(i+1) = Gampa(i+1) * SynExp_ampa;

                Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Sctx(i)./ts_nmda);
                Gnmda(i+1) = Gnmda(i+1) * SynExp_nmda;

                Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i)./ ts_gaba); % add the MS PSPs
                Ggaba(i+1) = Ggaba(i+1) * SynExp_gaba;

                
                % D2 intrinsic + synaptic: affects AMPA only...
                BD2all_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD2all(i)*0.062));   
                vD2all(i+1) = vD2all(i) + dt*(kD2*(vD2all(i)-vr)*(vD2all(i)-vt)-uD2all(i) + ...
                    (1-cD2*D2(j))*(Gampa(i+1) .* (Eampa - vD2all(i)))+ BD2all_nmda*(Gnmda(i+1) .* (Enmda - vD2all(i))) + (Ggaba(i+1) .* (Egaba - vD2all(i))))/C;

                uD2all(i+1) = uD2all(i) + dt*a*(b*(vD2all(i)-vr)-uD2all(i));
                % spikes?   
                if vD2all(i+1)>=vpeak
                    vD2all(i)=vpeak; vD2all(i+1)=c; 
                    uD2all(i+1)=uD2all(i+1)+d;
                end
            end
            brate(loop,j) = brate(loop,j) + (sum(S) / (T*1e-3)) / nbatches;
            ffD2_DA(loop,j) = ffD2_DA(loop,j) + (sum(vD2all(f_start:f_end) >= vpeak) ./ f_time) / nbatches;
            temp = find(vD2all == vpeak); isis = diff(temp)*dt;
            if isis rate =  1000./mean(isis); else rate = 0; end
            ffISID2_DA (loop,j) = ffISID2_DA (loop,j) + rate / nbatches;
            
        end
    end
end

%--------------------------------------------------------------------------
% plot results - f-f curves
figure(8); clf; hold on
plot(brate,ffD2_DA); 
plot(brate,ffISID2_DA,':'); 

strLgnd = cell(numel(D2),1);
for i = 1:numel(D2)
    strLgnd{i} = num2str(D2(i));
end
legend(strLgnd,'Location','Best')

save ffD2_effect