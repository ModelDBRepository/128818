% script to generate time-to-first-spike data for spiking input

clear all
load fit_model_results_NEWtuning


nbatches = 10;

% -------------------------------------------------------------------------
% spike-train parameters
N_ctx = 84; alpha_ctx = 0;  r_ctx = [4:0.5:8]; 
N_gaba = 84; alpha_gaba = 0; r_gaba = r_ctx;

% -------------------------------------------------------------------------
% all PSP parameters in saved file
Egaba = -60;
Enmda = 0;
Eampa = 0;

% these should stay in the same ratio
PSPampa = Xsyn; 
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
kD2 = k * (1-alpha*D2);

% synaptic
cD1 = Xd1all;
cD2 = Xd2all;

% -------------------------------------------------------------------------
% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
f_start = 1000/dt;
f_end = T/dt;
f_time = (f_end - f_start) * 1e-3 * dt;

nHz = numel(r_ctx);
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);

ttfs = zeros(nHz,1); ttfsD1all = zeros(nHz,1); ttfsD2all = zeros(nHz,1);
brate = zeros(nHz,1); frate = zeros(nHz,1); frateD1all = zeros(nHz,1); frateD2all = zeros(nHz,1);
for loop = 1:nHz
    loop
    n_spk = 0; n_D1allspk = 0; n_D2allspk = 0;
    for bloop=1:nbatches
        Ggaba = zeros(1,n);
        Gampa = zeros(1,n);
        Gnmda = zeros(1,n);
        v = vr*ones(1,n); u=0*v;
        vD1all = vr*ones(1,n); uD1all=0*v;
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

            B_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-v(i)*0.062));    

            %%% unmodified
            v(i+1) = v(i) + dt*(k*(v(i)-vr)*(v(i)-vt)-u(i) + (Gampa(i+1) .* (Eampa - v(i))) ...
                              + B_nmda*(Gnmda(i+1) .* (Enmda - v(i))) + (Ggaba(i+1) .* (Egaba - v(i))) )/C;
            u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));
            if v(i+1)>=vpeak
                v(i)=vpeak;
                v(i+1)=c;
                u(i+1)=u(i+1)+d;
            end

            %%% D1 effects        
            % D1 intrinsic + synaptic: affects NMDA only 
            BD1all_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD1all(i)*0.062));   
            vD1all(i+1) = vD1all(i) + dt*(k*(vD1all(i)-vrD1)*(vD1all(i)-vt)-uD1all(i) + ...
                (Gampa(i+1) .* (Eampa - vD1all(i)))+ (1+cD1*D1)*BD1all_nmda*(Gnmda(i+1) .* (Enmda - vD1all(i))) + (Ggaba(i+1) .* (Egaba - vD1all(i))))/C;

            uD1all(i+1) = uD1all(i) + dt*a*(b*(vD1all(i)-vrD1)-uD1all(i));
            % spikes?   
            if vD1all(i+1)>=vpeak
                vD1all(i)=vpeak; vD1all(i+1)=c; 
                uD1all(i+1)=uD1all(i+1)+dD1;
            end

            %%% D2 effects

            % D2 intrinsic + synaptic: affects AMPA only...
            BD2all_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD2all(i)*0.062));   
            vD2all(i+1) = vD2all(i) + dt*(kD2*(vD2all(i)-vr)*(vD2all(i)-vt)-uD2all(i) + ...
                (1-cD2*D2)*(Gampa(i+1) .* (Eampa - vD2all(i)))+ BD2all_nmda*(Gnmda(i+1) .* (Enmda - vD2all(i))) + (Ggaba(i+1) .* (Egaba - vD2all(i))))/C;

            uD2all(i+1) = uD2all(i) + dt*a*(b*(vD2all(i)-vr)-uD2all(i));
            % spikes?   
            if vD2all(i+1)>=vpeak
                vD2all(i)=vpeak; vD2all(i+1)=c; 
                uD2all(i+1)=uD2all(i+1)+d;
            end
        end
        brate(loop) = brate(loop) + (sum(S) / (T*1e-3)) / nbatches;
        frate(loop) = frate(loop) + (sum(v(f_start:f_end) >= vpeak) ./ f_time) / nbatches;
        frateD1all(loop) = frateD1all(loop) + (sum(vD1all(f_start:f_end) >= vpeak) ./ f_time) / nbatches;
        frateD2all(loop) = frateD2all(loop) + (sum(vD2all(f_start:f_end) >= vpeak) ./ f_time) / nbatches;

        % time to first spike
        temp = find(v == vpeak);
        if temp 
            ttfs(loop) = ttfs(loop) + (temp(1) * dt); 
            n_spk = n_spk + 1;
        end

        temp = find(vD1all == vpeak);
        if temp 
            ttfsD1all(loop) = ttfsD1all(loop) + (temp(1) * dt); 
            n_D1allspk = n_D1allspk + 1; 
        end

        temp = find(vD2all == vpeak);
        if temp 
            ttfsD2all(loop) = ttfsD2all(loop) + (temp(1) * dt);
            n_D2allspk = n_D2allspk + 1;  
        end  % time in ms 
    end
    if n_spk ttfs(loop) = ttfs(loop) / n_spk; else ttfs(loop) = nan; end
    if n_D1allspk ttfsD1all(loop) = ttfsD1all(loop) / n_D1allspk; else ttfsD1all(loop) = nan; end
    if n_D2allspk ttfsD2all(loop) = ttfsD2all(loop) / n_D2allspk; else ttfsD2all(loop) = nan; end
end


figure(2); clf
plot(brate,ttfs,'r'); hold on
plot(brate,ttfsD1all,'g');
plot(brate,ttfsD2all,'b');

save TTFS_spike_input_results
