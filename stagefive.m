function f = stagefive(pars)

% return fit of tuned Izi MSN basic neuron and intrinsic only DA effects to Moyer et al f-f plot data

% approximate r values to get same overall synaptic Hz...
N_ctx = 84; alpha_ctx = 0;  r_ctx = [4:0.5:8]; % r_ctx = 8;  
N_gaba = 84; alpha_gaba = 0; r_gaba = r_ctx; % r_gaba = 8;

% pars = [k,a, b, c, vr, vpeak, T, dt, f_start, f_end, f_time, Istore, C, vt,d];
k = pars(1); a = pars(2); b = pars(3); c = pars(4); vr = pars(5); vpeak = pars(6);
T = pars(7); dt = pars(8); f_start=pars(9); f_end=pars(10); f_time=pars(11); Istore = pars(12);
% newly tuned parameters from Stage 1 added to  pars list
C = pars(13); vt = pars(14); d = pars(15);

% fitted linear function to Moyer data
B1 = pars(16); B2 = pars(17);

% synaptic parameters - fixed values
r_ampa_nmda = pars(18);  r_ampa_gaba = pars(19); 
Egaba = pars(20); Enmda = pars(21); Eampa = pars(22); Mg = pars(23);
ts_ampa = pars(24); ts_nmda =pars(25); ts_gaba = pars(26);

% set conductances  - from Stage 3 search
PSPampa = pars(27);
PSPnmda = PSPampa / r_ampa_nmda;
PSPgaba = PSPampa / r_ampa_gaba;

% set intrinsic DA model parameters - from Stage 2 search
DAtype = pars(28);
DA = pars(29); 
switch DAtype
    case 1
        % last parameters are mKIR and mLCA
        vrD1 = vr*(1+DA*pars(30));
        dD1 = d*(1-DA*pars(31));
    case 2
        % last parameter is alpha
        kD2 = k * (1-pars(30)*DA);
end

% run simulation
t = 0:dt:T;         
n = length(t); % number of time points
nHz = numel(r_ctx);

SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);

for loop = 1:nHz
    v = vr*ones(1,n); u=0*v;
    Ggaba = zeros(1,n);
    Gampa = zeros(1,n);
    Gnmda = zeros(1,n);

    % generate the spike trains
    Sctx = spkgen([0:dt:T], N_ctx, r_ctx(loop), alpha_ctx);
    Sgaba = spkgen([0:dt:T], N_gaba, r_gaba(loop), alpha_gaba);
    S = Sctx + Sgaba;

    for i = 1:n-1
        %%% update synaptic input
        Gampa(i+1) = Gampa(i) + (PSPampa .* Sctx(i)./ts_ampa);
        Gampa(i+1) = Gampa(i+1) * SynExp_ampa;
        Gnmda(i+1) = Gnmda(i) + (PSPnmda .* Sctx(i)./ts_nmda);
        Gnmda(i+1) = Gnmda(i+1) * SynExp_nmda;
        Ggaba(i+1) = Ggaba(i) + (PSPgaba .* Sgaba(i) ./ ts_gaba); 
        Ggaba(i+1) = Ggaba(i+1) * SynExp_gaba;
        
        % update NMDA plug
        B_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-v(i)*0.062));  
        
        %%% switch D1 and D2 type
        switch DAtype
            case 1
                v(i+1) = v(i) + dt*(k*(v(i)-vrD1)*(v(i)-vt)-u(i) + ...
                    (Gampa(i+1) .* (Eampa - v(i)))+ B_nmda*(Gnmda(i+1) .* (Enmda - v(i))) + (Ggaba(i+1) .* (Egaba - v(i))))/C;
                u(i+1) = u(i) + dt*a*(b*(v(i)-vrD1)-u(i));

                % spikes?   
                if v(i+1)>=vpeak
                    v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+dD1;
                end
            case 2
                v(i+1) = v(i) + dt*(kD2*(v(i)-vr)*(v(i)-vt)-u(i) + ...
                    (Gampa(i+1) .* (Eampa - v(i)))+ B_nmda*(Gnmda(i+1) .* (Enmda - v(i))) + (Ggaba(i+1) .* (Egaba - v(i))))/C;       
                u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));

                % spikes?   
                if v(i+1)>=vpeak
                    v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+d;
                end
     
        end

    end
    brate(loop) = sum(S) / (T*1e-3);
    ff(loop) = sum(v(f_start:f_end) == vpeak) ./ f_time;
    temp = find(v == vpeak); isis = diff(temp)*dt;
    if isis ffISI(loop) = 1000./mean(isis); else ffISI(loop) = 0; end
    
end
%-------------- compute fit
% calculate expected value from fitted function
expff = B1 + B2 .* brate;
expff(expff < 0) = 0;   % rectify fit!! 

% f = sum((ff-expff).^2); % SSE
% SRE
norm = expff; norm(norm == 0) = 1;
f = sum(abs(ff-expff)./norm);

