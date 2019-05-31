% script to run MS neuron with DA modulation
% and test f-I and f-f curves....

clear all
% -------------------------------------------------------------------------
% Injection current parameters
I = 200:5:300;

% spike train parameters - need to change these later....
N_ctx = 84; alpha_ctx = 0;  r_ctx = [4:0.5:8]; % r_ctx = 8;  
N_gaba = 84; alpha_gaba = 0; r_gaba = r_ctx; % r_gaba = 8;

% dopamine
D1 = 0.8;
D2 = 0.8;

% -------------------------------------------------------------------------
% PSP parameters - all from Moyer et al (2007)
ts_ampa  = 6;    % time constant of the glutamate AMPA PSPs
ts_nmda  = 160;    % time constant of the glutamate NMDA PSPs
ts_gaba  = 4;   % time constant of the GABA PSPs

% these should stay in the same ratio?
PSPampa = 6.5; %6;    %?? 0.6 nS conductance of the cortical synapses
%PSPnmda = 3;    %?? 0.3 nS conductance of the FS synapses
%PSPgaba = 4;     %?? 0.44 nS conductance of the MS synapses
% so: AMPA:NMDA = 2:1; AMPA:GABA = 1.4
PSPnmda = PSPampa / 2; PSPgaba = PSPampa ./ 1.4;

Eampa = 0;       % glutamate reversal potential
Enmda = 0;
Egaba = -60;    % GABA reversal potential
Mg = 1;         % magnesium ion concentration in mM

% -------------------------------------------------------------------------
% simulation parameters
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step
tau = dt;

% -------------------------------------------------------------------------
% Izhikevich MSN parameters
C = 50; vr = -80;  vt = -25; 
k = 1; a = 0.01; b = -20; 
c = -55; d = 150; vpeak = 40; 

% modified MS neuron parameters
C = 15;     
vt = -30;    
d = 90;
 
% D1 model
KIR = 0.03;
LCA = 0.3;
vrD1 = vr*(1+D1*KIR);
dD1 = d*(1-D1*LCA);

% D2 - intrinsic
alpha = 0.04;

% synaptic
cD1 = 6; 
cD2 = 0.2; 

%--------------------------------------------------------------------------
% Moyer's data!

% f-I curves....
inj = [0.22	0.225 0.23 0.2350 0.2400 0.2450	0.2500	0.2550	0.2600	0.2650	0.2700	0.2750	0.2800	0.2850	0.2900	0.2950	0.3000];
injScaled = inj .* 1e3;
unmod = [0.0000	0.0000	0.0000	2.0000	4.0000	4.0000	6.0000	6.0000	8.0000	8.0000	10.0000	10.0000	12.0000	12.0000	14.0000	14.0000	16.0000];
d1int = [0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	2.0000	6.0000	8.0000	10.0000	12.0000	14.0000	16.0000	16.0000	18.0000];
d2int = [2.0000	2.0000	4.0000	6.0000	6.0000	8.0000	8.0000	10.0000	10.0000	12.0000	12.0000	14.0000	14.0000	16.0000	16.0000	18.0000	18.0000];

% f-f curves
synapticHz = [850.0800	900.1440	950.2080	1000.2720	1050.3360	1100.4000	1150.4640	1200.5280	1250.5920	1300.6560	1350.7200];
unmodff = [0 0.7778	1.4444	2.5556	3.0000	5.0000	5.3333	7.5556	7.8889	9.5556	10.2222];
d1intff = [0 0.2222	0.5556	1.0000	2.0000	3.6667	4.3333	7.2222	8.1111	11.1111	11.5556];
d2intff = [0.7778	1.5556	2.6667	4.0000	5.5556	7.0000	7.2222	9.7778	10.2222	12.3333	12.5556];
d1allff = [0.4444	1.2222	3.2222	4.2222	6.3333	8.7778	9.3333	13.1111	13.6667	16.4444	17.3333];
d2allff = [ 0 0 0.4444	0.5556	1.0000	2.0000	2.4444	4.0000	4.5556	6.4444	6.6667];

% -------------------------------------------------------------------------
% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
f_start = 1000/dt;
f_end = T/dt;
f_time = (f_end - f_start) * 1e-3 * dt;

nInj = numel(I);
Istore = 270;

% -------------------------------------------------------------------------
%% GO SIMULATIONS
% f-I curves
%V = repmat(v,3,1);
for loop = 1:nInj
    loop
    v = vr*ones(1,n); u=0*v;
    vD1 = vr*ones(1,n); uD1=0*vD1;
    vD2 = vr*ones(1,n); uD2=0*vD2;

    for i = 1:n-1
        %--- unmodulated
        v(i+1) = v(i) + dt*(k*(v(i)-vr)*(v(i)-vt)-u(i) + I(loop))/C;
        u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));
        % spikes?   
        if v(i+1)>=vpeak
            v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+d;
        end
     
        %%% new D1 model
        vD1(i+1) = vD1(i) + dt*(k*(vD1(i)-vrD1)*(vD1(i)-vt)-uD1(i) + I(loop))/C;
        uD1(i+1) = uD1(i) + dt*a*(b*(vD1(i)-vrD1)-uD1(i));
      
        % spikes?   
        if vD1(i+1)>=vpeak
            vD1(i)=vpeak; vD1(i+1)=c; 
            uD1(i+1)=uD1(i+1)+dD1;
        end
        
        %--- D2 type        
        kD2 = k * (1-alpha*D2);
        
        vD2(i+1) = vD2(i) + dt*(kD2*(vD2(i)-vr)*(vD2(i)-vt)-uD2(i) + I(loop))/C;
        
        uD2(i+1) = uD2(i) + dt*a*(b*(vD2(i)-vr)-uD2(i));
        % spikes?   
        if vD2(i+1)>=vpeak
            vD2(i)=vpeak; vD2(i+1)=c; 
            uD2(i+1)=uD2(i+1)+d;
        end
    end
    % time to first spike
    temp = find(v == vpeak); isis = diff(temp)*dt;
    if temp tfs(loop) = temp(1) * dt; else tfs(loop) = nan; end   % time in ms 
    temp = find(vD1 == vpeak); isisD1 = diff(temp)*dt;
    if temp tfsD1(loop) = temp(1) * dt; else tfsD1(loop) = nan; end   % time in ms 
    temp = find(vD2 == vpeak); isisD2 = diff(temp)*dt;
    if temp tfsD2(loop) = temp(1) * dt; else tfsD2(loop) = nan; end  % time in ms 
    
    % firing rate at this frequency
    fI(loop) = sum(v(f_start:f_end) == vpeak) ./ f_time;
    fID1(loop) = sum(vD1(f_start:f_end) == vpeak) ./ f_time;
    fID2(loop) = sum(vD2(f_start:f_end) == vpeak) ./ f_time;
    
    % instantaneous rate (first ISI)
    if isis fI1st(loop) = 1000./isis(1); else fI1st(loop) = 0; end
    if isisD1 fID11st(loop) = 1000./isisD1(1); else fID11st(loop) = 0; end
    if isisD2 fID21st(loop) = 1000./isisD2(1); else fID21st(loop) = 0; end
    
    % store one set of membrane potential traces
    if I(loop) == Istore
        V = [v; vD1; vD2];
    end
end
% -------------------------------------------------------------------------
% plot results - f-I curves
figure(1); clf; plot(I,fI,'+-'); hold on; plot(I,fID1,'r+-'); plot(I,fID2,'k+-');
plot(I,fI1st,'+--'); plot(I,fID11st,'r+--'); plot(I,fID21st,'k+--');
plot(injScaled,unmod,'+:'); plot(injScaled,d1int,'r+:'); plot(injScaled,d2int,'k+:');
xlabel('Current injection (pA)'); ylabel('Spiking frequency (spikes/s)'); title('Izhikevich MSN f-I curves')
legend('No dopamine','D1 intrinsic','D2 intrinsic','Location','Best')

% plot results - time to first spike
figure(2); clf; plot(I,tfs,'+-'); hold on; plot(I,tfsD1,'r+-'); plot(I,tfsD2,'k+-');
xlabel('Current injection (pA)'); ylabel('Time to first spike (milliseconds)'); title('Izhikevich MSN time-to-first-spike')
legend('No dopamine','D1 intrinsic','D2 intrinsic','Location','Best')

% plot results - membrane trace
figure(3); clf;
subplot(131),plot(t,V(1,:)); xlabel('Time (milliseconds)'); ylabel('Membrane potential (mV)'); title(['Response of no DA MSN to ' num2str(Istore) 'pA injection']);
axis([0 1000 vr vpeak+5])
subplot(132),plot(t,V(2,:)); xlabel('Time (milliseconds)'); ylabel('Membrane potential (mV)'); title(['Response of MSN D1 intrinsic to ' num2str(Istore) 'pA injection']);
axis([0 1000 vr vpeak+5])
subplot(133),plot(t,V(3,:)); xlabel('Time (milliseconds)'); ylabel('Membrane potential (mV)'); title(['Response of MSN D2 intrinsic to ' num2str(Istore) 'pA injection']);
axis([0 1000 vr vpeak+5])

%return

%--------------------------------------------------------------------------
%% f-f curves
nHz = numel(r_ctx);
SynExp_ampa = exp(-dt / ts_ampa);
SynExp_nmda = exp(-dt / ts_nmda);
SynExp_gaba = exp(-dt / ts_gaba);


for loop = 1:nHz
    loop
    Ggaba = zeros(1,n);
    Gampa = zeros(1,n);
    Gnmda = zeros(1,n);
    v = vr*ones(1,n); u=0*v;
    vD1int = vr*ones(1,n); uD1int=0*v;
    vD2int = vr*ones(1,n); uD2int=0*v;
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
        % D1 intrinsic only
        BD1int_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD1int(i)*0.062));    
        vD1int(i+1) = vD1int(i) + dt*(k*(vD1int(i)-vrD1)*(vD1int(i)-vt)-uD1int(i) + ...
            (Gampa(i+1) .* (Eampa - vD1int(i)))+ BD1int_nmda*(Gnmda(i+1) .* (Enmda - vD1int(i))) + (Ggaba(i+1) .* (Egaba - vD1int(i))))/C;
        uD1int(i+1) = uD1int(i) + dt*a*(b*(vD1int(i)-vrD1)-uD1int(i));
        % spikes?   
        if vD1int(i+1)>=vpeak
            vD1int(i)=vpeak; vD1int(i+1)=c; 
            uD1int(i+1)=uD1int(i+1)+dD1;
        end
        
        % D1 intrinsic + synaptic
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
        kD2 = k * (1-alpha*D2);
        
        % D2 intrinsic only
        BD2int_nmda  = 1 ./ (1 + (Mg/3.57) * exp(-vD2int(i)*0.062));   
        vD2int(i+1) = vD2int(i) + dt*(kD2*(vD2int(i)-vr)*(vD2int(i)-vt)-uD2int(i) + ...
            (Gampa(i+1) .* (Eampa - vD2int(i)))+ BD2int_nmda*(Gnmda(i+1) .* (Enmda - vD2int(i))) + (Ggaba(i+1) .* (Egaba - vD2int(i))))/C;
        
        uD2int(i+1) = uD2int(i) + dt*a*(b*(vD2int(i)-vr)-uD2int(i));
        % spikes?   
        if vD2int(i+1)>=vpeak
            vD2int(i)=vpeak; vD2int(i+1)=c; 
            uD2int(i+1)=uD2int(i+1)+d;
        end

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
    brate(loop) = sum(S) / (T*1e-3);
    frate(loop) = sum(v(f_start:f_end) >= vpeak) ./ f_time;
    frateD1int(loop) = sum(vD1int(f_start:f_end) >= vpeak) ./ f_time;
    frateD2int(loop) = sum(vD2int(f_start:f_end) >= vpeak) ./ f_time;
    frateD1all(loop) = sum(vD1all(f_start:f_end) >= vpeak) ./ f_time;
    frateD2all(loop) = sum(vD2all(f_start:f_end) >= vpeak) ./ f_time;
    
    temp = find(v == vpeak); isis = diff(temp)*dt;
    if isis ffISI(loop) = 1000./mean(isis); else ffISI(loop) = 0; end
    temp = find(vD1int == vpeak); isis = diff(temp)*dt;
    if isis ffISID1int (loop) = 1000./mean(isis); else ffISID1int(loop) = 0; end
    temp = find(vD2int == vpeak); isis = diff(temp)*dt;
    if isis ffISID2int (loop) = 1000./mean(isis); else ffISID2int(loop) = 0; end
    temp = find(vD1all == vpeak); isis = diff(temp)*dt;
    if isis ffISID1all (loop) = 1000./mean(isis); else ffISID1all(loop) = 0; end
    temp = find(vD2all == vpeak); isis = diff(temp)*dt;
    if isis ffISID2all (loop) = 1000./mean(isis); else ffISID2all(loop) = 0; end
end

%--------------------------------------------------------------------------
% plot results - f-f curves
%figure(5); clf; bar(t,S); xlabel('Time (milliseconds)'); ylabel('Number of spikes');
figure(6);clf; subplot(131), plot(t,Gampa); title('AMPA'); xlabel('Time (milliseconds)'); ylabel('G_{ampa}');
 subplot(132), plot(t,Gnmda); title('NMDA'); xlabel('Time (milliseconds)'); ylabel('G_{nmda}');
  subplot(133), plot(t,Ggaba); title('GABA'); xlabel('Time (milliseconds)'); ylabel('G_{gaba}');
figure(7); clf;
plot(t,v); xlabel('Time (milliseconds)'); ylabel('Membrane potential (mV)'); title(['Response to total background rate of ' num2str(brate(end)) ' spikes/s']);

figure(8); clf; hold on
plot(brate,frate,'+-'); 
plot(brate,frateD1int,'r+-'); 
plot(brate,frateD2int,'k+-'); 
plot(brate,frateD1all,'m+-'); 
plot(brate,frateD2all,'g+-'); 
axis([800 1400 0 20]);
xlabel('Total synaptic input (spikes/s)'); ylabel('Output (spikes/s)'); title('Izhikevich MSN neuron f-f curves')
legend('Unmodified','D1 intrinsic','D2 intrinsic','D1 all','D2 all','Location','Best')

figure(9); clf; hold on
plot(synapticHz,unmodff,'+:')
plot(synapticHz,d1intff,'r+:');
plot(synapticHz,d2intff,'k+:');
plot(synapticHz,d1allff,'m+:');
plot(synapticHz,d2allff,'g+:');
axis([800 1400 0 20]);
xlabel('Total synaptic input (spikes/s)'); ylabel('Output (spikes/s)'); title('Moyer et al MSN neuron f-f curves')
legend('Unmodified','D1 intrinsic','D2 intrinsic','D1 all','D2 all','Location','Best')

figure(10); clf; hold on
plot(brate,frate,'+-'); 
plot(brate,frateD1int,'r+-'); 
plot(brate,frateD2int,'k+-'); 
plot(brate,frateD1all,'m+-'); 
plot(brate,frateD2all,'g+-'); 
plot(synapticHz,unmodff,'+:')
plot(synapticHz,d1intff,'r+:');
plot(synapticHz,d2intff,'k+:');
plot(synapticHz,d1allff,'m+:');
plot(synapticHz,d2allff,'g+:');
axis([800 1400 0 20]);
xlabel('Total synaptic input (spikes/s)'); ylabel('Output (spikes/s)'); title('Izhikevich MSN neuron f-f curves vs data')
legend('Unmodified model','D1 intrinsic','D2 intrinsic','D1 all model','D2 all model','Location','Best')

%%% same again, but using mean ISI as firing rate
figure(11); clf; hold on
plot(brate,ffISI,'+-'); 
plot(brate,ffISID1int,'r+-'); 
plot(brate,ffISID2int,'k+-'); 
plot(brate,ffISID1all,'m+-'); 
plot(brate,ffISID2all,'g+-'); 
plot(synapticHz,unmodff,'+:')
plot(synapticHz,d1intff,'r+:');
plot(synapticHz,d2intff,'k+:');
plot(synapticHz,d1allff,'m+:');
plot(synapticHz,d2allff,'g+:');
axis([800 1400 0 20]);
xlabel('Total synaptic input (spikes/s)'); ylabel('Output (spikes/s)'); title('Izhikevich MSN neuron f-f (ISI) curves vs data')
legend('Unmodified model','D1 intrinsic','D2 intrinsic','D1 all model','D2 all model','Location','Best')

save handtune_results
