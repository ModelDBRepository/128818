%%% script to show properties of neuron models...

%--------------------------------------------------------------------------
clear all

dt = 0.1; % time step

% dopamine level
D1 = 0.8;
D2 = 0.8;

% original MSN parameters
C = 50; vr = -80;  vt = -25; 
k = 1; a = 0.01; b = -20; 
c = -55; d = 150; vpeak = 40; 

load fit_model_results_NEWtuning

% MS neuron parameters in saved file
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

% for control... increase a...
% a = 0.02;

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


% paired-pulse faciliation (see Mahon paper first for times and
% sizes...)
Ton1 = 50; Toff1 = 250;
ISIs = [200:100:1000];
Ion = 400;
I = 0;
ISIstore = 200;

%% set up sim
deltaT = zeros(numel(ISIs),1); deltaTD1 = zeros(numel(ISIs),1); deltaTD2 = zeros(numel(ISIs),1);
for loop=1:numel(ISIs)
    % set up on and off pulses
    Ton = [Ton1 Toff1+ISIs(loop)];
    Toff = [0 Toff1 Toff1+ISIs(loop)+(Toff1-Ton1)];
    T = Toff(end)+200;
    
    t = dt:dt:T;
    n = length(t); % number of time points
    v = vr*ones(1,n); u=0*v;
    vD1 = vr*ones(1,n); uD1=0*vD1;
    vD2 = vr*ones(1,n); uD2=0*vD2;
    
    
    for i = 1:n-1

        if any(i*dt == Ton)
            I = Ion;
        elseif any(i*dt == Toff)
            I = 0;
        end    
        v(i+1) = v(i) + dt*(k*(v(i)-vr)*(v(i)-vt)-u(i) + I)/C;
        u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));

        if v(i+1)>=vpeak
            v(i)=vpeak;
            v(i+1)=c;
            u(i+1)=u(i+1)+d;
        end
        
        %%% D1 type
        vD1(i+1) = vD1(i) + dt*(k*(vD1(i)-vrD1)*(vD1(i)-vt)-uD1(i) + I)/C;
        uD1(i+1) = uD1(i) + dt*a*(b*(vD1(i)-vrD1)-uD1(i));
      
        % spikes?   
        if vD1(i+1)>=vpeak
            vD1(i)=vpeak; vD1(i+1)=c; 
            uD1(i+1)=uD1(i+1)+dD1;
        end
        
        %--- D2 type                    
        vD2(i+1) = vD2(i) + dt*(kD2*(vD2(i)-vr)*(vD2(i)-vt)-uD2(i) + I)/C;

        uD2(i+1) = uD2(i) + dt*a*(b*(vD2(i)-vr)-uD2(i));
        % spikes?   
        if vD2(i+1)>=vpeak
            vD2(i)=vpeak; vD2(i+1)=c; 
            uD2(i+1)=uD2(i+1)+d;
        end

    end
    
    % spike times for baseline
    spkts = find(v == vpeak)*dt;
    spk1 = spkts(find(spkts < Toff(2),1,'first'));
    spk2 = spkts(find(spkts > Ton(2) & spkts <Toff(end),1,'first'));
    t1 = spk1 - Ton(1); t2 = spk2 - Ton(2);
    if ~isempty(t1)& ~isempty(t2)  deltaT(loop) = t1-t2; end
    
    % for D1 model
    spktsD1 = find(vD1 == vpeak)*dt;
    spk1 = spktsD1(find(spktsD1 < Toff(2),1,'first'));
    spk2 = spktsD1(find(spktsD1 > Ton(2) & spktsD1 <Toff(end),1,'first'));
    t1 = spk1 - Ton(1); t2 = spk2 - Ton(2);
    if ~isempty(t1) & ~isempty(t2) deltaTD1(loop) = t1-t2; end
    
    % for D2 model
    spktsD2 = find(vD2 == vpeak)*dt;
    spk1 = spktsD2(find(spktsD2 < Toff(2),1,'first'));
    spk2 = spktsD2(find(spktsD2 > Ton(2) & spktsD2 <Toff(end),1,'first'));
    t1 = spk1 - Ton(1); t2 = spk2 - Ton(2);
    if ~isempty(t1)& ~isempty(t2) deltaTD2(loop) = t1-t2; end
   
    % draw current pulse on figure
    poff = -95;
    pon = poff + 2*((Ion/1000)/0.1);
    pulse = [ones(Ton(1)/dt,1).*poff; ones((Toff(2)-Ton(1))/dt,1).*pon; ones((Ton(2)-Toff(2))/dt,1).*poff;...
                ones((Toff(3)-Ton(2))/dt,1).*pon; ones((T-Toff(3))/dt,1).*poff];


    figure(loop); clf
    subplot(131), plot(t,v,'k','LineWidth',2) 
    hold on, plot(t,pulse,'LineWidth',2); 
    xlabel('Time (ms)'); ylabel('Membrane potential (mV)')
    subplot(132), plot(t,vD1,'k','LineWidth',2) 
    hold on, plot(t,pulse,'LineWidth',2); 
    xlabel('Time (ms)'); ylabel('Membrane potential (mV)')
    subplot(133), plot(t,vD2,'k','LineWidth',2) 
    hold on, plot(t,pulse,'LineWidth',2); 
    xlabel('Time (ms)'); ylabel('Membrane potential (mV)')
    
    if ISIs(loop) == ISIstore
        baseoutput = [t' v' pulse];
    end
end

figure
plot(ISIs,deltaT,'k+-')
hold on
% plot D1 and D2 effects
plot(ISIs,deltaTD1,'g+-')
plot(ISIs,deltaTD2,'b+-')

save paired_pulse_responses

