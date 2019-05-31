%%%% script to show effects of D2 intrinsic model

clear all
cl_fig

% found values
% basic tuned model from stage 1....
load fit_model_results_NEWtuning
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

% found MS parameters: X = [C,vt,d]
C = X(1); vt =X(2); d = X(3);

% extra DA model parameters in saved file
alpha = XD2(1);   

% parameters for this test
I = 200:5:400;
nInj = numel(I);

D2 = 0:0.2:1;

% init simulation 
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step
t = 0:dt:T;
n = length(t); % number of time points
f_start = 1000/dt;
f_end = T/dt;
f_time = (f_end - f_start) * 1e-3 * dt;

%% f-I curves
fID2_DA = zeros(numel(D2),nInj);
fI1stD2_DA = zeros(numel(D2),nInj);
for j = 1:numel(D2)
    for loop = 1:nInj
        loop
        vD2 = vr*ones(1,n); uD2=0*vD2;
        kD2 = k * (1-alpha*D2(j));    
        
        for i = 1:n-1
            %--- D2 type                        
            vD2(i+1) = vD2(i) + dt*(kD2*(vD2(i)-vr)*(vD2(i)-vt)-uD2(i) + I(loop))/C;

            uD2(i+1) = uD2(i) + dt*a*(b*(vD2(i)-vr)-uD2(i));
            % spikes?   
            if vD2(i+1)>=vpeak
                vD2(i)=vpeak; vD2(i+1)=c; 
                uD2(i+1)=uD2(i+1)+d;
            end
        end
       % firing rate at this frequency
        fID2_DA(j,loop) = sum(vD2(f_start:f_end) == vpeak) ./ f_time;
        temp = find(vD2 == vpeak); isis = diff(temp)*dt;
        if isis fI1stD2_DA(j,loop) = 1000./isis(1); else fI1stD2_DA(j,loop) = 0; end
        
    end
end

figure(2); clf
plot(I,fID2_DA); hold on
plot(I,fI1stD2_DA,':')

