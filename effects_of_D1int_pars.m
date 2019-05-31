clear all
cl_fig 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% effect of changing D1 on current injection response of tuned MSN


% found values
load fit_model_results_NEWtuning

% MS neuron parameters in saved file
k = izipars(1); a = izipars(2); b = izipars(3); c = izipars(4); vr = izipars(5); vpeak = izipars(6);

% found MS parameters: X = [C,vt,d]
C = X(1); vt =X(2); d = X(3);

% extra DA model parameters in saved file
KIR = XD1(1);    % KIR modifier 
LCA = XD1(2);    % LCA modifier


% parameters for this test
I = 200:5:400;
nInj = numel(I);

D1 = 0:0.2:1;

% simulation parameters
T = 5000; % duration of simulation (milliseconds)
dt = 0.1; % time step in ms

% init simulation 
t = 0:dt:T;
n = length(t); % number of time points
f_start = 1000/dt;
f_end = T/dt;
f_time = (f_end - f_start) * 1e-3 * dt;

%% f-I curves
fID1_DA = zeros(numel(D1),nInj);
fI1stD1_DA = zeros(numel(D1),nInj);
for j = 1:numel(D1)
    for loop = 1:nInj
        loop
        vD1 = vr*ones(1,n); uD1=0*vD1;
        vrD1 = vr*(1+D1(j)*KIR);
        dD1 = d*(1-D1(j)*LCA);
        for i = 1:n-1
            % D1 model
                vD1(i+1) = vD1(i) + dt*(k*(vD1(i)-vrD1)*(vD1(i)-vt)-uD1(i) + I(loop))/C;

                uD1(i+1) = uD1(i) + dt*a*(b*(vD1(i)-vrD1)-uD1(i));
                % spikes?   
                if vD1(i+1)>=vpeak
                    vD1(i)=vpeak; vD1(i+1)=c; 
                    uD1(i+1)=uD1(i+1)+dD1;
                end
        end
        % firing rate at this frequency
        fID1_DA(j,loop) = sum(vD1(f_start:f_end) == vpeak) ./ f_time;
        temp = find(vD1 == vpeak); isis = diff(temp)*dt;
        if isis fI1stD1_DA(j,loop) = 1000./isis(1); else fI1stD1_DA(j,loop) = 0; end
    end
end

figure(2)
plot(I,fID1_DA)
hold on
plot(I,fI1stD1_DA,':')

