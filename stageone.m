function f = stageone(x,pars)

% return fit of Izi MSN neuron to Moyer et al basic data

% Moyer data - injection current and corresponding output
inj = [0.22	0.225 0.23 0.2350 0.2400 0.2450	0.2500	0.2550	0.2600	0.2650	0.2700	0.2750	0.2800	0.2850	0.2900	0.2950	0.3000];
injScaled = round(inj .* 1e3);
unmod = [0.0000	0.0000	0.0000	2.0000	4.0000	4.0000	6.0000	6.0000	8.0000	8.0000	10.0000	10.0000	12.0000	12.0000	14.0000	14.0000	16.0000];

% pars = [k,a, b, c, vr, vt, vpeak, T, dt, f_start, f_end, f_time,Istore];
k = pars(1); a = pars(2); b = pars(3); c = pars(4); vr = pars(5);  vpeak = pars(6);
T = pars(7); dt = pars(8); f_start=pars(9); f_end=pars(10); f_time=pars(11);Istore = pars(12);
% x = [C, vt, d]
C = x(1); vt =x(2); d = x(3); 

% run simulation
t = 0:dt:T;         
n = length(t); % number of time points
nInj = numel(inj);

for loop = 1:nInj
    v = vr*ones(1,n); u=0*v;
    for i = 1:n-1
        %--- unmodulated
        v(i+1) = v(i) + dt*(k*(v(i)-vr)*(v(i)-vt)-u(i) + injScaled(loop))/C;
        u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));
        % spikes?   
        if v(i+1)>=vpeak
            v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+d;
        end
    end
    % time to first spike
    temp = find(v == vpeak); isis = diff(temp)*dt;
    if temp tfs(loop) = temp(1) * dt; else tfs(loop) = nan; end   % time in ms 
    fI(loop) = sum(v(f_start:f_end) == vpeak) ./ f_time;
    if isis fI1st(loop) = 1000./isis(1); else fI1st(loop) = 0; end

end

%-------------- compute fit


% %% RMSE of fI
% fIerror = sqrt(sum((fI-unmod).^2)/nInj); % RMSE

norm = unmod; norm(norm == 0) = 1;
fIerror = sum(abs(fI-unmod)./norm)/nInj; % mean summed relative error
f = fIerror;


