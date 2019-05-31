function f = stagetwo(x,pars)

% return fit of Izi MSN D1 or D2 type neuron (intrinsic) to corresponding Moyer et al data

% Moyer data - injection current and corresponding output
inj = [0.22	0.225 0.23 0.2350 0.2400 0.2450	0.2500	0.2550	0.2600	0.2650	0.2700	0.2750	0.2800	0.2850	0.2900	0.2950	0.3000];
injScaled = round(inj .* 1e3);
d1int = [0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	2.0000	6.0000	8.0000	10.0000	12.0000	14.0000	16.0000	16.0000	18.0000];
d2int = [2.0000	2.0000	4.0000	6.0000	6.0000	8.0000	8.0000	10.0000	10.0000	12.0000	12.0000	14.0000	14.0000	16.0000	16.0000	18.0000	18.0000];


% pars = [k,a,b, c, vr, vt, vpeak, T, dt, f_start, f_end, f_time, Istore, C, d, k];
k = pars(1); a = pars(2); b = pars(3); c = pars(4); vr = pars(5); vpeak = pars(6);
T = pars(7); dt = pars(8); f_start=pars(9); f_end=pars(10); f_time=pars(11); Istore=pars(12);
% newly tuned parameters added to end of pars list % [C, vt, d]
C = pars(13); vt = pars(14); d = pars(15);

% dopamine stuff
DAtype = pars(16);
DA = pars(17); 

% multipliers
if DAtype == 1
    dD1 = d*(1-DA*x(2));
    vrD1 = vr*(1+DA*x(1));
else
    kD2 = k * (1-x(1)*DA);
end

% run simulation
t = 0:dt:T;         
n = length(t); % number of time points
nInj = numel(inj);

for loop = 1:nInj
    v = vr*ones(1,n); u=0*v;
    for i = 1:n-1
        switch DAtype
            case 1
                %--- D1
                v(i+1) = v(i) + dt*(k*(v(i)-vrD1)*(v(i)-vt)-u(i) + injScaled(loop))/C;
                u(i+1) = u(i) + dt*a*(b*(v(i)-vrD1)-u(i));
                % spikes?   
                if v(i+1)>=vpeak
                    v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+dD1;
                end
            case 2
                %--- D2
                v(i+1) = v(i) + dt*(kD2*(v(i)-vr)*(v(i)-vt)-u(i) + injScaled(loop))/C;
                u(i+1) = u(i) + dt*a*(b*(v(i)-vr)-u(i));
                % spikes?   
                if v(i+1)>=vpeak
                    v(i)=vpeak; v(i+1)=c; u(i+1)=u(i+1)+d;
                end
            otherwise
                error('Unknown receptor type')
        end
       
    end
     % time to first spike
    temp = find(v == vpeak); isis = diff(temp)*dt;
    if temp tfs(loop) = temp(1) * dt; else tfs(loop) = nan; end   % time in ms 
    fI(loop) = sum(v(f_start:f_end) == vpeak) ./ f_time;
    if isis fI1st(loop) = 1000./isis(1); else fI1st(loop) = 0; end
end

%-------------- compute fit
% %% SSE of f-I
% switch DAtype
%     case 1
%         fIerror = sum((fI-d1int).^2); % SSE
%     case 2
%         fIerror = sum((fI-d2int).^2); % SSE
% end

%% summed relative error f-I
switch DAtype
    case 1
        norm = d1int; norm(norm == 0) = 1;
        fIerror =  sum(abs(fI-d1int)./norm); % summed relative error
    case 2
        norm = d2int; norm(norm == 0) = 1;
        fIerror = sum(abs(fI-d2int)./norm); % summed relative error
end


f = fIerror; 

