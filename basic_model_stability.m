function [FP,JA,JB,Ev,VA,VB,class] = basic_model_stability(vr,vt,a,b,k,C,I)

% BASIC_MODEL_STABILITY linear stability analysis of Izhikevich neuron model
% [B,JA,JB,FP,Ev,VA,VB,CL] = BASIC_MODEL_STABILITY(vr,vt,a,b,k,C,I) given standard
% model parameters (vr,vt,b,k,C) and specified injection current I, computes
% the two fixed points FP (in (v,u) pairs per row), corresponding Jacobians
% JA and JB, and their corresponding eigenvalues Ev and eigenvector
% matrices VA, VB; will also classify both fixed points, returning types in
% cell array CL.
%
% Mark Humphries & Nathan Lepora 27/11/2008

% Report eq 29
B = sqrt((vr+vt+b/k)^2 - 4*(vr*vt+(b*vr+I)/k));

% Report eq 32
JA = [(b+k*B)/C -1/C; a*b -a];
JB = [(b-k*B)/C -1/C; a*b -a];

% fixed points - report eq 28
vA = 0.5*(vr+vt+b/k)+0.5*B; uA = b/2*(-vr+vt+b/k)+0.5*b*B;
vB = 0.5*(vr+vt+b/k)-0.5*B; uB = b/2*(-vr+vt+b/k)-0.5*b*B;
FP = [vA uA; vB uB];

% eigenvalues
[VA,DA] = eig(JA);
[VB,DB] = eig(JB);
eA = diag(DA)';
eB = diag(DB)';

Ev = [eA; eB];
    
% classify.... if fixed points exist
class{1} = 'No fixed points';
class{2} = 'No fixed points';
if isreal(vA) class{1} = classify(eA); end
if isreal(vB) class{2} = classify(eB); end


function class = classify(ev)
delta = ev(1)*ev(2); tau = ev(1) + ev(2);

if delta < 0
    class = 'saddle';      
elseif delta > 0
    div = tau^2 - 4*delta;
    if tau > 0
        if div > 0
            class = 'unstable node';
        elseif div < 0
            class = 'unstable spiral';
        else
            class = 'stars, degenerate nodes';
        end   
    elseif tau < 0
        if div > 0
            class = 'stable node';
        elseif div < 0
            class = 'stable spiral';
        else
            class = 'stars, degenerate nodes';
        end   
    else
        class = 'centre';
    end

else
    class = 'tiled fixed points';
end
