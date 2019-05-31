%%% script to look at MSN and FS bifurcation diagrams in I,v plane

load fit_model_results_NEWtuning

D1 = 0.8;
D2 = 0.8;
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
kD2 = k*(1-alpha*D2);

%%%% stability properties
% [FP,JA,JB,Ev,VA,VB,class] = basic_model_stability(vr,vt,a,b,k,C,I);
% 
% [FPd2,JAd2,JBd2,Evd2,VAd2,VBd2,classd2] = basic_model_stability(vr,vt,a,b,kD2,C,I);

%%%% get I,v curves
v = -90:1:-40;
I0_B = 0.25*k*(vr+vt+b/k)^2-k*vr*vt-b*vr;
I0_D1 = 0.25*k*(vrD1+vt+b/k)^2-k*vrD1*vt-b*vrD1;
I0_D2 = 0.25*kD2*(vr+vt+b/kD2)^2-kD2*vr*vt-b*vr;

IB = -k.*(v-vr).*(v-vt)+ b.*(v-vr);
ID1 = -k.*(v-vrD1).*(v-vt)+ b.*(v-vrD1);
ID2 = -kD2.*(v-vr).*(v-vt)+ b.*(v-vr);

% find bifurcation points quickly
IBix = find(diff(IB) <0,1); ID1ix = find(diff(ID1) <0,1); ID2ix = find(diff(ID2) <0,1); 

IB_stable = IB(1:IBix)'; IB_saddle = IB(IBix+1:end)'; vIBstable = v(1:IBix)'; vIBsaddle = v(IBix+1:end)';
ID1_stable = ID1(1:ID1ix)'; ID1_saddle = ID1(ID1ix+1:end)'; vID1stable = v(1:ID1ix)'; vID1saddle = v(ID1ix+1:end)';
ID2_stable = ID2(1:ID2ix)'; ID2_saddle = ID2(ID2ix+1:end)'; vID2stable = v(1:ID2ix)'; vID2saddle = v(ID2ix+1:end)';

figure(1); clf
title('Bifurcation diagrams for MSN model')
plot(IB_stable,vIBstable); hold on
plot(IB_saddle,vIBsaddle,'--')
plot(ID1_stable,vID1stable,'r');
plot(ID1_saddle,vID1saddle,'r--')
plot(ID2_stable,vID2stable,'k');
plot(ID2_saddle,vID2saddle,'k--')
xlabel('I (pA)'); ylabel('v (mV)')





