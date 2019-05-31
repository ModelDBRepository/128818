This ZIP contains all the Matlab code necessary for tuning and
studying the dopamine-modulated MSN models (Humphries et
al. 2009). The code should all run straight out of the box. It is a
little cleaner than the work-in-progress code; we've removed dead code
to avoid confusion over the various dead-ends and checks etc.

For questions and assistance contact: m.d.humphries@sheffield.ac.uk or
drmdhumphries@gmail.com


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

The key functions are:
handtune_Iz_MSneuron.m - testbed for the construction of the D1 and D2
model extensions, and handtuning of their parameters

fit_Moyer_model.m - fine-tuning of the fits to the f-I and f-f curves
for the base, D1, and D2 MSN models. This function calls further
functions (stageone.m, stagetwo.m, stagethree.m and stagefour.m) that
correspond to the flow-chart in Figure 1 of the paper.  It also calls
the validation fit function stagefive.m.

test_Moyer_model_fit.m - takes the results of the fine-tuning, and
assesses them against all the available data (Figure 2, Figure 3C)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Other functions:

effects_of_D1int_pars.m - effect on D1 intrinsic model f-I curves of
changing D1 level (Figure 3A)

effects_of_D2int_pars.m - effect on D2 intrinsic model f-I curves of
changing D2 level (Figure 3B)

effects_of_D1synaptic_pars.m - effect on D1 complete model f-f curves
				of changing D1 level (Figure 3D); also
				runs on D1 intrinsic model to show
				that "signal-to-noise" effect is still
				there if there is no action of D1 on
				NMDA input

effects_of_D2synaptic_pars.m - effect on D2 complete model f-f curves
of changing D2 level (Figure 3E)

TTFS_spike_input.m - time-to-first-spike estimates under synaptic
input for baseline/D1/D2 (Figure 3F)

paired_pulse_response.m - show the paired-pulse facilitation, its
dependence on inter-pulse interval, and effects of the dopamine models
(Figure 4)

bifurcation_diagrams.m - the bifurcation curves from Figure 5

effects_of_NMDA_on_single_trial_bimodality.m - run NMDA agonist
						simulations for
						different NMDA
						conductance levels,
						and assess
						distribution of
						membrane
						potential. (Figure 6)

Detailed_effects_of_NMDA_on_single_trial_bimodality.m - gathers the
data and runs tests for Figure 7

effects_of_NMDA_nogate_on_single_trial_bimodality.m - run NMDA agonist
simulations without voltage gate (Figure 6)

effects_of_AMPA_on_single_trial_bimodality.m - run NMDA agonist
simulations using AMPA multiplier instead (Figure 6)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Helper functions:
spkgen.m - generates the spike-event input for various other functions
basic_model_stability.m - computes the Jacobian and eigenvalues of a
given basic Izhikevich model, and determines the number and type of
fixed points.  fitcurves.m - fits selected family of functions to
passed data. Uses lsqcurvefit.m, hence requires Optimization Toolbox.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

We include the tuning results file that we obtained, from which much
of the paper was derived: fit_model_NEWtuning.mat (found MSN
parameters, from fit_Moyer_model.m)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NOTES:

(1) Given different MATLAB versions and platforms, it is likely that
running fit_Moyer_model.m will results in slightly different values to
the ones we found.
(2) we have also tested a version that does NOT include the 1/tau_s
scaling of the post-synaptic potential model; this fits the Moyer et
al data well (and also has all the same attributes as the model
described in the paper). The principle difference is that the synaptic
conductances are reduced by an order of magnitude

Humphries, M. D., Lepora, N., Wood, R. & Gurney, K. (2009) Capturing
dopaminergic modulation and bimodal membrane behaviour of striatal
medium spiny neurons in accurate, reduced models. Frontiers in
Computational Neuroscience, 3, 26
