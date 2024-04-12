import ct_hiv_model_econ_eval.input_data as data
import ct_hiv_model_econ_eval.model_classes as model
import ct_hiv_model_econ_eval.param_classes as param
import ct_hiv_model_econ_eval.support as support

# simulating mono therapy
# create a cohort
cohort_mono = model.Cohort(id=0,
                           pop_size=data.POP_SIZE,
                           parameters=param.Parameters(therapy=param.Therapies.MONO))
# simulate the cohort
cohort_mono.simulate(sim_length=data.SIM_LENGTH)

# simulating combination therapy
# create a cohort
cohort_combo = model.Cohort(id=1,
                            pop_size=data.POP_SIZE,
                            parameters=param.Parameters(therapy=param.Therapies.COMBO))
# simulate the cohort
cohort_combo.simulate(sim_length=data.SIM_LENGTH)

# print the estimates for the mean survival time and mean time to AIDS
support.print_outcomes(sim_outcomes=cohort_mono.cohortOutcomes,
                       therapy_name=param.Therapies.MONO)
support.print_outcomes(sim_outcomes=cohort_combo.cohortOutcomes,
                       therapy_name=param.Therapies.COMBO)

# draw survival curves and histograms
support.plot_survival_curves_and_histograms(sim_outcomes_mono=cohort_mono.cohortOutcomes,
                                            sim_outcomes_combo=cohort_combo.cohortOutcomes)


# print comparative outcomes
support.print_comparative_outcomes(sim_outcomes_mono=cohort_mono.cohortOutcomes,
                                   sim_outcomes_combo=cohort_combo.cohortOutcomes)

# report the CEA results
support.report_CEA_CBA(sim_outcomes_mono=cohort_mono.cohortOutcomes,
                       sim_outcomes_combo=cohort_combo.cohortOutcomes)
