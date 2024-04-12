import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path

import ct_hiv_model_econ_eval.input_data as data
import ct_hiv_model_econ_eval.model_classes as model
import ct_hiv_model_econ_eval.param_classes as param
import ct_hiv_model_econ_eval.support as support

# selected therapy
therapy = param.Therapies.COMBO

# create a cohort
myCohort = model.Cohort(id=1,
                        pop_size=data.POP_SIZE,
                        parameters=param.Parameters(therapy=therapy))

# simulate the cohort over the specified time steps
myCohort.simulate(sim_length=data.SIM_LENGTH)

# plot the sample path (survival curve)
path.plot_sample_path(
    sample_path=myCohort.cohortOutcomes.nLivingPatients,
    title='Survival Curve',
    x_label='Time-Step (Year)',
    y_label='Number Survived')

# plot the histogram of survival times
hist.plot_histogram(
    data=myCohort.cohortOutcomes.survivalTimes,
    title='Histogram of Patient Survival Time',
    x_label='Survival Time (Year)',
    y_label='Count',
    bin_width=1)

# print the outcomes of this simulated cohort
support.print_outcomes(sim_outcomes=myCohort.cohortOutcomes,
                       therapy_name=therapy)
