import deampy.econ_eval as econ
import deampy.plots.histogram as hist
import deampy.plots.sample_paths as path
import deampy.statistics as stats

import EconEvalInputData as D


def print_outcomes(sim_outcomes, therapy_name):
    """ prints the outcomes of a simulated cohort
    :param sim_outcomes: outcomes of a simulated cohort
    :param therapy_name: the name of the selected therapy
    """
    # mean and confidence interval of patient survival time
    survival_mean_CI_text = sim_outcomes.statSurvivalTime\
        .get_formatted_mean_and_interval(interval_type='c',
                                         alpha=D.ALPHA,
                                         deci=2)

    # mean and confidence interval text of time to AIDS
    time_to_HIV_death_CI_text = sim_outcomes.statTimeToAIDS\
        .get_formatted_mean_and_interval(interval_type='c',
                                         alpha=D.ALPHA,
                                         deci=2)

    # mean and confidence interval text of discounted total cost
    cost_mean_CI_text = sim_outcomes.statCost\
        .get_formatted_mean_and_interval(interval_type='c',
                                         alpha=D.ALPHA,
                                         deci=0,
                                         form=',')

    # mean and confidence interval text of discounted total utility
    utility_mean_CI_text = sim_outcomes.statUtility\
        .get_formatted_mean_and_interval(interval_type='c',
                                         alpha=D.ALPHA,
                                         deci=2)

    # print outcomes
    print(therapy_name)
    print("  Estimate of mean survival time and {:.{prec}%} confidence interval:".format(1 - D.ALPHA, prec=0),
          survival_mean_CI_text)
    print("  Estimate of mean time to AIDS and {:.{prec}%} confidence interval:".format(1 - D.ALPHA, prec=0),
          time_to_HIV_death_CI_text)
    print("  Estimate of discounted cost and {:.{prec}%} confidence interval:".format(1 - D.ALPHA, prec=0),
          cost_mean_CI_text)
    print("  Estimate of discounted utility and {:.{prec}%} confidence interval:".format(1 - D.ALPHA, prec=0),
          utility_mean_CI_text)
    print("")


def plot_survival_curves_and_histograms(sim_outcomes_mono, sim_outcomes_combo):
    """ draws the survival curves and the histograms of time until HIV deaths
    :param sim_outcomes_mono: outcomes of a cohort simulated under mono therapy
    :param sim_outcomes_combo: outcomes of a cohort simulated under combination therapy
    """

    # get survival curves of both treatments
    survival_curves = [
        sim_outcomes_mono.nLivingPatients,
        sim_outcomes_combo.nLivingPatients
    ]

    # graph survival curve
    path.plot_sample_paths(
        sample_paths=survival_curves,
        title='Survival curve',
        x_label='Simulation time step (year)',
        y_label='Number of alive patients',
        legends=['Mono Therapy', 'Combination Therapy'],
        color_codes=['green', 'blue']
    )

    # histograms of survival times
    set_of_survival_times = [
        sim_outcomes_mono.survivalTimes,
        sim_outcomes_combo.survivalTimes
    ]

    # graph histograms
    hist.plot_histograms(
        data_sets=set_of_survival_times,
        title='Histogram of patient survival time',
        x_label='Survival time (year)',
        y_label='Counts',
        bin_width=1,
        legends=['Mono Therapy', 'Combination Therapy'],
        color_codes=['green', 'blue'],
        transparency=0.6
    )


def print_comparative_outcomes(sim_outcomes_mono, sim_outcomes_combo):
    """ prints average increase in survival time, discounted cost, and discounted utility
    under combination therapy compared to mono therapy
    :param sim_outcomes_mono: outcomes of a cohort simulated under mono therapy
    :param sim_outcomes_combo: outcomes of a cohort simulated under combination therapy
    """

    # increase in mean survival time under combination therapy with respect to mono therapy
    increase_survival_time = stats.DifferenceStatIndp(
        name='Increase in mean survival time',
        x=sim_outcomes_combo.survivalTimes,
        y_ref=sim_outcomes_mono.survivalTimes)

    # estimate and CI
    estimate_CI = increase_survival_time.get_formatted_mean_and_interval(interval_type='c',
                                                                         alpha=D.ALPHA,
                                                                         deci=2)
    print("Increase in mean survival time and {:.{prec}%} confidence interval:"
          .format(1 - D.ALPHA, prec=0),
          estimate_CI)

    # increase in mean discounted cost under combination therapy with respect to mono therapy
    increase_discounted_cost = stats.DifferenceStatIndp(
        name='Increase in mean discounted cost',
        x=sim_outcomes_combo.costs,
        y_ref=sim_outcomes_mono.costs)

    # estimate and CI
    estimate_CI = increase_discounted_cost.get_formatted_mean_and_interval(interval_type='c',
                                                                           alpha=D.ALPHA,
                                                                           deci=2,
                                                                           form=',')
    print("Increase in mean discounted cost and {:.{prec}%} confidence interval:"
          .format(1 - D.ALPHA, prec=0),
          estimate_CI)

    # increase in mean discounted utility under combination therapy with respect to mono therapy
    increase_discounted_utility = stats.DifferenceStatIndp(
        name='Increase in mean discounted utility',
        x=sim_outcomes_combo.utilities,
        y_ref=sim_outcomes_mono.utilities)

    # estimate and CI
    estimate_CI = increase_discounted_utility.get_formatted_mean_and_interval(interval_type='c',
                                                                              alpha=D.ALPHA,
                                                                              deci=2)
    print("Increase in mean discounted utility and {:.{prec}%} confidence interval:"
          .format(1 - D.ALPHA, prec=0),
          estimate_CI)


def report_CEA_CBA(sim_outcomes_mono, sim_outcomes_combo):
    """ performs cost-effectiveness and cost-benefit analyses
    :param sim_outcomes_mono: outcomes of a cohort simulated under mono therapy
    :param sim_outcomes_combo: outcomes of a cohort simulated under combination therapy
    """

    # define two strategies
    mono_therapy_strategy = econ.Strategy(
        name='Mono Therapy',
        cost_obs=sim_outcomes_mono.costs,
        effect_obs=sim_outcomes_mono.utilities,
        color='green'
    )
    combo_therapy_strategy = econ.Strategy(
        name='Combination Therapy',
        cost_obs=sim_outcomes_combo.costs,
        effect_obs=sim_outcomes_combo.utilities,
        color='blue'
    )

    # do CEA
    # (the first strategy in the list of strategies is assumed to be the 'Base' strategy)
    CEA = econ.CEA(
        strategies=[mono_therapy_strategy, combo_therapy_strategy],
        if_paired=False
    )

    # plot cost-effectiveness figure
    CEA.plot_CE_plane(
        title='Cost-Effectiveness Analysis',
        x_label='Additional QALYs',
        y_label='Additional Cost',
        x_range=(-1, 5),
        y_range=(-10000, 90000),
        interval_type='c'  # to show confidence intervals for cost and effect of each strategy
    )

    # report the CE table
    CEA.build_CE_table(
        interval_type='c',
        alpha=D.ALPHA,
        cost_digits=0,
        effect_digits=2,
        icer_digits=2,
        file_name='CETable.csv')

    # CBA
    NBA = econ.CBA(
        strategies=[mono_therapy_strategy, combo_therapy_strategy],
        wtp_range=[0, 50000],
        if_paired=False
    )
    # show the net monetary benefit figure
    NBA.plot_marginal_nmb_lines(
        title='Cost-Benefit Analysis',
        x_label='Willingness-to-pay per QALY ($)',
        y_label='Marginal Net Monetary Benefit ($)',
        interval_type='c',
        show_legend=True,
        figure_size=(6, 5)
    )