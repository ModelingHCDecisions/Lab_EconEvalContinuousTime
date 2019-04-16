import ParameterClasses as P
import SimPy.RandomVariantGenerators as RVGs
import SimPy.SamplePathClasses as Path
import SimPy.EconEvalClasses as Econ
import SimPy.StatisticalClasses as Stat
import SimPy.MarkovClasses as Markov


class Patient:
    def __init__(self, id, parameters):
        """ initiates a patient
        :param id: ID of the patient
        :param parameters: an instance of the parameters class
        """
        self.id = id
        self.rng = RVGs.RNG(seed=id)  # random number generator for this patient
        self.params = parameters
        # gillespie algorithm
        self.gillespie = Markov.Gillespie(transition_rate_matrix=parameters.rateMatrix)
        self.stateMonitor = PatientStateMonitor(parameters=parameters)  # patient state monitor

    def simulate(self, sim_length):
        """ simulate the patient over the specified simulation length """

        t = 0  # simulation time
        if_stop = False

        while not if_stop:

            # find time to next event, and next state
            dt, new_state_index = self.gillespie.get_next_state(
                current_state_index=self.stateMonitor.currentState.value,
                rng=self.rng)

            # stop if time to next event (dt) is None
            if dt is None:
                if_stop = True

            # else if  the next event occurs beyond simulation length
            elif dt + t > sim_length:
                if_stop = True
                # collect cost and health outcomes
                self.stateMonitor.costUtilityMonitor.update(time=sim_length,
                                                            current_state=self.stateMonitor.currentState,
                                                            next_state=self.stateMonitor.currentState)
            else:
                # advance time to the time of next event
                t += dt
                # update health state
                self.stateMonitor.update(time=t, new_state=P.HealthStates(new_state_index))


class PatientStateMonitor:
    """ to update patient outcomes (years survived, cost, etc.) throughout the simulation """
    def __init__(self, parameters):

        self.currentState = parameters.initialHealthState   # initial health state
        self.survivalTime = None      # survival time
        self.timeToAIDS = None        # time to develop AIDS
        self.ifDevelopedAIDS = False  # if the patient developed AIDS
        # patient's cost and utility monitor
        self.costUtilityMonitor = PatientCostUtilityMonitor(parameters=parameters)

    def update(self, time, new_state):
        """
        update the current health state to the new health state
        :param time: current time
        :param new_state: new state
        """

        # update survival time
        if new_state == P.HealthStates.HIV_DEATH or P.HealthStates.NATUAL_DEATH:
            self.survivalTime = time

        # update time until AIDS
        if self.currentState != P.HealthStates.AIDS and new_state == P.HealthStates.AIDS:
            self.ifDevelopedAIDS = True
            self.timeToAIDS = time

        # update cost and utility
        self.costUtilityMonitor.update(time=time,
                                       current_state=self.currentState,
                                       next_state=new_state)

        # update current health state
        self.currentState = new_state


class PatientCostUtilityMonitor:

    def __init__(self, parameters):

        self.tLastRecorded = 0  # time when the last cost and outcomes got recorded

        # model parameters for this patient
        self.params = parameters

        # total cost and utility
        self.totalDiscountedCost = 0
        self.totalDiscountedUtility = 0

    def update(self, time, current_state, next_state):
        """ updates the discounted total cost and health utility
        :param time: simulation time
        :param current_state: current health state
        :param next_state: next health state
        """

        # cost and utility (per unit of time) during the period since the last recording until now
        cost = self.params.annualStateCosts[current_state.value] + self.params.annualTreatmentCost
        utility = self.params.annualStateUtilities[current_state.value]

        # discounted cost and utility (continuously compounded)
        discounted_cost = Econ.pv_continuous_payment(payment=cost,
                                                     discount_rate=self.params.discountRate,
                                                     discount_period=(self.tLastRecorded, time))
        discounted_utility = Econ.pv_continuous_payment(payment=utility,
                                                        discount_rate=self.params.discountRate,
                                                        discount_period=(self.tLastRecorded, time))

        # update total discounted cost and utility
        self.totalDiscountedCost += discounted_cost
        self.totalDiscountedUtility += discounted_utility

        # update the time since last recording to the current time
        self.tLastRecorded = time


class Cohort:
    def __init__(self, id, pop_size, parameters):
        """ create a cohort of patients
        :param id: cohort ID
        :param pop_size: population size of this cohort
        :param parameters: parameters
        """
        self.id = id
        self.popSize = pop_size
        self.params = parameters
        self.patients = []  # list of patients
        self.cohortOutcomes = CohortOutcomes()  # outcomes of the this simulated cohort

    def simulate(self, sim_length):
        """ simulate the cohort of patients over the specified number of time-steps
        :param sim_length: simulation length
        """

        # populate the cohort
        for i in range(self.popSize):
            # create a new patient (use id * pop_size + n as patient id)
            patient = Patient(id=self.id * self.popSize + i, parameters=self.params)
            # add the patient to the cohort
            self.patients.append(patient)

        # simulate all patients
        for patient in self.patients:
            # simulate
            patient.simulate(sim_length)

        # store outputs of this simulation
        self.cohortOutcomes.extract_outcomes(self.patients)

        # clear patients
        self.patients.clear()


class CohortOutcomes:
    def __init__(self):

        self.survivalTimes = []         # patients' survival times
        self.timesToAIDS = []           # patients' times to AIDS
        self.costs = []                 # patients' discounted costs
        self.utilities =[]              # patients' discounted utilities
        self.nLivingPatients = None     # survival curve (sample path of number of alive patients over time)

        self.statSurvivalTime = None    # summary statistics for survival time
        self.statTimeToAIDS = None      # summary statistics for time to AIDS
        self.statCost = None            # summary statistics for discounted cost
        self.statUtility = None         # summary statistics for discounted utility

    def extract_outcomes(self, simulated_patients):
        """ extracts outcomes of a simulated cohort
        :param simulated_patients: a list of simulated patients"""

        # record survival time and time until AIDS
        for patient in simulated_patients:
            # survival time
            if not (patient.stateMonitor.survivalTime is None):
                self.survivalTimes.append(patient.stateMonitor.survivalTime)
            # time until AIDS
            if patient.stateMonitor.ifDevelopedAIDS:
                self.timesToAIDS.append(patient.stateMonitor.timeToAIDS)
            # discounted cost and discounted utility
            self.costs.append(patient.stateMonitor.costUtilityMonitor.totalDiscountedCost)
            self.utilities.append(patient.stateMonitor.costUtilityMonitor.totalDiscountedUtility)

        # summary statistics
        self.statSurvivalTime = Stat.SummaryStat('Survival time', self.survivalTimes)
        self.statTimeToAIDS = Stat.SummaryStat('Time until AIDS', self.timesToAIDS)
        self.statCost = Stat.SummaryStat('Discounted cost', self.costs)
        self.statUtility = Stat.SummaryStat('Discounted utility', self.utilities)

        # survival curve
        self.nLivingPatients = Path.PrevalencePathBatchUpdate(
            name='# of living patients',
            initial_size=len(simulated_patients),
            times_of_changes=self.survivalTimes,
            increments=[-1]*len(self.survivalTimes)
        )
