from enum import Enum
import numpy as np
import InputData as Data
import SimPy.MarkovClasses as Markov


class HealthStates(Enum):
    """ health states of patients with HIV """
    CD4_200to500 = 0
    CD4_200 = 1
    AIDS = 2
    HIV_DEATH = 3
    NATUAL_DEATH = 4


class Therapies(Enum):
    """ mono vs. combination therapy """
    MONO = 0
    COMBO = 1


class ParametersFixed:
    def __init__(self, therapy):

        # selected therapy
        self.therapy = therapy

        # initial health state
        self.initialHealthState = HealthStates.CD4_200to500

        # annual treatment cost
        if self.therapy == Therapies.MONO:
            self.annualTreatmentCost = Data.Zidovudine_COST
        else:
            self.annualTreatmentCost = Data.Zidovudine_COST + Data.Lamivudine_COST

        # transition probability matrix of the selected therapy
        self.rateMatrix = []

        # calculate transition probabilities between hiv states

        prob_matrix_mono = get_prob_matrix(trans_matrix=Data.TRANS_MATRIX)

        if self.therapy == Therapies.MONO:
            # calculate transition probability matrix for the mono therapy
            self.rateMatrix = get_rate_matrix_mono(trans_prob_matrix=prob_matrix_mono)

        elif self.therapy == Therapies.COMBO:
            # calculate transition probability matrix for the combination therapy
            self.rateMatrix = get_rate_matrix_combo(
                rate_matrix_mono=get_rate_matrix_mono(trans_prob_matrix=prob_matrix_mono),
                combo_rr=Data.TREATMENT_RR)

        # annual state costs and utilities
        self.annualStateCosts = Data.ANNUAL_STATE_COST
        self.annualStateUtilities = Data.ANNUAL_STATE_UTILITY

        # discount rate
        self.discountRate = Data.DISCOUNT


def get_prob_matrix(trans_matrix):
    """
    :param trans_matrix: transition matrix containing counts of transitions between states
    :return: transition probability matrix
    """

    # initialize transition probability matrix
    trans_prob_matrix = []

    # for each row in the transition matrix
    for row in trans_matrix:
        # calculate the transition probabilities
        prob_row = np.array(row)/sum(row)
        # add this row of transition probabilities to the transition probability matrix
        trans_prob_matrix.append(prob_row)

    return trans_prob_matrix


def get_rate_matrix_mono(trans_prob_matrix):

    # find the transition rate matrix
    trans_rate_matrix = Markov.discrete_to_continuous(
        prob_matrix=trans_prob_matrix,
        delta_t=1)

    # calculate background mortality rate
    mortality_rate = -np.log(1 - Data.ANNUAL_PROB_BACKGROUND_MORT)

    # add background mortality rate
    for row in trans_rate_matrix:
        row.append(mortality_rate)

    # add 2 rows for HIV death and natural death
    trans_rate_matrix.append([0] * len(HealthStates))
    trans_rate_matrix.append([0] * len(HealthStates))

    return trans_rate_matrix


def get_rate_matrix_combo(rate_matrix_mono, combo_rr):
    """
    :param rate_matrix_mono: (list of lists) transition rate matrix under mono therapy
    :param combo_rr: relative risk of the combination treatment
    :returns (list of lists) transition rate matrix under combination therapy """

    # create an empty list of lists
    matrix_combo = []
    for row in rate_matrix_mono:
        matrix_combo.append([0]*len(row))  # adding a row [0, 0, 0]

    # populate the combo matrix
    # calculate the effect of combo-therapy on non-diagonal elements
    for s in range(len(matrix_combo)):
        for next_s in range(s + 1, len(HealthStates)):
            matrix_combo[s][next_s] = combo_rr * rate_matrix_mono[s][next_s]

    return matrix_combo


# # tests
# matrix_mono = get_rate_matrix(Data.TRANS_MATRIX)
# matrix_combo = get_rate_matrix_combo(matrix_mono, Data.TREATMENT_RR)
#
# print(matrix_mono)
# print(matrix_combo)
