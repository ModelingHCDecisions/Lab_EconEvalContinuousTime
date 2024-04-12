from enum import Enum

import deampy.markov as markov
import numpy as np

import ct_hiv_model_econ_eval.input_data as data
from ct_hiv_model_econ_eval.input_data import HealthStates


class Therapies(Enum):
    """ mono vs. combination therapy """
    MONO = 0
    COMBO = 1


class Parameters:
    def __init__(self, therapy):

        # selected therapy
        self.therapy = therapy

        # initial health state
        self.initialHealthState = HealthStates.CD4_200to500

        # annual treatment cost
        if self.therapy == Therapies.MONO:
            self.annualTreatmentCost = data.Zidovudine_COST
        else:
            self.annualTreatmentCost = data.Zidovudine_COST + data.Lamivudine_COST

        # calculate transition probabilities between hiv states
        prob_matrix_mono = get_trans_prob_matrix(trans_matrix=data.TRANS_MATRIX)

        # transition probability matrix of the selected therapy
        self.transRateMatrix = []

        if self.therapy == Therapies.MONO:
            # calculate transition rate matrix for the mono therapy
            self.transRateMatrix = get_trans_rate_matrix(trans_prob_matrix=prob_matrix_mono)

        elif self.therapy == Therapies.COMBO:
            # calculate transition probability matrix for the combination therapy
            self.transRateMatrix = get_trans_rate_matrix_combo(
                rate_matrix_mono=get_trans_rate_matrix(trans_prob_matrix=prob_matrix_mono),
                combo_rr=data.TREATMENT_RR)

        # annual state costs and utilities
        self.annualStateCosts = data.ANNUAL_STATE_COST
        self.annualStateUtilities = data.ANNUAL_STATE_UTILITY

        # discount rate
        self.discountRate = data.DISCOUNT


def get_trans_prob_matrix(trans_matrix):
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


def get_trans_rate_matrix(trans_prob_matrix):

    # find the transition rate matrix
    trans_rate_matrix = markov.discrete_to_continuous(
        trans_prob_matrix=trans_prob_matrix,
        delta_t=1)

    # calculate background mortality rate
    mortality_rate = -np.log(1 - data.ANNUAL_PROB_BACKGROUND_MORT)

    # add background mortality rate
    for row in trans_rate_matrix:
        row.append(mortality_rate)

    # add 2 rows for HIV death and natural death
    trans_rate_matrix.append([0] * len(HealthStates))
    trans_rate_matrix.append([0] * len(HealthStates))

    return trans_rate_matrix


def get_trans_rate_matrix_combo(rate_matrix_mono, combo_rr):
    """
    :param rate_matrix_mono: (list of lists) transition rate matrix under mono therapy
    :param combo_rr: relative risk of the combination treatment
    :returns (list of lists) transition rate matrix under combination therapy """

    # create an empty list of lists
    matrix_combo = []
    for row in rate_matrix_mono:
        matrix_combo.append([0]*len(row))  # adding a row [0, 0, 0, 0, 0]

    # populate the combo matrix
    # calculate the effect of combo-therapy on non-diagonal elements
    for s in range(len(matrix_combo)):
        # rates to HIV states
        for next_s in range(s + 1, len(HealthStates)-1):
            matrix_combo[s][next_s] = combo_rr * rate_matrix_mono[s][next_s]

        # rates of background mortality
        matrix_combo[s][-1] = rate_matrix_mono[s][-1]

    return matrix_combo


# tests
if __name__ == '__main__':
    probMatrix = get_trans_prob_matrix(data.TRANS_MATRIX)
    rateMatrixMono = get_trans_rate_matrix(probMatrix)
    rateMatrixCombo = get_trans_rate_matrix_combo(rateMatrixMono, data.TREATMENT_RR)

    print(rateMatrixMono)
    print(rateMatrixCombo)
