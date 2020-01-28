"""
Athena election auditing calculations for Arlo.

Based on work by Poorvi Vora, Filip Zagorski, Neal McBurnett, Sarah Morin and Grant McClearn.

TODO: change name and API from aurror.py to athena.py here and in calling functions
"""

import math
from scipy.stats import binom
from audits.audit import RiskLimitingAudit


# For now, use the following two functions from https://github.com/filipzz/aurror/blob/master/code/aurror.py
# TODO: Replace with more general Athena module.


def next_round_prob(margin, round_size_prev, round_size, kmin, prob_table_prev):
    "Returns the probability distribution of being at a given risk level at the end of a given round."
    prob_table = [0] * (round_size + 1)
    for i in range(kmin + 1):
        for j in range(round_size + 1):
            prob_table[j] = prob_table[j] + binom.pmf(j-i, round_size - round_size_prev, (1+margin)/2) * prob_table_prev[i]

    return prob_table


def aurror(margin, alpha, round_schedule):
    round_schedule = [0] + round_schedule
    number_of_rounds = len(round_schedule)
    prob_table_prev = [1]
    prob_tied_table_prev = [1]
    kmins = [0] * number_of_rounds
    prob_sum = [0] * number_of_rounds
    prob_tied_sum = [0] * number_of_rounds

    for round in range(1,number_of_rounds):
        prob_table = next_round_prob(margin, round_schedule[round - 1], round_schedule[round], kmins[round - 1], prob_table_prev)
        prob_tied_table = next_round_prob(0, round_schedule[round - 1], round_schedule[round], kmins[round - 1], prob_tied_table_prev)

        kmin_found = False
        kmin_candidate = math.floor(round_schedule[round]/2)
        while kmin_found == False and kmin_candidate <= round_schedule[round]:
            if alpha * (sum(prob_table[kmin_candidate:len(prob_table)]) + prob_sum[round - 1]) >= (sum(prob_tied_table[kmin_candidate:len(prob_tied_table)]) + prob_tied_sum[round - 1]):
                kmin_found = True
                kmins[round] = kmin_candidate
                prob_sum[round] = sum(prob_table[kmin_candidate:len(prob_table)]) + prob_sum[round - 1]
                prob_tied_sum[round] = sum(prob_tied_table[kmin_candidate:len(prob_tied_table)]) + prob_tied_sum[round - 1]
            else:
                kmin_candidate = kmin_candidate + 1

        # cleaning prob_table/prob_tied_table
        for i in range(kmin_candidate, round_schedule[round] + 1):
            prob_table[i] = 0
            prob_tied_table[i] = 0

        prob_table_prev = prob_table
        prob_tied_table_prev = prob_tied_table

    return {"kmins" : kmins[1:len(kmins)], "prob_sum" : prob_sum[1:len(prob_sum)], "prob_tied_sum" : prob_tied_sum[1:len(prob_tied_sum)]}


# This class is currently based on functions from bravo.py, with new Athena / Aurror code being incorporated

class Aurror(RiskLimitingAudit):
    def __init__(self, risk_limit):
        super().__init__(risk_limit)


    # Based on BRAVO code for now.
    def get_expected_sample_sizes(self, margins, contests, sample_results):
        """
        Returns the expected sample size for a Aurror audit of each contest in contests.
        For now, estimate by calculating BRAVO's ASN and multiplying by AURROR_RATIO

        Input:
            margins - a dict of the margins for each contest passed from Sampler
            contests - a dict of the contests passed in from Sampler
            sample_results - a dict of the sample results from the Sampler

        Output:
            expected sample sizes - dict of computed expected sample size for each contest:
                {
                    contest1: asn1,
                    contest3: asn2,
                    ...
                }
        """

        # Temporary estimate of Aurror estimate as ration of BRAVO ASN
        AURROR_RATIO = 0.5

        asns = {}
        for contest in contests:
            margin = margins[contest]

            p_w = 10**7
            s_w = 0 
            p_l = 0
            # Get smallest p_w - p_l
            for winner in margin['winners']:
                if margin['winners'][winner]['p_w'] < p_w:
                    p_w = margin['winners'][winner]['p_w']

            if not margin['losers']:
                asns[contest] = -1
                continue

            for loser in margin['losers']:
                if margin['losers'][loser]['p_l'] > p_l:
                    p_l = margin['losers'][loser]['p_l']


            s_w = p_w/(p_w + p_l) 
            s_l = 1 - s_w

            if p_w == 1:
                # Handle single-candidate or crazy landslides
                asns[contest] = -1
            elif p_w == p_l:
                asns[contest] = contests[contest]['ballots']
            else: 
                z_w = math.log(2 * s_w)
                z_l = math.log(2 - 2 * s_w)
                p = p_w + s_l

                T = min(self.get_test_statistics(margins[contest], sample_results[contest]).values())

                weighted_alpha = math.log((1.0/self.risk_limit)/T)  
                asns[contest] = math.ceil((AURROR_RATIO * (weighted_alpha + (z_w / 2.0)) / (p_w*z_w + p_l*z_l)))

        return asns


    # Adapted from bravo.py
    def get_sample_sizes(self, contests, margins, sample_results):
        """
        Computes initial sample sizes parameterized by likelihood that the
        initial sample will confirm the election result, assuming no
        discrpancies.

        Inputs:
            sample_results - if a sample has already been drawn, this will
                             contain its results. 

        Outputs:
            samples - dictionary mapping confirmation likelihood to sample size, plus ASN info:
                    {
                       contest1:  {
                            'asn': {'size': asn, 'prob': stopping_probability},
                            likelihood1: sample_size,
                            likelihood2: sample_size,
                            ...
                        },
                        ...
                    }
        """

        quants = [.7, .8, .9]  # FIXME - still needed for error reports
        # For now, use multiples of expected size
        multiples = [1.2, 1.6, 2.1]

        samples = {}

        asns = self.get_expected_sample_sizes(margins, contests, sample_results)

        for contest in contests:
            samples[contest] = {}

            p_w = 10**7
            s_w = 0 
            p_l = 0
            best_loser = ''
            worse_winner = ''


            # For multi-winner, do nothing
            if 'numWinners' not in contests[contest] or contests[contest]['numWinners'] != 1:
                samples[contest] = {
                    'asn': {
                        'size': asns[contest],
                        'prob': None
                    }
                }
                return samples  # FIXME: doesn't complete loop over contests

            margin = margins[contest]
            # Get smallest p_w - p_l
            for winner in margin['winners']:
                if margin['winners'][winner]['p_w'] < p_w:
                    p_w = margin['winners'][winner]['p_w']
                    worse_winner = winner

            for loser in margin['losers']:
                if margin['losers'][loser]['p_l'] > p_l:
                    p_l = margin['losers'][loser]['p_l']
                    best_loser = loser

            # TODO: refactor out duplicated code to deal with ties and single-candidate race from bravo.
            # Put in audits.audit?

            # If we're in a single-candidate race, set sample to 0
            if not margin['losers']:
                samples[contest]['asn'] = {
                    'size': -1,
                    'prob': -1
                }
                for quant in quants:
                    samples[contest][quant] = -1

                continue

            s_w = p_w/(p_w + p_l) 
            s_l = 1 - s_w

            num_ballots = contests[contest]['ballots']

            # Handles ties
            if p_w == p_l:
                samples[contest]['asn'] = {
                    'size': num_ballots,
                    'prob': 1,
                }

                for quant in quants:
                    samples[contest][quant] = num_ballots
                continue

            # FIXME: need to look at all the sample margins, not just the closest reported pair
            sample_w = sample_results[contest][worse_winner]
            sample_l = sample_results[contest][best_loser]

            closest_margin = .1 # FIXME
            samples[contest]['asn'] = {
                'size': asns[contest],
                'prob': .52 # FIXME for Aurror.  self.expected_prob(p_w, p_l, sample_w, sample_l, asns[contest])
                }

            for multiple in multiples:
                sample_size = math.ceil(multiple * asns[contest])
                rounds = [sample_size]
                plan = aurror(closest_margin, self.risk_limit, rounds)
                prob = round(plan['prob_sum'][0], 2)
                samples[contest][prob] = sample_size
                # FIXME: feed in actual samples self.aurror_sample_sizes(p_w, p_l, sample_w, sample_l, quant)

        return samples


    def get_test_statistics(self, margins, sample_results):
        """
        Computes T*, the test statistic from an existing sample. 

        Inputs: 
            margins        - the margins for the contest being audited
            sample_results - mapping of candidates to votes in the (cumulative)
                             sample:

                    {
                        candidate1: sampled_votes,
                        candidate2: sampled_votes,
                        ...
                    }

        Outputs:
            T - Mapping of (winner, loser) pairs to their test statistic based
                on sample_results
        """

        winners = margins['winners']
        losers = margins['losers']

        T = {}

        # Setup pair-wise Ts:
        for winner in winners:
            for loser in losers:
                T[(winner, loser)] = 1

        # Handle the no-losers case
        if not losers:
            for winner in winners:
                T[(winner,)] = 1

        for cand, votes in sample_results.items():
            if cand in winners:
                for loser in losers:
                    T[(cand, loser)] *= (winners[cand]['swl'][loser]/0.5)**votes
            elif cand in losers:
                for winner in winners:
                    T[(winner, cand)] *= ((1 - winners[winner]['swl'][cand])/0.5)**votes

        return T


    def compute_risk(self, margins, sample_results):
        """
        Computes the risk-value of <sample_results> based on results in <contest>.

        Inputs: 
            margins        - the margins for the contest being audited
            sample_results - mapping of candidates to votes in the (cumulative)
                             sample:

                    {
                        candidate1: sampled_votes,
                        candidate2: sampled_votes,
                        ...
                    }

        Outputs:
            measurements    - the p-value of the hypotheses that the election
                              result is correct based on the sample, for each winner-loser pair. 
            confirmed       - a boolean indicating whether the audit can stop
        """

        T = self.get_test_statistics(margins, sample_results)

        measurements = {}
        finished = True
        for pair in T:
            measurements[pair] = 1/T[pair]
            if measurements[pair] > self.risk_limit:
                finished = False
        return measurements, finished 
