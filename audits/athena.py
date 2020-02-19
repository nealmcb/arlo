"""
Athena election auditing calculations for Arlo.

Based on work by Poorvi Vora, Filip Zagorski, Neal McBurnett, Sarah Morin and Grant McClearn.
"""

import math
from scipy.stats import binom
from audits.audit import RiskLimitingAudit

from audits.r2b2 import athena

# This class is currently based on functions from bravo.py, with new Athena code being incorporated

class Athena(RiskLimitingAudit):
    def __init__(self, risk_limit):
        super().__init__(risk_limit)


    # Based on BRAVO code for now.
    def get_expected_sample_sizes(self, margins, contests, sample_results):
        """
        Returns the expected sample size for a Athena audit of each contest in contests.
        For now, estimate by calculating BRAVO's ASN and multiplying by ATHENA_RATIO

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

        # Temporary estimate of Athena estimate as ration of BRAVO ASN
        ATHENA_RATIO = 0.5

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
                asns[contest] = math.ceil((ATHENA_RATIO * (weighted_alpha + (z_w / 2.0)) / (p_w*z_w + p_l*z_l)))

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

        quants = [.7, .8, .9]

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

            closest_margin = p_w - p_l
            samples[contest]['asn'] = {
                'size': asns[contest],
                'prob': .52 # FIXME for Athena.  self.expected_prob(p_w, p_l, sample_w, sample_l, asns[contest])
                }

            athena_audit = athena.AthenaAudit()

            for quant in quants:
                #sample_size = math.ceil(multiple * asns[contest])
                #rounds = [sample_size]
                ## plan = athena.athena(closest_margin, self.risk_limit, rounds)
                # prob = round(plan['prob_sum'][0], 2)
                samples[contest][quant] = athena_audit.find_next_round_size(closest_margin, self.risk_limit, [], quant, 10)['size']
                print(f"quant: {samples[contest][quant]}")
                # import pdb; pdb.set_trace()
                # FIXME: where does 10 come from (besides current find_next_round_sizes()...)
                # samples[contest][prob] = sample_size
                # FIXME: feed in actual samples self.athena_sample_sizes(p_w, p_l, sample_w, sample_l, quant)

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
