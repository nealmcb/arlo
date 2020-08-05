"""shim.py: Shim code to interface between the calling conventions expected by the current
bravo_sample_sizes() code with the API currently provided by the athena module.

Over time we expect both the Arlo and the Athena calling conventions to change, so
this is a very temporary solution.

TODO: Is is worth finding a way to keep the Audit objects cached?
Or is it better to make them up for each pairwise estimate as we go?
"""

import logging
import math
from athena.audit import Audit

def get_athena_test_statistics(
    risk_limit: float, p_w: float, p_r: float, sample_w: int, sample_r: int,
) -> (float, float):
    """
    Return Athena p-value and delta.
    TODO: refactor to pass in integer vote shares to allow more exact calculations, incorporate or
    track round schedule over time, and handle sampling without replacement.

    Inputs:
        risk_limit      - the risk-limit for this audit
        p_w             - the fraction of vote share for the winner
        p_r             - the fraction of vote share for the loser
        sample_w        - the number of votes for the winner that have already
                          been sampled
        sample_r        - the number of votes for the runner-up that have
                          already been sampled

    Outputs:
        p_value        - p-value for given circumstances
        delta:         Athena delta value (1/likelihood ratio)

    Test cases from https://github.com/gwexploratoryaudits/brla_explore/pull/10/files/988f068e65fd955c8e5d1512865ef5e95a1d7b3c..94693c67aa33a1c642a98336ca5b7fcd32c1ce33#
    test26: pass
    >>> get_athena_test_statistics(0.1, 0.224472184613, 0.12237580158, 50, 36)
    (0.08761871736, 0.7050552943)

    test27: fail
    >>> get_athena_test_statistics(0.1, 0.224472184613, 0.12237580158, 49, 37)
    (0.1245044151, 1.293272854)
    """

    # calculate the undiluted "two-way" share of votes for the winner
    p_wr = p_w + p_r
    p_w2 = p_w / p_wr

    ballots_cast = 100000
    a = int(ballots_cast * p_w2)
    b = ballots_cast - a
    if sample_w or sample_r:
        round_sizes = [sample_w + sample_r]
    else:
        round_sizes = []

    election = {
        "alpha": risk_limit,
        "delta": 1.0,
        "candidates": ["A", "B"],
        "results": [a, b],
        "total_ballots": ballots_cast,
        "winners": 1,
        "name": "pure_pair",
        "model": "bin",
    }

    a = Audit("athena", election["alpha"], election["delta"])
    a.add_election(election)
    a.add_round_schedule(round_sizes)

    if round_sizes:
        r = a.find_risk([sample_w])
    else:
        r = None

    logging.info(
        f"shim get_athena_test_statistics: margin {(p_w2 - 0.5) * 2} (pw {p_w} pr {p_r}) (sw {sample_w} sr {sample_w}) r {r}"
    )

    return r


def athena_sample_sizes(
    risk_limit: float,
    p_w: float,
    p_r: float,
    sample_w: int,
    sample_r: int,
    p_completion: float,
) -> int:
    """
    Return Athena round size based on completion probability, assuming the election outcome is correct.
    TODO: refactor to pass in integer vote shares to allow more exact calculations, incorporate or
    track round schedule over time, and handle sampling without replacement.

    Inputs:
        risk_limit      - the risk-limit for this audit
        p_w             - the fraction of vote share for the winner
        p_r             - the fraction of vote share for the loser
        sample_w        - the number of votes for the winner that have already
                          been sampled
        sample_r        - the number of votes for the runner-up that have
                          already been sampled
        p_completion    - the desired chance of completion in one round,
                          if the outcome is correct

    Outputs:
        sample_size     - the expected sample size for the given chance
                          of completion in one round

    >>> athena_sample_sizes(0.1, 0.6, 0.4, 56, 56, 0.7)
    244
    """

    # calculate the undiluted "two-way" share of votes for the winner
    p_wr = p_w + p_r
    p_w2 = p_w / p_wr

    ballots_cast = 100000
    a = int(ballots_cast * p_w2)
    b = ballots_cast - a
    pstop_goal = [p_completion]
    if sample_w or sample_r:
        round_sizes = [sample_w + sample_r]
    else:
        round_sizes = []

    election = {
        "alpha": risk_limit,
        "delta": 1.0,
        "candidates": ["A", "B"],
        "results": [a, b],
        "total_ballots": ballots_cast,
        "winners": 1,
        "name": "pure_pair",
        "model": "bin",
    }

    a = Audit("athena", election["alpha"], election["delta"])
    a.add_election(election)
    a.add_round_schedule(round_sizes)

    if round_sizes:
        r = a.find_risk([sample_w])
        below_kmin = max(r["required"]) - max(r["observed"])
    else:
        below_kmin = 0
    x = a.find_next_round_size(pstop_goal)
    next_round_size_0 = x["future_round_sizes"][0]

    next_round_size = next_round_size_0 + 2 * below_kmin

    size_adj = math.ceil(next_round_size / p_wr)

    logging.info(
        f"shim sample sizes: margin {(p_w2 - 0.5) * 2} (pw {p_w} pr {p_r}) (sw {sample_w} sr {sample_w}) pstop {p_completion} below_kmin {below_kmin} raw {next_round_size} scaled {size_adj}"
    )

    return size_adj


if __name__ == '__main__':
     import doctest
     doctest.testmod(verbose=True)
