from quadratic_sieve_factorization import QuadraticSieveFactorizer
import argparse

"""
Factorization Driver Code
"""

parser = argparse.ArgumentParser()

parser.add_argument('--num', type=int, required=True, help='Number to factor')
parser.add_argument('--B', type=int, default=3000, required=False, help='B-smoothness value')
parser.add_argument('--M', type=int, default=50_000, required=False, help='High bound for factor candidate list to be sieved')
parser.add_argument('--qs_thresh', type=float, default=2.5, required=False, help='log-space score threshold for sieving')

msg = 'Bypass the pseudoprime test for bases: 2, 3, 5, 7. Only do this if you think you have a carmicheal number'
parser.add_argument('--bypass_psprime', type=bool, default=False, required=False, help=msg)

args = parser.parse_args()

factorizer = QuadraticSieveFactorizer(args.B, args.M, args.qs_thresh, args.bypass_psprime)

factor = factorizer(args.num)

if factor is not None:
    print("Factor found {factor}")
else:
    print("Factor not found")
