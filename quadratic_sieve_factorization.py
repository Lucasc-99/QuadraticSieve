"""
Quadratic Sieve Algorithm - Written by Lucas Cecchi
"""
import math


def jacobi(n, k):
    if k == 2:
        return 1
    assert (k > 0 and k % 2 == 1)
    n = n % k
    t = 1
    while n != 0:
        while n % 2 == 0:
            n /= 2
            r = k % 8
            if r == 3 or r == 5:
                t = -t
        n, k = k, n
        if n % 4 == k % 4 == 3:
            t = -t
        n %= k
    if k == 1:
        return t
    else:
        return 0


def euclid(a, b):
    while b:
        a, b = b, a % b
    return a


def extended_euclid(a, b):
    """
    :param a: an integer
    :param b: an integer
    :return:
    old_s: n
    old_t: m
    old_r: gcd(a,b)
    memo: r_i computation memo
    """
    memo = []

    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1

    while r != 0:
        q = old_r // r

        prov = r
        r = old_r - q * prov
        old_r = prov

        prov = s
        s = old_s - q * prov
        old_s = prov

        prov = t
        t = old_t - q * prov
        old_t = prov

        memo.append(old_r)

    return old_s, old_t, old_r, memo


"""
Trial Division, Tonelli's, Strong Pseudo-Prime Test
"""


def trial_division(primes, n):
    out = []
    for p in primes:
        exp = 0
        while p ** exp <= n and n % p ** exp == 0:
            exp += 1
        exp -= 1 if exp > 0 else 0
        out.append(exp)
        n /= p ** exp

    return out, int(n)


def strong_pseudo_prime(bases, n):
    if n == 2:
        return True
    elif n < 2 or n % 2 == 0:
        return False

    a, t = trial_division([2], n - 1)

    def base_b_pseudo(base):
        if base == n:
            return True

        t_test = t
        b_test = pow(base, t_test, n)

        if b_test == n - 1 or b_test == 1:
            return True

        for _ in range(a[0]):
            b_test = pow(b_test, 2, n)
            if b_test == n - 1:
                return True

        return False

    is_pseudo_prime = True
    for b in bases:
        is_pseudo_prime = is_pseudo_prime and base_b_pseudo(b)

    return is_pseudo_prime


def inverse_modulo(a, n):
    s, t, r, memo = extended_euclid(a, n)
    return (s % n + n) % n if r == 1 else None


def find_non_residue(p, bases=(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31,
                               37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)):
    for b in bases:
        if jacobi(b, p) == -1:
            return b


def tonelli(a, p):
    if jacobi(a, p) != 1:
        return None
    s_l, t = trial_division([2], p - 1)
    s = s_l[0]
    b = find_non_residue(p)
    i = 2
    i_l = [2]

    for k in range(2, s + 1):
        if pow(a * inverse_modulo(b, p) ** i, t * (2 ** (s - k)), p) != 1:
            i = i + (2 ** (k - 1))
        i_l.append(i)
    return ((b ** (i // 2)) * ((a * inverse_modulo(b, p) ** i) ** ((t + 1) // 2))) % p, i_l


"""
Linear Algebra Tools
"""


def reduce_mat(M, mod=2):
    for row_i in range(len(M)):
        for col_i in range(len(M[row_i])):
            M[row_i][col_i] = int(M[row_i][col_i] % mod)


def rref_mod2(M):
    reduce_mat(M)

    lead = 0
    rowCount = len(M)
    columnCount = len(M[0])
    for r in range(rowCount):
        if lead >= columnCount:
            return
        i = r
        while M[i][lead] == 0:
            i += 1
            if i == rowCount:
                i = r
                lead += 1
                if columnCount == lead:
                    return

        # Swap rows i and r
        if i != r:
            M[i], M[r] = M[r], M[i]

        for i in range(rowCount):
            if i != r:
                # Subtract M[i, lead] multiplied by row r from row i
                lead_val = M[i][lead]
                M[i] = [int((m_i - (lead_val * m_r)) % 2)
                        for m_r, m_i in zip(M[r], M[i])
                        ]
        lead += 1


def get_solution_basis(V):
    # Reduce
    rref_mod2(V)

    # Find free variables
    pivot_col = set()
    for i in range(len(V)):
        for j in range(len(V[0])):
            if V[i][j] == 1:
                pivot_col.add(j)
                break
    free_list = set([i for i in range(len(V[0]))]) - pivot_col
    sol_vecs = []

    for free_var in free_list:
        sol_vec = []
        sol_idx = 0
        row_idx = 0
        while sol_idx < len(V[0]):
            if sol_idx in free_list:
                sol_vec.append(int(sol_idx == free_var))
            else:
                sol_vec.append(V[row_idx][free_var])
                row_idx += 1
            sol_idx += 1

        sol_vecs.append(sol_vec)

    return sol_vecs


"""
Quadratic Sieve tools
"""


def generate_factor_base(B, n):
    # Note: this may allow composites for very large B
    for b in reversed(range(3, B)):
        if strong_pseudo_prime([2, 3, 5, 7, 11], b) and jacobi(n, b) == 1:
            yield b
    yield 2 if B > 2 else None



def quadratic_sieve(n, B, M, qs_threshold, bypass_psuedoprime_test=False):

    if not bypass_psuedoprime_test and strong_pseudo_prime([2, 3, 5, 7, 11], n):
        print(f"{n} is likely prime")
        return

    # Preliminaries
    k = math.ceil(math.sqrt(n))
    r_l = [r for r in range(k, k + M + 1)]
    s_r = [math.log((r ** 2) - n) for r in r_l]
    factor_base = list(reversed([b for b in generate_factor_base(B, n)]))
    residue_base = [1] + [tonelli(n, p)[0] for p in factor_base[1:]]

    sub_count = 0
    # Sieve:
    # subtract log(p) from divisible r_i
    for t, p in zip(residue_base, factor_base):
        r_i = 0
        seen = set()
        while r_i < len(r_l) and r_l[r_i] % p != t:
            r_i += 1
        while r_i < len(r_l) and r_l[r_i] % p == t:
            seen.add(r_l[r_i])
            sub_count += 1
            s_r[r_i] -= math.log(p)
            r_i += p
        r_i = 0
        while r_i < len(r_l) and r_l[r_i] % p != (p - t):
            r_i += 1
        while r_i < len(r_l) and r_l[r_i] not in seen and r_l[r_i] % p == (p - t):
            sub_count += 1
            s_r[r_i] -= math.log(p)
            r_i += p

    sieved_r = []
    for s_i in range(len(s_r)):
        if abs(s_r[s_i]) < qs_threshold:
            sieved_r.append(r_l[s_i])

    if len(sieved_r) < len(factor_base) + 1:
        print("Not enough r values, sieve may fail")

    # Get factorizations, reverse and transpose
    M_tup = [trial_division(factor_base, r ** 2 - n) for r in sieved_r]

    # Make sure all the values are B-smooth
    for row in M_tup:
        if row[1] != 1:
            print("Sieved r values are not all B-smooth, adjust parameters or sieve may fail")
            return

    M = [m[0] for m in M_tup]

    M_t = [[M[j][i] for j in range(len(M))] for i in range(len(M[0]))]

    # Solve
    solutions = get_solution_basis(M_t)

    for sol in solutions:
        # Calculate x
        x = 1
        for i in range(len(sol)):
            if sol[i]:
                x *= sieved_r[i]
                x %= n

        y_exp = [0] * len(factor_base)
        for i in range(len(sol)):
            if sol[i] == 1:
                for exponent_i in range(len(M[i])):
                    y_exp[exponent_i] += M[i][exponent_i]
        y = 1
        for i in range(len(y_exp)):
            y *= factor_base[i] ** (y_exp[i] // 2)
            y %= n

        assert (y * y) % n == (x * x) % n

        # Kraitchik's
        f = euclid(abs(x - y), n)
        if f != 1 and f != n:
            assert n % f == 0
            return f


class QuadraticSieveFactorizer:

    def __init__(self, B, M, qs_threshold, bypass_pseudo_prime_test) -> None:
        self.B = B
        self.M = M
        self.qs_threshold = qs_threshold
        self.bypass_pseudo_prime_test = bypass_pseudo_prime_test

    def __call__(self, n):
        return quadratic_sieve(n, self.B, self.M, self.qs_threshold, self.bypass_pseudo_prime_test)
        
    
        