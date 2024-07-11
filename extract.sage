from bisect import bisect

# too messy to consider p = 2 so I just ignored it
def root_extract(a, r, arr_p_k):
    res = []
    for (p, k) in arr_p_k:
        d = gcd(r, p - 1)
        # d solutions in total, and x^d ≡ x^(d*(y*r/d+z*(p-1)/d)) ≡ x^(y*r) ≡ a^y (mod p), y ≡ (r/d)^(-1) (mod (p-1)/d)
        a_ = power_mod(a, inverse_mod(r // d, (p - 1) // d), p)
        if d != 1:
            sol = AMM_recurr(a_, d, p)
        else:
            sol = [a_]
        # even for d = 1, it's much faster than directly calculate x ≡ a^y (mod p^k), y ≡ r^(-1) (mod p^k*(p-1)) when p or k is large
        if k != 1:
            sol = hensel_lift(sol, r, a, p, k)
        res.append(sol)
    return res

# can be greatly improved by parting k into a sum of l numbers such that all the numbers divides t and minimizing l, but now I'm lazy
def AMM_recurr(a, d, p):
    # r is not large, a piece of cake
    factors = d.factor()
    sol = [a]
    q = randint(2, p)

    prod = 1
    for (x, y) in factors:
        prod *= x
    if prod == 2:
        # Jacobi symbol is much faster than Euler's criterion
        while kronecker(q, p) == 1:
            q = randint(2, p)
    else:
        power = (p - 1) // prod
        while power_mod(q, power, p) == 1:
            q = randint(2, p)

    for (d0, k) in factors:
        s = p - 1
        t = 0
        while s % d0 == 0:
            s //= d0
            t += 1
        b = inverse_mod(d0, s)
        qs = power_mod(q, s, p)
        for i in range(k):
            next_sol = []
            for c in sol:
                next_sol += AMM(c, b, d0, p, qs, t)
            sol = next_sol
    return sol

def AMM(a, b, d, p, qs, t):
    # a^((d*b-1)*d^(t-1)) ≡ 1 (mod p)
    # c ≡ a^(d*b-1) (mod p)
    c = power_mod(a, d * b - 1, p)
    sol = power_mod(a, b, p)

    # q^(s*d^(t-1))
    qsd = power_mod(qs, d ^ (t - 1), p)
    # q^(s*i*d^(t-1))
    qsd_tab = [(1, 0)]
    for i in range(1, d):
        qsd_tab.append(((qsd_tab[i - 1][0] * qsd) % p, i))
    qsd_tab.sort()
    qsds = [x[0] for x in qsd_tab]
    for i in range(t - 1):
        tmp = power_mod(c, d ^ (t - 2 - i), p)
        next_qs = power_mod(qs, d, p)
        if tmp != 1:
            index = d - qsd_tab[bisect(qsds, tmp) - 1][1]
            c = (c * power_mod(next_qs, index, p)) % p
            sol = (sol * power_mod(qs, index, p)) % p
        qs = next_qs
    return [(sol * x) % p for x in qsds]

# simplified, only to solve x^d ≡ a (mod p^k) by known x0's s.t. x0^r ≡ a (mod p), with (d, p) = 1
def hensel_lift(init_sol, d, a, p, k):
    sol = []
    p_tab = [p]
    for i in range(k - 1):
        p_tab.append(p_tab[i] * p)
    for x0 in init_sol:
        x = x0
        inv = inverse_mod((d * power_mod(x0, d - 1, p)) % p, p)
        for i in range(1, k):
            x += (inv * (a - power_mod(x, d, p_tab[i]))) % p_tab[i]
        sol.append(x)
    return sol

# only works for d = 2 and will be much slower if p or k is large
def dickson(x0, a, p, k):
    pk1 = p ^ (k - 1)
    pk = pk1 * p
    x = (power_mod(x0, pk1, pk) * power_mod(a, (pk - 2 * pk1 + 1) // 2, pk)) % pk
    return [x, pk - x]
