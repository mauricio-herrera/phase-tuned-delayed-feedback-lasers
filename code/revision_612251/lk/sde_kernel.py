"""Euler-Maruyama integrator of the stochastic Lang-Kobayashi equations in the frame
rotating with a selected external-cavity mode (ECM).

    dE/dt = (1 + i alpha) N E - i Omega E + kappa e^{-i theta} E(t - tau) + sigma xi(t)
    dN/dt = eps [P - N - (1 + 2N)|E|^2]

The ECM is specified by its total return phase theta (rev_core.ecm_theta), so the chosen
mode is followed even when several ECMs coexist. Same scheme and noise convention as the
original simulate_nonlinear_dde.py (each field quadrature receives sigma*sqrt(dt)*N(0,1)).

Per realization the kernel returns
  q   : unwrapped optical phase, sampled every `stride` steps after burn-in
  sI, sI2, sN, sN2, n : running sums of (|E|^2 - A^2) and (N - Ns) after burn-in (every step)
Statistics are accumulated about the ECM values to avoid cancellation errors.
"""
import numpy as np

try:
    from numba import njit, prange
    HAVE_NUMBA = True
except ImportError:  # pragma: no cover
    HAVE_NUMBA = False


def _kernel_py(A, Ns, Om, theta, kappa, alpha, eps, P, sigma, dt, nd, nburn, nsteps, stride, seeds):
    n_ens = len(seeds)
    nrec = (nsteps - nburn) // stride
    q = np.zeros((nrec, n_ens)); stats = np.zeros((n_ens, 5))
    rng = np.random.default_rng(int(seeds[0]))
    hist = np.full((nd + 1, n_ens), A + 0j)
    z = np.full(n_ens, A + 0j); N = np.full(n_ens, Ns)
    fb = kappa * np.exp(-1j * theta); rot = 1 + 1j * alpha
    sq = sigma * np.sqrt(dt); ptr = nd; ph = np.zeros(n_ens); k = 0
    for step in range(1, nsteps + 1):
        zdel = hist[(ptr + 1) % (nd + 1)] if nd > 0 else z
        znew = z + dt * ((rot * N - 1j * Om) * z + fb * zdel) + sq * (rng.standard_normal(n_ens) + 1j * rng.standard_normal(n_ens))
        I = z.real ** 2 + z.imag ** 2
        N = N + dt * eps * (P - N - (1 + 2 * N) * I)
        ph += np.angle(znew * np.conj(z))
        ptr = (ptr + 1) % (nd + 1); hist[ptr] = znew; z = znew
        if step > nburn:
            dI = z.real ** 2 + z.imag ** 2 - A * A; dN = N - Ns
            stats[:, 0] += dI; stats[:, 1] += dI * dI; stats[:, 2] += dN; stats[:, 3] += dN * dN; stats[:, 4] += 1
            if (step - nburn) % stride == 0 and k < nrec:
                q[k] = ph; k += 1
    return q, stats


if HAVE_NUMBA:
    @njit(parallel=True, fastmath=False, cache=True)
    def _kernel_nb(A, Ns, Om, theta, kappa, alpha, eps, P, sigma, dt, nd, nburn, nsteps, stride, seeds):
        n_ens = len(seeds)
        nrec = (nsteps - nburn) // stride
        q = np.zeros((nrec, n_ens)); stats = np.zeros((n_ens, 5))
        fbr = kappa * np.cos(theta); fbi = -kappa * np.sin(theta)
        sq = sigma * np.sqrt(dt)
        for j in prange(n_ens):
            np.random.seed(seeds[j])
            L = nd + 1
            hr = np.full(L, A); hi = np.zeros(L)
            xr = A; xi = 0.0; N = Ns; ph = 0.0; ptr = nd; k = 0
            s0 = 0.0; s1 = 0.0; s2 = 0.0; s3 = 0.0; cnt = 0.0
            for step in range(1, nsteps + 1):
                if nd > 0:
                    idx = (ptr + 1) % L; dr = hr[idx]; di = hi[idx]
                else:
                    dr = xr; di = xi
                # (1 + i alpha) N - i Om  times z
                ar = N; ai = alpha * N - Om
                fr = ar * xr - ai * xi + fbr * dr - fbi * di
                fi = ar * xi + ai * xr + fbr * di + fbi * dr
                nr = xr + dt * fr + sq * np.random.standard_normal()
                ni = xi + dt * fi + sq * np.random.standard_normal()
                I = xr * xr + xi * xi
                N = N + dt * eps * (P - N - (1.0 + 2.0 * N) * I)
                # exact phase increment arg(znew * conj(z))
                ph += np.arctan2(ni * xr - nr * xi, nr * xr + ni * xi)
                ptr = (ptr + 1) % L; hr[ptr] = nr; hi[ptr] = ni
                xr = nr; xi = ni
                if step > nburn:
                    dI = xr * xr + xi * xi - A * A; dN = N - Ns
                    s0 += dI; s1 += dI * dI; s2 += dN; s3 += dN * dN; cnt += 1.0
                    if (step - nburn) % stride == 0 and k < nrec:
                        q[k, j] = ph; k += 1
            stats[j, 0] = s0; stats[j, 1] = s1; stats[j, 2] = s2; stats[j, 3] = s3; stats[j, 4] = cnt
        return q, stats


def integrate(e, kappa, tau, alpha=4.0, eps=0.005, P=1.0, sigma=5e-4, dt=0.02,
              burn=1000.0, T=10000.0, rec=1.0, seeds=(1,), backend='auto'):
    """Integrate len(seeds) independent realizations; returns (q, stats, rec_dt)."""
    nd = int(round(tau / dt)) if kappa > 0 else 0
    if kappa > 0 and abs(nd * dt - tau) > 1e-9:
        raise ValueError('tau/dt must be an integer')
    nburn = int(round(burn / dt)); nsteps = nburn + int(round(T / dt)); stride = max(1, int(round(rec / dt)))
    seeds = np.asarray(seeds, dtype=np.int64)
    args = (float(e['A']), float(e['Ns']), float(e['Om']), float(e['theta']), float(kappa), float(alpha),
            float(eps), float(P), float(sigma), float(dt), nd, nburn, nsteps, stride, seeds)
    use_nb = HAVE_NUMBA if backend == 'auto' else (backend == 'numba')
    if use_nb and not HAVE_NUMBA:
        raise RuntimeError('numba not installed')
    q, stats = (_kernel_nb if use_nb else _kernel_py)(*args)
    return q, stats, stride * dt


def estimators(q, stats, A, rec_dt, lags):
    """Per-realization phase-diffusion, RIN and carrier-variance estimates.
    D_j: half the least-squares slope of the phase mean-square displacement over `lags`
    (overlapping increments). RIN_j = var(|E|^2)/A^4; N_j = var(N)."""
    L = np.unique(np.maximum(1, np.rint(np.asarray(lags) / rec_dt).astype(int)))
    L = L[L < q.shape[0] // 4]
    lt = L * rec_dt
    msd = np.array([np.mean((q[l:] - q[:-l]) ** 2, axis=0) for l in L])  # (nlag, n_ens)
    X = np.vstack([lt, np.ones_like(lt)]).T
    coef, *_ = np.linalg.lstsq(X, msd, rcond=None)
    D = coef[0] / 2
    n = stats[:, 4]
    varI = stats[:, 1] / n - (stats[:, 0] / n) ** 2
    varN = stats[:, 3] / n - (stats[:, 2] / n) ** 2
    return D, varI / A ** 4, varN, lt
