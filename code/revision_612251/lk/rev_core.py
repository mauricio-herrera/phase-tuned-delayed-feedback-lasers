"""Core routines for the JOSA B revision.

ECMs are parameterized by the total return phase theta, which enumerates *all*
external-cavity modes and antimodes: phi = theta + kappa*tau*(alpha cos theta + sin theta).
Stability is decided by an independent Mikhailov (argument-principle) count along
Re(lambda) = -mu, applied to g(lambda) = det Delta(lambda)/lambda, which removes the
neutral phase root.
"""
import numpy as np
from scipy.optimize import brentq

B_FIELD = np.array([[1., 0.], [0., 1.], [0., 0.]])


def ecm_theta(theta, kappa, tau, alpha=4.0, P=1.0):
    c, s = np.cos(theta), np.sin(theta)
    Ns = -kappa * c
    A2 = (P - Ns) / (1 + 2 * Ns)
    Om = -kappa * (alpha * c + s)
    phi = theta - Om * tau
    F = 1 + kappa * tau * (c - alpha * s)
    return dict(theta=theta, phi=phi, Ns=Ns, A=np.sqrt(A2), Om=Om, F=F)


def mats(e, kappa, alpha=4.0, eps=0.005):
    th, A, Ns = e['theta'], e['A'], e['Ns']
    cr, ci = kappa * np.cos(th), -kappa * np.sin(th)
    A0 = np.array([[-cr, ci, A], [-ci, -cr, alpha * A],
                   [-2 * eps * (1 + 2 * Ns) * A, 0., -eps * (1 + 2 * A * A)]])
    A1 = np.array([[cr, -ci, 0.], [ci, cr, 0.], [0., 0., 0.]])
    return A0, A1


def det_delta(lam, A0, A1, tau):
    lam = np.asarray(lam, complex)
    z = np.exp(-lam * tau)
    M = lam[:, None, None] * np.eye(3) - A0 - z[:, None, None] * A1
    return np.linalg.det(M)


def count_right(A0, A1, tau, mu=0.0, R=6.0):
    """Number of non-neutral characteristic roots with Re(lambda) > -mu."""
    n = 2 * (int(max(4000, 40 * R * tau / np.pi + 2000)) // 2)  # even: avoids lambda = 0 on the line
    for _ in range(6):
        w = np.linspace(-R, R, n)
        lam = -mu + 1j * w
        g = det_delta(lam, A0, A1, tau) / lam
        ph = np.unwrap(np.angle(g))
        if np.max(np.abs(np.diff(ph))) < 0.4:
            break
        n *= 2
    N = 1 - (ph[-1] - ph[0]) / (2 * np.pi)
    return int(np.rint(N))


def margin(A0, A1, tau, mu_max=None, tol=2e-5):
    """m = -max_{lambda != 0} Re(lambda), by bisection on the Mikhailov count.
    Returns (m, capped) where capped=True means m >= mu_max."""
    if mu_max is None:
        mu_max = min(0.06, 5.0 / tau)
    if count_right(A0, A1, tau, 0.0) > 0:
        # unstable: bisect on positive real part
        lo, hi = 0.0, 0.2
        if count_right(A0, A1, tau, -hi) > 0:
            return -hi, True
        while hi - lo > tol:
            mid = .5 * (lo + hi)
            if count_right(A0, A1, tau, -mid) > 0: lo = mid
            else: hi = mid
        return -.5 * (lo + hi), False
    if count_right(A0, A1, tau, mu_max) == 0:
        return mu_max, True
    lo, hi = 0.0, mu_max
    while hi - lo > tol:
        mid = .5 * (lo + hi)
        if count_right(A0, A1, tau, mid) == 0: lo = mid
        else: hi = mid
    return .5 * (lo + hi), False


def omega_grid(tau, wmax=200.0):
    dw = min(2e-4, 0.03 / tau)
    w1 = np.arange(1e-7, 0.6, dw)
    w2 = np.geomspace(0.6, wmax, 1200)[1:]
    return np.concatenate([w1, w2])


def variance(A0, A1, tau, C, B=B_FIELD, om=None, wmax=None):
    if om is None:
        om = omega_grid(tau)
    if wmax is not None:
        om = om[om <= wmax]
    D = 1j * om[:, None, None] * np.eye(3) - A0 - np.exp(-1j * om * tau)[:, None, None] * A1
    H = np.linalg.solve(D, np.broadcast_to(B, (len(om),) + B.shape))
    G = np.einsum('ab,wbc->wac', np.atleast_2d(C), H)
    psd = np.real(np.einsum('wac,wac->w', G, G.conj()))
    return np.trapezoid(psd, om) / np.pi


class Model:
    def __init__(self, tau, alpha=4.0, eps=0.005, P=1.0):
        self.tau, self.alpha, self.eps, self.P = tau, alpha, eps, P
        self.om = omega_grid(tau)
        e0 = ecm_theta(0.0, 0.0, tau, alpha, P)
        A0, A1 = mats(e0, 0.0, alpha, eps)
        self.R0 = variance(A0, A1, tau, [2 / e0['A'], 0, 0], om=self.om)
        self.N0 = variance(A0, A1, tau, [0, 0, 1], om=self.om)

    def state(self, theta, kappa, stability=True, carrier=False):
        e = ecm_theta(theta, kappa, self.tau, self.alpha, self.P)
        A0, A1 = mats(e, kappa, self.alpha, self.eps)
        e['Dratio'] = self.P / (e['A'] ** 2 * e['F'] ** 2)
        e['Rratio'] = variance(A0, A1, self.tau, [2 / e['A'], 0, 0], om=self.om) / self.R0
        if carrier:
            e['Nratio'] = variance(A0, A1, self.tau, [0, 0, 1], om=self.om) / self.N0
        if stability:
            e['nunst'] = count_right(A0, A1, self.tau, 0.0)
        return e

    def theta_of_phi(self, phi, kappa, branch_center=None):
        """All ECM/antimode solutions theta for a given phi (mod 2 pi)."""
        kt = kappa * self.tau
        th = np.linspace(-np.pi, np.pi, 4001)
        def f(t):
            d = t + kt * (self.alpha * np.cos(t) + np.sin(t)) - phi
            return (d + np.pi) % (2 * np.pi) - np.pi
        v = f(th); out = []
        for i in range(len(th) - 1):
            if v[i] == 0: out.append(th[i])
            elif v[i] * v[i + 1] < 0 and abs(v[i] - v[i + 1]) < np.pi:
                out.append(brentq(f, th[i], th[i + 1], xtol=1e-13))
        return np.array(out)


def wrap(x):
    return (x + np.pi) % (2 * np.pi) - np.pi
