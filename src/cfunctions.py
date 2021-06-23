# This file was automatically generated by SWIG (http://www.swig.org).
# Version 4.0.2
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info < (2, 7, 0):
    raise RuntimeError("Python 2.7 or later required")

# Import the low-level C/C++ module
if __package__ or "." in __name__:
    from . import _cfunctions
else:
    import _cfunctions

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_instance_variable(set):
    def set_instance_attr(self, name, value):
        if name == "thisown":
            self.this.own(value)
        elif name == "this":
            set(self, name, value)
        elif hasattr(self, name) and isinstance(getattr(type(self), name), property):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add instance attributes to %s" % self)
    return set_instance_attr


def _swig_setattr_nondynamic_class_variable(set):
    def set_class_attr(cls, name, value):
        if hasattr(cls, name) and not isinstance(getattr(cls, name), property):
            set(cls, name, value)
        else:
            raise AttributeError("You cannot add class attributes to %s" % cls)
    return set_class_attr


def _swig_add_metaclass(metaclass):
    """Class decorator for adding a metaclass to a SWIG wrapped class - a slimmed down version of six.add_metaclass"""
    def wrapper(cls):
        return metaclass(cls.__name__, cls.__bases__, cls.__dict__.copy())
    return wrapper


class _SwigNonDynamicMeta(type):
    """Meta class to enforce nondynamic attributes (no new attributes) for a class"""
    __setattr__ = _swig_setattr_nondynamic_class_variable(type.__setattr__)



def calcMSD(r, ct, x, y, z, N, n_t):
    return _cfunctions.calcMSD(r, ct, x, y, z, N, n_t)

def logMSD(r, x, y, z, N, n_t):
    return _cfunctions.logMSD(r, x, y, z, N, n_t)

def calcM4D(r, ct, x, y, z, N, n_t):
    return _cfunctions.calcM4D(r, ct, x, y, z, N, n_t)

def logM4D(r, x, y, z, N, n_t):
    return _cfunctions.logM4D(r, x, y, z, N, n_t)

def RDF1(g, r, nbins, N, rho, L):
    return _cfunctions.RDF1(g, r, nbins, N, rho, L)

def RDF2(g, r, nbins, N, rho, L, Nc):
    return _cfunctions.RDF2(g, r, nbins, N, rho, L, Nc)

def RDF3(g, r, nbins, N, rho, L, Nc):
    return _cfunctions.RDF3(g, r, nbins, N, rho, L, Nc)

def RDF4(g, r, nbins, N, rho, L, N2):
    return _cfunctions.RDF4(g, r, nbins, N, rho, L, N2)

def fcn(f, r_RDF, g, nbins, L, rho):
    return _cfunctions.fcn(f, r_RDF, g, nbins, L, rho)

def lam(rt, r, ct, x, y, z, q, N, n_t):
    return _cfunctions.lam(rt, r, ct, x, y, z, q, N, n_t)

def loglam(r, x, y, z, q, N, n_t):
    return _cfunctions.loglam(r, x, y, z, q, N, n_t)

def lam_cross1(rt, r, ct, x, y, z, q, N, n_t):
    return _cfunctions.lam_cross1(rt, r, ct, x, y, z, q, N, n_t)

def loglam_cross1(r, x, y, z, q, N, n_t):
    return _cfunctions.loglam_cross1(r, x, y, z, q, N, n_t)

def lam_cross2(rt, r, ct, x, y, z, q, N, n_t):
    return _cfunctions.lam_cross2(rt, r, ct, x, y, z, q, N, n_t)

def loglam_cross2(r, x, y, z, q, N, n_t):
    return _cfunctions.loglam_cross2(r, x, y, z, q, N, n_t)

def lam_cross3(rt, r, ct, x, y, z, q, N, n_t):
    return _cfunctions.lam_cross3(rt, r, ct, x, y, z, q, N, n_t)

def loglam_cross3(r, x, y, z, q, N, n_t):
    return _cfunctions.loglam_cross3(r, x, y, z, q, N, n_t)

def fqt(fqt, ct, x, y, z, qv, N, n_t, num_q):
    return _cfunctions.fqt(fqt, ct, x, y, z, qv, N, n_t, num_q)

def logfqt(fqt, x, y, z, qv, N, n_t, num_q):
    return _cfunctions.logfqt(fqt, x, y, z, qv, N, n_t, num_q)

def SSF1(ssf, x, y, z, qv, nq, N, n_t, num_p):
    return _cfunctions.SSF1(ssf, x, y, z, qv, nq, N, n_t, num_p)

def SSF2(ssf, x, y, z, qv, nq, N1, N2, n_t, num_p):
    return _cfunctions.SSF2(ssf, x, y, z, qv, nq, N1, N2, n_t, num_p)

def logvanHove(g, x, y, z, N1, N2, nbins, n_t, L):
    return _cfunctions.logvanHove(g, x, y, z, N1, N2, nbins, n_t, L)

def vanHove(g, ct, x, y, z, N1, N2, nbins, n_t, L):
    return _cfunctions.vanHove(g, ct, x, y, z, N1, N2, nbins, n_t, L)

def vanHove_short(g, x, y, z, N1, N2, nbins, n_t, n_t2, L):
    return _cfunctions.vanHove_short(g, x, y, z, N1, N2, nbins, n_t, n_t2, L)

def MUAC(muac, mu_int, ct, N, n_t):
    return _cfunctions.MUAC(muac, mu_int, ct, N, n_t)

def MUAC_log(muac, mu_int, N, n_t):
    return _cfunctions.MUAC_log(muac, mu_int, N, n_t)

def pofn(p1, r, N1, N2, lo, hi, L):
    return _cfunctions.pofn(p1, r, N1, N2, lo, hi, L)

def pofM(p1, r, N1, N2, lo, hi, L, Nc):
    return _cfunctions.pofM(p1, r, N1, N2, lo, hi, L, Nc)

def pair_lj_dipole(u, r, mu, ep, sig, L, N):
    return _cfunctions.pair_lj_dipole(u, r, mu, ep, sig, L, N)

def pair_lj_dipole_cat(u, r, mu, rc, ep, sig, L, N, q):
    return _cfunctions.pair_lj_dipole_cat(u, r, mu, rc, ep, sig, L, N, q)

def pair_fene(u, r, ep, k_fene, R0, L, N, M):
    return _cfunctions.pair_fene(u, r, ep, k_fene, R0, L, N, M)

def pair_angle(u, r, K, L, N, M):
    return _cfunctions.pair_angle(u, r, K, L, N, M)

def find_shell(index, r, rc, L, N):
    return _cfunctions.find_shell(index, r, rc, L, N)

def pair_fene_shell(u, r, index, ep, k_fene, R0, L, N, M, bnd):
    return _cfunctions.pair_fene_shell(u, r, index, ep, k_fene, R0, L, N, M, bnd)

def pair_angle_shell(u, r, index, K, L, N, M, ang):
    return _cfunctions.pair_angle_shell(u, r, index, K, L, N, M, ang)

def pair_lj_dipole_shell(u, r, index, mu, ep, sig, L, N):
    return _cfunctions.pair_lj_dipole_shell(u, r, index, mu, ep, sig, L, N)

