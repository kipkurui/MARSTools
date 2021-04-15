from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import map
from builtins import range
from builtins import object
class DerivVar(object):
    """Variable with derivatives

    Constructor: DerivVar(|value|, |index| = 0)

    Arguments:

    |value| -- the numerical value of the variable

    |index| -- the variable index (an integer), which serves to
               distinguish between variables and as an index for
               the derivative lists. Each explicitly created
               instance of DerivVar must have a unique index.

    Indexing with an integer yields the derivatives of the corresponding
    order.
    """

    def __init__(self, value, index=0, order=1):
        if order > 1:
            raise ValueError('Only first-order derivatives')
        self.value = value
        if order == 0:
            self.deriv = []
        elif type(index) == type([]):
            self.deriv = index
        else:
            self.deriv = index * [0] + [1]

    def __getitem__(self, item):
        if item < 0 or item > 1:
            raise ValueError('Index out of range')
        if item == 0:
            return self.value
        else:
            return self.deriv

    def __repr__(self):
        return repr((self.value, self.deriv))

    def __str__(self):
        return str((self.value, self.deriv))

    def __coerce__(self, other):
        if isDerivVar(other):
            return self, other
        else:
            return self, DerivVar(other, [])

    def __cmp__(self, other):
        return ((self.value > other.value) - (self.value < other.value)) 

    def __neg__(self):
        return DerivVar(-self.value, [-a for a in self.deriv])

    def __pos__(self):
        return self

    def __abs__(self):  # cf maple signum # derivate of abs
        absvalue = abs(self.value)
        return DerivVar(absvalue, list(map(lambda a, d=self.value / absvalue:
                                      d * a, self.deriv)))

    def __bool__(self):
        return self.value != 0

    def __add__(self, other):
        return DerivVar(self.value + other.value,
                        _mapderiv(lambda a, b: a + b, self.deriv, other.deriv))

    __radd__ = __add__

    def __sub__(self, other):
        return DerivVar(self.value - other.value,
                        _mapderiv(lambda a, b: a - b, self.deriv, other.deriv))

    def __rsub__(self, other):
        return DerivVar(other.value - self.value,
                        _mapderiv(lambda a, b: a - b, other.deriv, self.deriv))

    def __mul__(self, other):
        print("We got here")
        print(other)
        return DerivVar(self.value * other.value,
                        _mapderiv(lambda a, b: a + b,
                                  list(map(lambda x, f=other.value: f * x, self.deriv)),
                                  list(map(lambda x, f=self.value: f * x, other.deriv))))

    __rmul__ = __mul__

    def __div__(self, other):
        if not other.value:
            raise ZeroDivisionError('DerivVar division')
        inv = 1. / other.value
        return DerivVar(self.value * inv,
                        _mapderiv(lambda a, b: a - b,
                                  list(map(lambda x, f=inv: f * x, self.deriv)),
                                  list(map(lambda x, f=self.value * inv * inv: f * x,
                                      other.deriv))))

    def __rdiv__(self, other):
        return other / self

    def __pow__(self, other, z=None):
        if z is not None:
            raise TypeError('DerivVar does not support ternary pow()')
        val1 = pow(self.value, other.value - 1)
        val = val1 * self.value
        deriv1 = list(map(lambda x, f=val1 * other.value: f * x, self.deriv))
        if isDerivVar(other) and len(other.deriv) > 0:
            deriv2 = list(map(lambda x, f=val * Numeric.log(self.value): f * x,
                         other.deriv))
            return DerivVar(val, _mapderiv(lambda a, b: a + b, deriv1, deriv2))
        else:
            return DerivVar(val, deriv1)

    def __rpow__(self, other):
        return pow(other, self)

    def exp(self):
        v = Numeric.exp(self.value)
        return DerivVar(v, list(map(lambda x, f=v: f * x, self.deriv)))

    def log(self):
        v = Numeric.log(self.value)
        d = 1. / self.value
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def log10(self):
        v = Numeric.log10(self.value)
        d = 1. / (self.value * Numeric.log(10))
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def sqrt(self):
        v = Numeric.sqrt(self.value)
        d = 0.5 / v
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def sign(self):
        if self.value == 0:
            raise ValueError("can't differentiate sign() at zero")
        return DerivVar(Numeric.sign(self.value), 0)

    def sin(self):
        v = Numeric.sin(self.value)
        d = Numeric.cos(self.value)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def cos(self):
        v = Numeric.cos(self.value)
        d = -Numeric.sin(self.value)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def tan(self):
        v = Numeric.tan(self.value)
        d = 1. + pow(v, 2)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def sinh(self):
        v = Numeric.sinh(self.value)
        d = Numeric.cosh(self.value)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def cosh(self):
        v = Numeric.cosh(self.value)
        d = Numeric.sinh(self.value)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def tanh(self):
        v = Numeric.tanh(self.value)
        d = 1. / pow(Numeric.cosh(self.value), 2)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def arcsin(self):
        v = Numeric.arcsin(self.value)
        d = 1. / Numeric.sqrt(1. - pow(self.value, 2))
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def arccos(self):
        v = Numeric.arccos(self.value)
        d = -1. / Numeric.sqrt(1. - pow(self.value, 2))
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def arctan(self):
        v = Numeric.arctan(self.value)
        d = 1. / (1. + pow(self.value, 2))
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))

    def arctan2(self, other):
        den = self.value * self.value + other.value * other.value
        s = self.value / den
        o = other.value / den
        return DerivVar(Numeric.arctan2(self.value, other.value),
                        _mapderiv(lambda a, b: a - b,
                                  list(map(lambda x, f=o: f * x, self.deriv)),
                                  list(map(lambda x, f=s: f * x, other.deriv))))

    def gamma(self):
        from transcendental import gamma, psi
        v = gamma(self.value)
        d = v * psi(self.value)
        return DerivVar(v, list(map(lambda x, f=d: f * x, self.deriv)))


# Type check

def isDerivVar(x):
    "Returns 1 if |x| is a DerivVar object."
    return hasattr(x, 'value') and hasattr(x, 'deriv')


# Map a binary function on two first derivative lists

def _mapderiv(func, a, b):
    nvars = max(len(a), len(b))
    a = a + (nvars - len(a)) * [0]
    b = b + (nvars - len(b)) * [0]
    return list(map(func, a, b))


# Define vector of DerivVars

def DerivVector(x, y, z, index=0):
    """Returns a vector whose components are DerivVar objects.

    Arguments:

    |x|, |y|, |z| -- vector components (numbers)

    |index| -- the DerivVar index for the x component. The y and z
               components receive consecutive indices.
    """

    from Scientific.Geometry.Vector import Vector
    if isDerivVar(x) and isDerivVar(y) and isDerivVar(z):
        return Vector(x, y, z)
    else:
        return Vector(DerivVar(x, index),
                      DerivVar(y, index + 1),
                      DerivVar(z, index + 2))


def newtonRaphson(function, lox, hix, xacc, lista):
    """Finds the root of |function| which is bracketed by values
    |lox| and |hix| to an accuracy of +/- |xacc|. The algorithm
    used is a safe version of Newton-Raphson (see page 366 of NR in
    C, 2ed). |function| must be a function of one variable, and may
    only use operations defined for the DerivVar objects in the
    module FirstDerivatives.

    Example:

      >>>from Scientific.Functions.FindRoot import newtonRaphson
      >>>from Numeric import pi, sin, cos
      >>>def func(x):
      >>>    return (2*x*cos(x) - sin(x))*cos(x) - x + pi/4.0
      >>>newtonRaphson(func, 0.0, 1.0, 1.0e-12)

      yields '0.952847864655'.
    """

    maxit = 500
    tmp = function(DerivVar(lox), lista)
    fl = tmp[0]
    tmp = function(DerivVar(hix), lista)
    fh = tmp[0]
    if ((fl > 0.0 and fh > 0.0) or (fl < 0.0 and fh < 0.0)):
        ##        print "Calculo de lambda. No se encuentra solucion en [%f %f]" %(lox, hix)
        return None
    if (fl == 0.0): return lox
    if (fh == 0.0): return hix
    if (fl < 0.0):
        xl = lox
        xh = hix
    else:
        xh = lox
        xl = hix
    rts = 0.5 * (lox + hix)
    dxold = abs(hix - lox)
    dx = dxold
    tmp = function(DerivVar(rts), lista)
    f = tmp[0]
    df = tmp[1][0]
    for j in range(maxit):
        if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0)
            or (abs(2.0 * f) > abs(dxold * df))):
            dxold = dx
            dx = 0.5 * (xh - xl)
            rts = xl + dx
            if (xl == rts): return rts
        else:
            dxold = dx
            dx = f / df
            temp = rts
            rts = rts - dx
            if (temp == rts): return rts
        if (abs(dx) < xacc): return rts
        tmp = function(DerivVar(rts), lista)
        f = tmp[0]
        df = tmp[1][0]
        if (f < 0.0):
            xl = rts
        else:
            xh = rts
    print("Maximum number of iterations exceeded in newtonRaphson()")
    return 0.0
