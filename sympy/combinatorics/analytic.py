#!/usr/bin/env python
'''
Analytic Combinatorics
Sedgewick, R. and P. Flajolet. Analytic Combinatorics. Cambridge University
Press, 2007.
'''

from sympy.core.assumptions import ManagedProperties
from sympy.abc import z
from sympy.core import Basic, Dummy
from sympy import oo, exp, log, Sum, totient


class Multiton(ManagedProperties):
    '''
    Like a singleton, except there can be different labels on the class. Only
    one object with each name can be instantiated.

    use it as a __metaclass__
    '''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if not args:
            args = (0, )
        key = (cls, str(args[0]))

        if key not in cls._instances:
            cls._instances[key] = super(Multiton, cls).__call__()
            cls._instances[key].name = str(args[0])

        return cls._instances[key]

class CombinatorialClass(Basic):
    @property
    def gf(self):
        return self._gf

    def __add__(self, other):
        return CombinatorialSum(self, other)

    def __mul__(self, other):
        return CombinatorialProduct(self, other)

    def __pow__(self, n):
        if n == 2:
            return self * self
        elif n == 3:
            return self * self * self
        elif n == 4:
            return self * self * self * self
        else:
            raise NotImplementedError('pow')

class CombinatorialAtom(CombinatorialClass):
    __metaclass__ = Multiton
    _gf = z

    def __str__(self):
        return 'Z_' + self.name

    def __repr__(self):
        return str(self)

class Empty(CombinatorialAtom):
    _gf = 1

class CombinatorialSum(CombinatorialClass):
    def __new__(cls, a, b):
        c = CombinatorialClass.__new__(cls)
        c._gf = a.gf + b.gf
        c.a = a
        c.b = b
        return c

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b)

class CombinatorialProduct(CombinatorialClass):
    def __new__(cls, a, b):
        c = CombinatorialClass.__new__(cls)
        c._gf = a.gf * b.gf
        c.a = a
        c.b = b
        return c

    def __str__(self):
        return str(self.a) + ' x ' + str(self.b)

class SEQ(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        c._gf = 1/(1-c.a.gf)
        return c

    def __str__(self):
        return 'SEQ(' + str(self.a) + ')'

class MSET(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = exp(Sum(c.a.gf.subs(z, z**k)/k, (k, 1, oo)))
        return c

    def __str__(self):
        return 'MSET(' + str(self.a) + ')'

class PSET(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = exp(Sum((-1)**(k-1)*c.a.gf.subs(z, z**k)/k, (k, 1, oo)))
        return c

    def __str__(self):
        return 'PSET(' + str(self.a) + ')'

class CYC(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = Sum(totient(k)/k) * log(1/(1-c.a.gf.subs(z,z**k)), (k, 1, oo))
        return c

    def __str__(self):
        return 'CYC(' + str(self.a) + ')'

# restricted constructions (p.30)
