#!/usr/bin/env python
'''
Analytic Combinatorics
Flajolet, P. and R. Sedgewick. Analytic Combinatorics. Electronic edition.
Cambridge University Press, 2009.
http://algo.inria.fr/flajolet/Publications/book.pdf
'''
from sympy.core.assumptions import ManagedProperties
from sympy.core import Basic, Dummy, Function
from sympy import oo, exp, log, Sum, solve, Integer
from sympy.ntheory import totient
from sympy.abc import z

class Multiton(ManagedProperties):
    '''
    Like a singleton, except there can be different labels on the class. Only
    one object with each name can be instantiated.

    use it as a __metaclass__
    '''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if not args:
            key = (cls, '0')
        else:
            key = (cls, str(args[0]))

        if key not in cls._instances:
            cls._instances[key] = super(Multiton, cls).__call__()
            if args:
                cls._instances[key].name = str(args[0])

        return cls._instances[key]

class CombinatorialClass(Basic):
    def __new__(cls, arg='C', gf=None):
        obj = Basic.__new__(cls)
        obj._gf = Function(arg)(z)
        obj.name = arg
        obj._did_compute = False
        return obj

    @property
    def gf(self):
        '''
        Generating function for the combinatorial class
        '''
        if not self._did_compute:
            self._gf_compute()
        return self._gf

    def _gf_compute(self):
        self._did_compute = True
        if self._gf.has(Function):
            functions = self._gf.find(Function)
            if len(functions) > 1:
                raise NotImplementedError("Systems of %d generating functions are not supported" % len(functions))
            # print 'F', self.name, functions
            if not isinstance(self._gf, Function):
                f = functions.pop()
                # print 'solve', self._gf
                solutions = solve(self._gf - f, f)
                if solutions:
                    # print 'hello', solutions
                    for func in solutions:
                        good = True
                        ser = func.series(n=10)
                        if ser.coeff(1/z) != 0:
                            good = False
                            continue
                        for n in range(10):
                            if ser.coeff(z**n) < 0:
                                good = False
                                break

                        if good:
                            self._gf = func

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

    def __str__(self):
        return self.name

class CombinatorialAtom(CombinatorialClass):
    __metaclass__ = Multiton

    def __new__(cls, arg='0'):
        obj = CombinatorialClass.__new__(cls, arg)
        obj._gf = z
        return obj

    def __str__(self):
        return 'Z_' + self.name
        ## TODO: would like to have just plain 'Z'

class NeutralClass(CombinatorialAtom):
    def __new__(cls):
       obj = CombinatorialAtom.__new__(cls, 'E')
       obj._gf = Integer(1)
       return obj

class CombinatorialSum(CombinatorialClass):
    def __new__(cls, a, b):
        c = CombinatorialClass.__new__(cls, 'SUM')
        # print 'pre-sum', c._gf
        c._gf = a._gf + b._gf
        # print 'new sum', a._gf, b._gf
        c.a = a
        c.b = b
        return c

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b)

class CombinatorialProduct(CombinatorialClass):
    def __new__(cls, a, b):
        c = CombinatorialClass.__new__(cls, 'PROD')
        # print 'pre-product', c._gf
        c._gf = a._gf * b._gf
        # print 'new product', a._gf, b._gf
        c.a = a
        c.b = b
        return c

    def __str__(self):
        return str(self.a) + ' x ' + str(self.b)

class SEQ(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        c._gf = 1/(1-c.a._gf)
        return c

    def __str__(self):
        return 'SEQ(' + str(self.a) + ')'

class MSET(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = exp(Sum(c.a._gf.subs(z, z**k)/k, (k, 1, oo)))
        return c

    def __str__(self):
        return 'MSET(' + str(self.a) + ')'

class PSET(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = exp(Sum((-1)**(k-1)*c.a._gf.subs(z, z**k)/k, (k, 1, oo)))
        return c

    def __str__(self):
        return 'PSET(' + str(self.a) + ')'

class CYC(CombinatorialClass):
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = Sum(totient(k)/k * log(1/(1-c.a._gf.subs(z,z**k))), (k, 1, oo))
        return c

    def __str__(self):
        return 'CYC(' + str(self.a) + ')'

# restricted constructions (p.30)
# labelled structures
