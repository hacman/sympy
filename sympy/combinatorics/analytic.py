#!/usr/bin/env python
'''
Analytic Combinatorics
Flajolet, P. and R. Sedgewick. Analytic Combinatorics. Electronic edition.
Cambridge University Press, 2009.

Available for free in PDF form from
http://algo.inria.fr/flajolet/Publications/AnaCombi/anacombi.html
'''
from sympy.core.assumptions import ManagedProperties
from sympy.core import Basic, Dummy, Function
from sympy import oo, exp, log, Sum, solve, Integer
from sympy.ntheory import totient
from sympy.abc import z

class Multiton(ManagedProperties):
    '''
    Like a singleton, except there can be different labels on the class. Only
    one object with each label can be instantiated.

    Use it as a ``__metaclass__``.
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
    '''
    This is the basic construct in analytic combinatorics. As the name
    indicates, it represents a combinatorial class. Every other class
    in the analytic combinatorics module inherits from this. In particular,
    all classes have the ``.gf`` attribute to access their generating
    function.
    '''
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
                # FIXME: this eliminates use of SymPy functions in generating
                # function equations, too
                print
                print functions
                print self._gf
                print
                raise NotImplementedError("Systems of %d generating functions are not supported" % len(functions))
            # print 'F', self.name, functions
            if not isinstance(self._gf, Function):
                f = functions.pop()
                # print 'solve', self._gf
                solutions = solve(self._gf - f, f)
                if solutions:
                    # print 'hello', solutions

                    ## try to find "good" solutions: no terms z**n for n < 0
                    ## and all positive integer coefficients
                    for func in solutions:
                        good = True
                        ser = func.series(n=10)
                        _, exponent = ser.compute_leading_term(z).as_coeff_exponent(z)
                        if exponent < 0:
                            good = False
                            continue
                        for n in range(10):
                            c = ser.coeff(z, n)
                            if c < 0 or not c.is_Integer:
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
    '''
    An atom is a combinatorial class with a single element of size `1`.

    Examples
    ========
    >>> from sympy.combinatorics.analytic import CombinatorialAtom
    >>> a = CombinatorialAtom()
    >>> a.gf
    z
    '''

    __metaclass__ = Multiton

    def __new__(cls, arg='0'):
        obj = CombinatorialClass.__new__(cls, arg)
        obj._gf = z
        return obj

    def __str__(self):
        return 'Z_' + self.name
        ## TODO: would like to have just plain 'Z'

class NeutralClass(CombinatorialAtom):
    r'''
    A neutral class is a combinatorial class `\mathcal{E}` with a
    single element of size `0`.

    Examples
    ========
    >>> from sympy.combinatorics.analytic import NeutralClass
    >>> A = NeutralClass()
    >>> A.gf
    1
    '''

    def __new__(cls):
       obj = CombinatorialAtom.__new__(cls, 'E')
       obj._gf = Integer(1)
       return obj

class CombinatorialSum(CombinatorialClass):
    r'''
    A combinatorial sum is a class `\mathcal{C}` constructed from two existing
    classes by taking a disjoint union of their constituent elements. It is
    represented as `\mathcal{C} = \mathcal{A} + \mathcal{B}`. The corresponding
    operation on generating functions is pointwise addition: `A(z) + B(z)` which
    is the same as term-by-term addition:

    .. math ::
        C_N = A_N + B_N

    Examples
    ========
    >>> from sympy.combinatorics.analytic import NeutralClass, CombinatorialAtom
    >>> from sympy.combinatorics.analytic import CombinatorialSum
    >>> A = NeutralClass()
    >>> B = NeutralClass()
    >>> (A + B).gf
    2
    >>> C = CombinatorialAtom()
    >>> (C + A).gf
    z + 1
    '''

    def __new__(cls, *args):
        c = CombinatorialClass.__new__(cls, 'SUM')
        # print 'pre-sum', c._gf
        c._gf = reduce(lambda x, y: x._gf + y._gf, args)
        # print 'new sum', a._gf, b._gf
        c.summands = args
        return c

    def __str__(self):
        return str(self.a) + ' + ' + str(self.b)

class CombinatorialProduct(CombinatorialClass):
    r'''
    A combinatorial product is a class `\mathcal{C}` constructed from two
    existing classes by taking the cartesian product of their constituent
    elements. It is represented as `\mathcal{A} \times \mathcal{B}`. The
    corresponding operation on generating functions is pointwise multiplication
    `A(z) \cdot B(z)`, which is a convolution of the terms:

    .. math ::
        C_N = \sum_{k=0}^{N}{ A_k C_{N-k} }

    .. TODO: example
    '''
    def __new__(cls, a, b):
        ## TODO: accept any SymPy Integer
        if a == 1:
            return b
        elif b == 1:
            return a
        elif isinstance(a, int):
            b_list = [b] * a
            return CombinatorialSum(*b_list)
        elif isinstance(b, int):
            a_list = [a] * b
            return CombinatorialSum(*a_list)

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
    r'''
    A sequence is a compound combinatorial class, defined as finite ordered
    sequences of objects from another combinatorial class. Unlike with
    mathematical sets, two sequences with the same elements present in a
    different order from each other are considered distinct sequences.
    Sequences of a given length `k` are simply the product of `k` copies of
    a class with itself:

    .. math ::
        \operatorname{SEQ}_k (\mathcal{A}) = \underbrace{\mathcal{A} \times \cdots \times \mathcal{A}}_{k\: \mathrm{ copies}}

    succinctly noted as `\mathcal{A}^k`. The corresponding operation on
    generating functions is exponentiation: `A(z)^k`. The class of all sequences
    is defined as

    .. math ::
        \operatorname{SEQ} (\mathcal{A}) = \mathcal{E} + \mathcal{A} + \mathcal{A}^2 + \mathcal{A}^3 + \cdots

    with the generating function

    .. math ::
         1 + A(z) + A(z)^2 + A(z)^3 + \cdots

    which is the geometric series in `A(z)`. This can be rewritten as

    .. math ::
         \frac{1}{1- A(z)}

    It should be noted that in the sequence construction, the argument
    `\mathcal{A}` can contain no element of size zero, otherwise the class fails
    to be properly defined (it does not satisfy the finiteness condition).
    Unfortunately, there is no simple formula relating the coefficients of a
    sequence class to those of class from which it was constructed.

     .. TODO: example
    '''
    def __new__(cls, other, k=None):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        c.k = k
        if not c.k:
            c._gf = 1/(1-c.a.gf)
        else:
            c._gf = c.a.gf ** c.k
        return c

    def __str__(self):
        return 'SEQ(' + str(self.a) + ')'

class MSET(CombinatorialClass):
    r'''
    A *multi-set* is a generalization of the mathematical set which allows
    multiple copies of the same element of the set. Alternatively, one could
    define it as a set whose elements have a multiplicity (which would always
    be `1` for ordinary sets). Unlike with sequences (which also admit multiple
    copies of elements), order does not matter in multi-sets. The multi-set is
    defined as a set of equivalence classes of sequences, where two sequences
    are in the same class if they can be mapped to each other by a permutation.
    The generating function for a multi-set formed from elements of a
    combinatorial class `\mathcal{A}` is

    .. math ::
        \exp\left(\sum_{k=1}^{\infty}{\frac{A(z^k)}{k}}\right).

    An equivalent construction of the generating function is

    .. math ::
        \prod_{N\geq 1}{(1-z^N)^{-A_N}}

    where `A_N` are the coefficients of `A(z)`.

     .. TODO: example
    '''
    def __new__(cls, other, k=None):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        c.k = k
        if not c.k:
            m = Dummy('m', integer=True)
            c._gf = exp(Sum(c.a._gf.subs(z, z**m)/m, (m, 1, oo)))
        else:
            c._gf = exp(c.a._gf.subs(z, z**k)/k)
            # raise NotImplementedError('is this right?')

        return c

    def __str__(self):
        return 'MSET(' + str(self.a) + ')'

class PSET(CombinatorialClass):
    r'''
    A *power set* (or simply a *set*) is a combinatorial class formed from sets
    of elements of another class. Like the multi-set, order does not matter
    here, but unlike the multi-set, no multiplicity is allowed in a power set.
    The generating function of `\operatorname{PSET}(\mathcal{A})` is

    .. math ::
        \exp\left(\sum_{k=1}^{\infty}{\frac{(-1)^{k-1} A(z^k)}{k}}\right)

    The generating function transformation is similar to that of a multi-set,
    but the sign of the terms in the series alternates. Growth is slower,
    which reflects the fact that `\operatorname{PSET}(\mathcal{A}) \subset \operatorname{MSET}(\mathcal{A})`
    so there are fewer power sets of a given size than there are multi-sets.

    An equivalent construction of the generating function is

    .. math ::
        \prod_{N\geq 1}{(1+z^N)^{A_N}}

    where `A_N` are the coefficients of `A(z)`.

     .. TODO: example
    '''
    def __new__(cls, other):
        c = CombinatorialClass.__new__(cls)
        c.a = other
        k = Dummy('k', integer=True)
        c._gf = exp(Sum((-1)**(k-1)*c.a._gf.subs(z, z**k)/k, (k, 1, oo)))
        return c

    def __str__(self):
        return 'PSET(' + str(self.a) + ')'

class CYC(CombinatorialClass):
    r'''
    Cycles are constructed from sequences in a similar fashion to multi-sets.
    Whereas multi-sets were equivalent if they were equal under any permutation,
    cycles are equivalent only under cyclic permutations. The transformation of
    generating functions is
    
    .. math ::
        \sum_{k = 1}^{\infty}{\frac{\phi(k)}{k} \log\frac{1}{1-A(z^k)}}

    where `\phi(k)` is Euler's :class:`totient function <sympy.ntheory.factor_.totient>`.

     .. TODO: examples
    '''
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
