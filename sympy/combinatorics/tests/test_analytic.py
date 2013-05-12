from sympy.combinatorics.analytic import (Multiton,
            CombinatorialClass,
            CombinatorialAtom, CombinatorialSum, SEQ,
            CombinatorialProduct, CYC, MSET,
            NeutralClass,
            )
from sympy.abc import z
from sympy import catalan, fibonacci
from sympy.utilities.pytest import XFAIL

def test_multiton():
    class MyAtom(object):
        __metaclass__ = Multiton
    assert MyAtom() is MyAtom()
    assert MyAtom(0) is MyAtom()
    assert MyAtom(1) is MyAtom(1)
    assert MyAtom('a') is MyAtom('a')

def test_Atom():
    assert CombinatorialAtom(0) is CombinatorialAtom(0)
    assert CombinatorialAtom(1) is CombinatorialAtom(1)

    one = CombinatorialAtom(1)
    zero = CombinatorialAtom(0)
    assert CombinatorialAtom(0) is zero
    assert CombinatorialAtom(1) is one
    assert CombinatorialAtom() is zero
    assert str(one) == 'Z_1', str(one)

    assert zero.gf == z
    assert one.gf == z

def test_Sum():
    one = CombinatorialAtom(1)
    zero = CombinatorialAtom(0)
    s = one + zero
    assert s.gf == 2*z

@XFAIL
def test_int_Product():
    Z = CombinatorialAtom()
    CC = Z*10

@XFAIL
def test_Product_from_power():
    Z = CombinatorialAtom()
    CC = Z**9

def test_Product():
    Z = CombinatorialAtom()
    CC = Z*Z
    assert CC.gf == z**2
    NN= NeutralClass()
    assert (CC * NN).gf == CC.gf
    assert (NN * CC).gf == CC.gf

def test_Seq():
    Z = CombinatorialAtom(0)
    seq = SEQ(Z)
    assert seq.gf == 1/(1-z)

def test_Cyc():
    Z = CombinatorialAtom(0)
    cyc = CYC(Z)

def test_restricted_Seq():
    Z = CombinatorialAtom(0)
    seq = SEQ(Z, 2)
    assert seq.gf == z**2
    seq = SEQ(Z, 3)
    assert seq.gf == z**3
    seq = SEQ(Z, 3000)
    assert seq.gf == z**3000

## examples from the documentation's "preview" of the symbolic method

def test_integers():
    '''
    A sequence of atoms defines the integers.
    '''
    Z = CombinatorialAtom(0)
    I = SEQ(Z)
    ser = I.gf.series(z, 25)
    [ser.coeff(z, k) for k in range(25)] == range(1, 25)

def test_binary_strings():
    '''
    Counting binary strings gives 2**n
    '''
    Z0 = CombinatorialAtom('0')
    Z1 = CombinatorialAtom('1')
    B = CombinatorialClass('B')
    ## a binary string is either Empty or it's a zero or a one followed by
    ## a binary string
    B = NeutralClass() + (Z0 + Z1) * B
    ser = B.gf.series(n=32)
    for k in range(32):
        assert ser.coeff(z, k) == 2**k

def test_binary_strings_no_00():
    '''
    Counting binary strings with no consecutive 0s gives the Fibonacci numbers
    '''
    Z0 = CombinatorialAtom('0')
    Z1 = CombinatorialAtom('1')
    B00 = CombinatorialClass('B_00')
    ## a binary string with no 00 is either:
    ##   - Empty
    ##   - a zero
    ##   - a one or zero-one followed by a binary string with no 00
    B00 = NeutralClass() + Z0 + (Z1 + Z0 * Z1) * B00
    ser = B00.gf.series(n=32)
    for k in range(32):
        assert ser.coeff(z, k) == fibonacci(k+2)

def test_Catalan():
    '''
    Counting rooted plane binary trees gives the Catalan numbers
    '''
    Z = CombinatorialAtom()
    T = CombinatorialClass('T')
    ## a binary Tree is either empty or it is a node and two attached trees
    T = NeutralClass() + Z * T * T
    ser = T.gf.series(n=32)
    for k in range(32):
        assert ser.coeff(z, k) == catalan(k)

def test_dollar_change():
    '''
    This problem was orignally posed by P\'olya, according to Sedgewick &
    Flajolet. One wishes to determine the number of ways to make change for $1
    using coins of size $0.01, $0.05, $0.10, $0.25. Here, the order is not
    important, so we use a set construction.
    '''
    Z = CombinatorialAtom()
    P = MSET(SEQ(Z), 1)*MSET(SEQ(Z), 5) * MSET(SEQ(Z), 10) * MSET(SEQ(Z), 25)
    assert P.gf.series(z, 0, 101).coeff(z, 100) == 242
