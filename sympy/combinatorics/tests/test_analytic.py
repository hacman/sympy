from sympy.combinatorics.analytic import (Multiton,
            CombinatorialClass,
            CombinatorialAtom, CombinatorialSum, SEQ,
            CombinatorialProduct, CYC,
            NeutralClass,
            )
from sympy.abc import z
from sympy import catalan

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

def test_Catalan():
    '''
    Counting rooted plane binary trees gives the Catalan numbers
    '''
    Z = CombinatorialAtom()
    T = CombinatorialClass('T')
    ## a binary Tree is either empty or it is a node and two attached trees
    T = NeutralClass() + Z * T * T
    ser = T.gf.series(n=25)
    for k in range(25):
        assert ser.coeff(z, k) == catalan(k)

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
