from sympy.combinatorics.analytic import (Multiton,
            CombinatorialAtom, CombinatorialSum, SEQ,
            CombinatorialProduct,
            Empty,
            )
from sympy.abc import z

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

def test_Seq():
    Z = CombinatorialAtom(0)
    seq = SEQ(Z)
    assert seq.gf == 1/(1-z)
