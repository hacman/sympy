.. _combinatorics-analytic:

Analytic Combinatorics
======================

Analytic combinatorics is an approach to studying combinatorial
structures where the central object of study is the generating function.
A *combinatorial class* is a set, `\mathcal{C}`, together with a
function `|\cdot|:\mathcal{C} \to \mathbb{Z}_{\geq 0}` mapping items in
`\mathcal{C}` to non-negative integers. This is called the *size function*
or simply the size. Combinatorics generally concerns itself with counting
all of the elements of a class of a given size, `N`, which is denoted
`C_N`. Define the generating function as the formal power series

.. math ::
    C(z) = \sum_{N \geq 0}{C_N z^N}.

Notational Aside: The convention is to use a script upper-case letter
for the class, a lower-case letter for an element of the class, a
regular upper-case letter for the generating function, and a subscripted
upper-case letter for the coefficients of the generating function. Square
brackets surrounding a power of `z` preceding a generating function,
`[z^N]C(z)`  denotes "the coefficient of `z^n` in `C`.

The basic identity that allows analytic combinatorics to work is

.. math ::
    \sum_{N \geq 0}{C_N z^N} = \sum_{c \in \mathcal{C}}{z^{|c|}}.

That is, the generating function can be recovered by summing over all
elements of the class, raising `z` to the size of the element. A moment's
reflection should convince you that this will work.

Analytic combinatorics has two collections of powerful *transfer
theorems*, mapping definitions of combinatorial structures to the
appropriate generating function and then allowing extraction of asymptotic
representations of the coefficients of the generating function. Full
explanation and proofs can be found in Flajolet & Sedgewick
[Flajolet2009]_.

Alternately, we can forego the second class of transfer theorems and
use SymPy to get the coefficients of the series representation of the
generating function.

The two basic combinatorial classes are the *neutral class*,
`\mathcal{E}`, and the *atomic class*, `\mathcal{Z}`, which have a
single element of size zero and one, respectively. By definition, their
generating functions are `E(z)=1` and `Z(z)=z`. From these, one can
build more elaborate structures using the operations of "disjoint sum",
"Cartesian product", and "sequence", among others.

Constructions
-------------

Disjoint Union
~~~~~~~~~~~~~~

Product
~~~~~~~

Sequence
~~~~~~~~

`\operatorname{S\small{EQ}}(\mathcal{A})`

Multiset
~~~~~~~~

Power Set
~~~~~~~~~

Cycle
~~~~~

Summary of Symbolic Transfer Theorems
-------------------------------------

    +--------------------------------------------+------------------------------------------------------------------------------------+
    | Construction                               | Generating Function                                                                |
    +============================================+====================================================================================+
    | `\mathcal{A} + \mathcal{B}`                | `\displaystyle A(z) + B(z)`                                                        |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\mathcal{A} \times \mathcal{B}`           | `\displaystyle A(z) B(z)`                                                          |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{S\small{EQ}}(\mathcal{A})`  | `\displaystyle \frac{1}{1-A(z)}`                                                   |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{M\small{SET}}(\mathcal{A})` | `\displaystyle \exp\left( \sum_{k\geq 1}{\frac{z^k}{k}} \right)`                   |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{P\small{SET}}(\mathcal{A})` | `\displaystyle \exp\left( \sum_{k\geq 1}{\frac{(-1)^{k-1}}{k}} A(z^k)\right)`      |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{C\small{YC}}(\mathcal{A})`  | `\displaystyle \sum_{k\geq 1}{\frac{\phi(k)}{k}}\log{\frac{1}{1-A(z^k)}}`          |
    +--------------------------------------------+------------------------------------------------------------------------------------+

Here, `\phi(k)` is Euler's totient function.

API
---

.. module:: sympy.combinatorics.analytic

.. autoclass:: CombinatorialClass
   :members:

.. autoclass:: CombinatorialAtom
   :members:

.. autoclass:: NeutralClass
   :members:

.. autoclass:: CombinatorialSum
   :members:

.. autoclass:: CombinatorialProduct
   :members:

.. autoclass:: SEQ
   :members:

.. autoclass:: MSET
   :members:

.. autoclass:: PSET
   :members:

.. autoclass:: CYC
   :members:

References
----------

.. [Flajolet2009] Flajolet, Philippe and Robert Sedgewick. Analytic
        Combinatorics.  New York: Cambridge University Press, 2009. Print.
        Online version: <http://algo.inria.fr/flajolet/Publications/book.pdf>

.. [Flajolet2013] Flajolet, Philippe and Robert Sedgewick. An Introduction to
        the Analysis of Algorithms, 2nd Edition.  Westford, Massachusetts:
        Addison-Wesley. 2013.
