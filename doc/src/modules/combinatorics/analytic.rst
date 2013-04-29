.. _combinatorics-analytic:

Analytic Combinatorics
======================

.. module:: sympy.combinatorics.analytic

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
upper-case letter for the values of the size funciton, which are
also the coefficients of the generating function. Square
brackets surrounding a power of `z` preceding a generating function,
`[z^N]C(z)` denotes "the coefficient of `z^N` in `C`".

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


Preview of Symbolic Transfer Theorems
-------------------------------------

Analytic combinatorics provides a formal language for defining combinatorial
structures, which then translates immediately to generating function equations.
As a preview, consider the following examples of increasingly complicated
structures that can be defined:

  1. Non-negative Integers, `\mathcal{I}`: one of the simplest constructions is
     the set of positive integers, `\{0, 1, 2, 3, ...\}` with the size function
     being the number of positive integers of magnitude `N`, of which there is
     obviously just `1`. In the language of analytic combinatorics, the
     positive integers are a sequence of atoms,
     `\mathcal{I} = \operatorname{S\small{EQ}}(\mathcal{Z})`. The symbolic
     transfer theorem immediately gives

     .. math ::
         I(z) = \frac{1}{1-z} = \sum_{N\geq0}{z^N}

    and the coefficient of `z^N` is `1` as expected.

  2. Finite-length binary strings, `\mathcal{B}`: the size function of interest
     here is the number of distinct (ordered) binary strings of length `N`.
     There are several ways to define binary strings, of which two follow:

        - A binary string is either empty, or it is a zero or a one followed
          by another binary string. The symbolic representation of this is
          `\mathcal{B} = \mathcal{E} + (\mathcal{Z}_0 + \mathcal{Z}_1) \times \mathcal{B}`
          which gives the generating function equation `B(z) = 1 + 2zB(z)`.
          Solving for `B(z)` gives

          .. math ::
              B(z) = \frac{1}{1-2z} = \sum_{N\geq0}{2^N z^N}

          so analytic combinatorics gives the answer that we knew: there are
          `2^N` binary strings of length `N`.

        - Another way to construct binary strings is as a sequence of zeros and
          ones: `\mathcal{B} = \operatorname{S\small{EQ}}(\mathcal{Z}_0 + \mathcal{Z}_1)`.
          The symbolic method gives the generating function equation

          .. math ::
              B(z) = \frac{1}{1-2z}

          as before, so we have shown that different but equivalent
          constructions give the same result in this case, verifying the
          consistency of the method.

  3. Binary trees, `\mathcal{T}`: As an example of a combinatorial question to
     which the answer is not immediately obvious, consider the problem of
     enumerating binary trees of a given size. Define a *binary tree* to be a
     tree where every node has zero or two children.  Additionally, there must
     be a distinguished root node. An *internal node* is a node with two
     children, and the counting function for the class is the number of
     internal nodes in the tree. A binary tree could be empty, as well.

     .. TODO: provide graphics

     Therefore, a binary tree can be empty, or it is a node with two binary
     trees attached.  The symbolic definition of binary trees is
     `\mathcal{T} = \mathcal{E} + \mathcal{Z} \times \mathcal{T} \times \mathcal{T}`.
     This gives the generating function equation

     .. math ::
         T(z) = 1 + zT(z)^2

     which we can solve with the quadratic equation (or SymPy):

     .. math ::
         T(z) = \frac{1-\sqrt{1-4z}}{2z}.

     With some algebra, one can show that the coefficients have an explicit
     formula, and that in fact they are equal to the famous Catalan numbers:

     .. math ::
         T_N = \frac{1}{N+1}\binom{2N}{N}

     but they can also be computed using SymPy:

     >>> from sympy import sqrt
     >>> from sympy.abc import z
     >>> T = (1-sqrt(1-4*z))/(2*z)
     >>> T.series()
                2      3       4       5    ⎛ 6⎞
     1 + z + 2⋅z  + 5⋅z  + 14⋅z  + 42⋅z  + O⎝z ⎠

Even this final example is a very basic and well-known combinatorial structure,
which may seem unimpressive to one familiar classical combinatorics.  The
examples above are intentionally simple to show the approach of analytic
combinatorics in familiar (or simple) situations.  The constructions of analytic
combinatorics can be used for far more than novel things of deriving results
about make familiar combinatorial structures and classes.  Variations on the
well known cases and combinations of the constructions listed here provide and
endless mix of new combinatorial structures, many of which would be difficult or
impossible to study without the formalism of analytic combinatorics.

Constructions
=============

This section provides an overview of the constructions that are available in
analytic combinatorics and demonstrates how to use them in SymPy.

Summary of Symbolic Transfer Theorems
-------------------------------------

Here is a summary of the symbolic transfer theorems between constructions and
their generating functions. Each of these constructions is explained in some
detail below, but the results are summarized here.

    +--------------------------------------------+------------------------------------------------------------------------------------+
    | Construction                               | Generating Function                                                                |
    +============================================+====================================================================================+
    | `\mathcal{A} + \mathcal{B}`                | `\displaystyle A(z) + B(z)`                                                        |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\mathcal{A} \times \mathcal{B}`           | `\displaystyle A(z) B(z)`                                                          |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{S\small{EQ}}(\mathcal{A})`  | `\displaystyle \frac{1}{1-A(z)}`                                                   |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{M\small{SET}}(\mathcal{A})` | `\displaystyle \exp\left( \sum_{k\geq 1}{\frac{A(z^k)}{k}} \right)`                |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{P\small{SET}}(\mathcal{A})` | `\displaystyle \exp\left( \sum_{k\geq 1}{\frac{(-1)^{k-1}}{k}} A(z^k)\right)`      |
    +--------------------------------------------+------------------------------------------------------------------------------------+
    | `\operatorname{C\small{YC}}(\mathcal{A})`  | `\displaystyle \sum_{k\geq 1}{\frac{\phi(k)}{k}}\log{\frac{1}{1-A(z^k)}}`          |
    +--------------------------------------------+------------------------------------------------------------------------------------+

Here, `\phi(k)` is Euler's totient function.

Basic Constructions
-------------------

The two basic combinatorial classes are the *neutral class*,
`\mathcal{E}`, and the *atomic class*, `\mathcal{Z}`, which have a
single element of size zero and one, respectively. By definition, their
generating functions are `E(z)=1` and `Z(z)=z`. From these, one can
build more elaborate structures using the operations of "disjoint sum",
"Cartesian product", and "sequence", among others.

.. autoclass:: CombinatorialClass
   :members:

.. autoclass:: CombinatorialAtom
   :members:

.. autoclass:: NeutralClass
   :members:

Disjoint Sum
------------

.. autoclass:: CombinatorialSum
   :members:

Product
-------

.. autoclass:: CombinatorialProduct
   :members:

Sequence
--------

.. autoclass:: SEQ
   :members:

Multiset
--------

.. autoclass:: MSET
   :members:

Power Set
---------

.. autoclass:: PSET
   :members:

Cycle
-----

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
