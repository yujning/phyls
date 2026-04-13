STP Computation
===============

The basic Semi-Tensor Product (STP) Computation based on Eigen Library 

Basic STP Computation
----------------------

Native Definition
^^^^^^^^^^^^^^^^^^^^^

Given two matrices :math:`A_{m \times n}` and :math:`B_{p \times q}`, the STP
of :math:`A` and :math:`B` is defined as 

.. math::

  A \ltimes B = (A \bigotimes I_{t/n}) \cdot (B \bigotimes I_{t/p})

where :math:`\bigotimes` is Kronecker product, :math:`I` is identity
matrix, and :math:`t` is least common multiple (LCM) of :math:`n` and :math:`p`.
As an example, suppose 

.. math::

  A = \begin{bmatrix}
  1 & 0 & 0 & 0 \\
  0 & 1 & 1 & 1
  \end{bmatrix}

and

.. math::

  B = \begin{bmatrix}
  1 & 1 & 0 & 1 \\
  0 & 0 & 1 & 0
  \end{bmatrix},

then :math:`m=2, n=4` and :math:`p=2, q=4`, the LCM :math:`t` is thus 4.
According to the definition, 

.. math::

  \begin{align}
  A \ltimes B = (A \bigotimes I_{4/4}) \cdot (B \bigotimes I_{4/2}) &=
  \begin{bmatrix}
  1 & 0 & 0 & 0 \\
  0 & 1 & 1 & 1
  \end{bmatrix} \cdot
  \begin{bmatrix}
  1 & 0 & 1 & 0 & 0 & 0 & 1 & 0 \\
  0 & 1 & 0 & 1 & 0 & 0 & 0 & 1 \\
  0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 
  \end{bmatrix} \notag \\ 
  &= 
  \begin{bmatrix}
  1 & 0 & 1 & 0 & 0 & 0 & 1 & 0 \\
  0 & 1 & 0 & 1 & 1 & 1 & 0 & 1
  \end{bmatrix}. \notag
  \end{align}

Copy Definition
^^^^^^^^^^^^^^^^^^^^^
Continue with the above example, we can view :math:`A` be composed by two
submatrices :math:`A_{left}=\begin{bmatrix}1 & 0 \\ 0 & 1\end{bmatrix}` 
and :math:`A_{right} = \begin{bmatrix}0 & 0 \\ 1 & 1\end{bmatrix}`. When we
compute the STP :math:`A \bigotimes B`, we can examine :math:`B` column by
column from leftside to rightside, if the column in :math:`B` is
:math:`\begin{bmatrix} 1 \\ 0 \end{bmatrix}`, we copy :math:`A_{left}` as a
partial result; otherwise, we copy :math:`A_{right}`. One can verify the
results are exactly the same as the ones computed by native definition.

Functions
^^^^^^^^^^^^^^^^^^^^^
In header file ``stp/stp_eigen.hpp``, we provide function::

  matrix semi_tensor_product( const matrix& A, const matrix& B 
                              const bool verbose = false,
                              const stp_method method = stp_method::copy_method )

to compute the STP of matrices :math:`A` and :math:`B`, where toggle ``verbose`` is off and toggle ``stp_method``
is used by the copy definition by default.

Example

.. code-block:: c++
  
  matrix A;
  matrix B;
  
  //default
  auto result = stp::semi_tensor_product( A, B );

  //print verbose information
  auto result = stp::semi_tensor_product( A, B, true );
  
  //use native definition for STP computation
  auto result = stp::semi_tensor_product( A, B, true, stp::stp_method::native_method );

One can find more examples or test cases in ``examples/stp_eigen.cpp`` and ``test/stp_eigen.cpp``.

Matrix Chain STP Computation
----------------------------
When we have :math:`n` matrices multiplication and :math:`n \ge 3`, we call
this as matrix chain STP computation. 

Sequence
^^^^^^^^^^^^^^^^^^^^^
The matrix chain STP is calculated by sequential order. For example, we have 4
matrices :math:`A`, :math:`B`, :math:`C`, and :math:`D`. The parenthesized
expression of the matrix chain multiplication is as follows: 

.. math::
  ABCD = (((AB)C)D).

Dynamic Programming
^^^^^^^^^^^^^^^^^^^^^
Since the computational complexity varies depending on the method of
parenthesizing, we also propose a dynamic programming approach for matrix
chain STP computation. This approach allows us to determine the optimal
parenthesization of the matrix chain.

.. math::
  ABCD = ((AB)(CD)).

Multi-threads
^^^^^^^^^^^^^^^^^^^^^
Once we obtained the computation orders based on dynamic programming, the
computation can also invoke multi-threads to accerlerate.

Functions
^^^^^^^^^^^^^^^^^^^^^
In header file ``stp/stp_eigen.hpp``, we provide function::

  matrix matrix_chain_multiply( const matrix_chain& mc, 
                                const bool verbose = false,
                                const mc_multiply_method method = mc_multiply_method::dynamic_programming )

to compute the STP of matrix chain :math:`mc`, where toggle ``verbose`` is off and toggle ``mc_multiply_method``
is used by the dynamic programming by default.

Example

.. code-block:: c++
  
  matrix_chain mc;
  
  //default
  auto result = stp::matrix_chain_multiply( mc );

  //print verbose information
  auto result = stp::matrix_chain_multiply( mc, true );
  
  //use sequence method for matrix chain STP computation
  auto result = stp::matrix_chain_multiply( mc, false, mc_multiply_method::sequence );

One can find more test cases in ``test/stp_eigen.cpp``.
