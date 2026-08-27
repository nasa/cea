Theory
******

This section provides a minimal overview of the theoretical concepts underlying CEA.
The full mathematical background of the CEA can be found in RP 1311 Part 1 [1]_.
The goal here is to highlight the core assumptions and numerical strategy so that
users understand where the solver is reliable (and where it is not).

Minimization of Gibbs Energy
----------------------------

The first equation in the original reference document [1]_ is:

.. math::

   PV = n_{g} R T

This basic equation leads to a fundamental assumption of the code: condensed species
are assumed to occupy zero volume. There are no guardrails in the code to prevent you
from violating this assumption, so assess the realism of the results whenever large
amounts of condensed matter are expected.

With that out of the way: the equilibrium problem is based on the minimization of Gibbs free-energy:

.. math::

   g = \sum_{j}^{NS} \mu_{j} n_{j}

Importantly, there is a constraint on this minimization problem that the equilibrium product mixture
has to contain the same amount of each element that we started with:

.. math::

   \sum_{j}^{NS} a_{ij} n_{j} - b_{i}^{\circ} = 0, \qquad i = 1, \ldots, NE

where :math:`a_{ij}` is the stoichiometric coefficient matrix, i.e. a matrix that translates product species
amounts to element amounts based on the product species formulas. The product :math:`a_{ij} n_{j}` equals :math:`b_{i}`,
so the expression can be simplified if needed. For example, if the product species list
simply contained :math:`\mathrm{H}_{2}`, :math:`\mathrm{O}_{2}`, :math:`\mathrm{H}_{2}\mathrm{O}`, and :math:`\mathrm{OH}`,
then the matrix :math:`a_{ij}` would be:

.. math::

   \begin{array}{c|cc}
                & H & O \\
     \hline
     H_{2}      & 2 & 0 \\
     O_{2}      & 0 & 2 \\
     H_{2}O     & 2 & 1 \\
     OH         & 1 & 1
   \end{array}

Now we have a minimization problem with a constraint, so we can define a Lagrangian in the form of:

.. math::

   \mathcal{L} = g + \sum_{i}^{NS} \lambda_{i} \left( b_{i} - b_{i}^{\circ} \right)

And then we find the stationary point of this Lagrangian to solve the constrained equilibrium problem:

.. math::

   \delta \mathcal{L} = 0

This results in a nonlinear system of equations, and a Newton method with a damped update is used to converge the solution.

The solution variables in the problem are:

* Species concentrations: :math:`\mathbf{n} = [n_{1}, ..., n_{j}]^{T}`
* Total mixture moles: :math:`n`
* Mixture temperature: :math:`T`

An initial guess is required for the iterative solution procedure:

* :math:`\mathbf{n} = [\frac{0.1}{NG}, ..., \frac{0.1}{NG}]^{T}` kg-mol/kg
* :math:`n = 0.1` kg-mol/kg
* :math:`T=3800` K

Using these initial conditions, the residual system of equations is solved iteratively and the solution
variables are updated using a damped update procedure to improve stability. The matrix system
of equations computes updates to four variables:

1. Change in log of gas species concentrations: :math:`\Delta \ln (n_{j}), j = 1, ..., NG`
2. Change in condensed species concentration: :math:`\Delta n_{j}, j = NG+1, ..., NS`
3. Change in log-moles: :math:`\Delta \ln (n)`
4. Change in log-temperature: :math:`\Delta \ln (T)`

This describes the general case, but some update variables are not included for certain problem types.
The damped update factor :math:`\lambda` (not to be confused with the Lagrange multiplier) is computed
using heuristics not described here, but typically :math:`\lambda` is less than 1 when the
solution is far from convergence and equal to 1 when the solution is close to convergence. This prevents
overshooting the solution initially and switches to the calculated step size near convergence when
the linear approximation is more accurate. For example, the updated value of log-concentrations is computed as:

.. math::

    \ln (n_{j})^{k+1} = \ln (n_{j})^{k} + \lambda (\Delta \ln (n_{j})), j = 1, ..., NG

.. note::

   Internally, :math:`\ln (n_{j})` and :math:`n_{j}` are stored separately. This seems redundant, but there
   is an important reason for this: a threshold is applied for species with :math:`\ln (n_{j}) < \ln (10^{-8})`,
   setting :math:`n_{j}=0` for those species. However, even after a species has been zeroed out, :math:`\ln (n_{j})`
   will still be continually updated, allowing the species to come in and out of being included in the mixture
   without loss of information from previous iterations. The default threshold value of :math:`\ln (10^{-8})=-18.420681`
   can be updated in the code by modifying the ``trace`` parameter.

Convergence Criteria
^^^^^^^^^^^^^^^^^^^^

After each iteration, a series of convergence checks are evaluated. Convergence is based on the relative size of
the update variables; if the updates are sufficiently small for all variables, we assume that
an equilibrium point has been reached. There are tests for each of:

* Gas species concentrations: :math:`n_{j}, j=1, ..., NG`
* Condensed species concentrations: :math:`n_{j}, j=NG+1, ..., NS`
* Total moles: :math:`n`
* Temperature: :math:`T`
* Element balance: :math:`b_{i}-b_{i}^{\circ}`
* Entropy: :math:`s`
* Modified Lagrange multipliers: :math:`\pi_{i}, i=1,...,NE`

For positive element inventories above :math:`10^{-6}` kmol/kg, the element-balance
residual is scaled by that element's own inventory, not the largest inventory.
The sum of absolute gas-species updates carrying that element must also be small
relative to its inventory. This prevents premature convergence of minor elements.

Equilibrium sensitivities use the species-activation threshold saved at convergence.
Final reporting can expose smaller species without another equilibrium solve; its
threshold must not replace the one used by the differentiated residual. With hard
activation, derivatives describe a fixed activity branch and need not agree with
finite differences that cross an activation boundary.

Next we will address the treatment for a number of special cases in the equilibrium problem:

Condensed Species
-----------------

First and foremost: **the inclusion of condensed species violates the ideal gas assumption.**
As long as the included condensed species take up approximately zero volume this is a fair assumption,
but the degree to which it is violated will affect the results; the user must gauge the error introduced.
With that disclaimer out of the way, the solution procedure for adding and removing condensed species is:

1. Start the initial guess with no condensed species included in the product mixture.
2. Compute equilibrium with only gas-phase species.
3. After convergence, test if adding any condensed species lowers the Gibbs energy of the mixture.
4. If any meet the criteria, add the species that lowers the Gibbs energy the most, and converge the system again.
5. Repeat steps 3 and 4 until final convergence.

There are also additional considerations for handling phase changes or the inclusion of multiple phases
of the same species.

*To be continued.*

.. Ionized Species
.. ---------------

.. Condensed Species
.. -----------------

.. Negative Reactant Amounts
.. -------------------------

.. Inert Species
.. -------------

.. Thermodynamic Properties
.. ------------------------

.. Transport Properties
.. --------------------

.. Frozen versus Equilibrium Transport Properties
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. Rocket Analysis
.. ---------------

.. Shock Problems
.. --------------

.. Detonation Problems
.. -------------------

.. Analytic Total Derivatives
.. --------------------------

.. [1] Gordon, S., McBride, B.J., "Computer program for calculation of complex chemical equilibrium compositions and applications. Part 1: Analysis",
    NASA RP-1311, 1994. [NTRS](https://ntrs.nasa.gov/citations/19950013764)
