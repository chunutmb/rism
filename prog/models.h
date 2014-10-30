#ifndef MODELS_H__
#define MODELS_H__



/* built-in models from literature
 *
 * References:
 *
 * Solution of a new integral equation for pair correlation function in molecular liquids
 * Lawrence J. Lowden and David Chandler
 * J. Chem. Phys. 59 (12) 6587 (1973)
 *
 * Applications of the RISM equation to diatomic fluids:
 * the liquids nitrogen, oxygen and bromine
 * C.S. Hsu, David Chandler and L.J. Lowden
 * Chem. Phys. 14 213-228 (1976)
 *
 * Comparisons of Monte Carlo and RISM calculations of pair correlation functions
 * David Chandler, C. S. Hsu and William B. Streett
 * J. Chem. Phys. 66 (11) 5231 (1977)
 *
 * Computation of the correlation functions for fluids composed of
 * diatomic molecules by means of the method of integration equations
 * Kazumitsu Kojima and Kiyoshi Arakawa
 * Bulletin of the Chemical Society of Japan 51(7) 1977-1981 (1978)
 *
 * An extended RISM equation for molecular polar fluids
 * Fumio Hirata and Peter J. Rossky
 * Chem. Phys. Lett. 83(2) 329-334 (1981)
 *
 * Application of an extended RISM equation to dipolar and quadrupolar fluids
 * Fumio Hirata, B. Montgomery Pettitt, Peter J. Rossky
 * J. Chem. Phys. 77(1) 509-520 (1982)
 *
 * Integral equation prediction of liquid state structure for
 * waterlike intermolecular potentials
 * B. Montgomery Pettitt and Peter J. Rossky
 * J. Chem. Phys. 77(3) 1451-1457 (1982)
 *
 * The interionic potential of mean force in a molecular polar
 * solvent from an extended RISM equation
 * Fumio Hirata, Peter J. Rossky, and B. Montgomery Pettitt
 * J. Chem. Phys. 78(6) 4133-4144 (1983) */
model_t models[] =
{
  {0}, /* empty model, place holder */
  /* 1. LC 1973, and Model I of CHS 1977 */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.500, 0.500}, {0.600}, 1.000, HARD_SPHERE,
    {}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    0, 0.5
  },
  /* 2. Model II of CHS 1977 */
  {2, {0.790, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.686, 0.686}, {0.490}, 1.000, HARD_SPHERE,
    {}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    0, 0.5
  },
  /* 3. Model III of CHS 1977 */
  {2, {0.675, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.825, 0.825}, {0.346}, 1.000, HARD_SPHERE,
    {}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    0, 0.5
  },
  /* 4. LC1973, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.696, 0.696}, {1.1/3.341}, 1/1.83, LJ_REPULSIVE,
    {}, 0, 0,
    IE_PY, 20.48, 1024,
    1, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    0, 0.5
  },
  /* 5. KA1978, liquid nitrogen */
  {2, {1.000, 1.000}, {{1, 1}, {1, 1}}, {{0}},
    {0.6964, 0.6964}, {1.1/3.341}, 1/1.61, LJ_FULL,
    {}, 0, 0,
    IE_PY, 20.48, 1024,
    5, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    20, 0.5
  },
  /* 6. HR1981, liquid nitrogen, neutral */
  {2, {3.341, 3.341}, {{1, 1}, {1, 1}}, {{0}},
    {0.01867, 0.01867}, {1.1}, 1/1.636, LJ_FULL,
    {}, 0, 0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.2
  },
  /* 7. HR1981, liquid nitrogen, charged, also HPR1982, model I */
  {2, {3.341, 3.341}, {{44.0, 44.0}, {44.0, 44.0}}, {{0}},
    {0.01867, 0.01867}, {1.1}, 1./72, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.2
  },
  /* 8. HPR1982, HCl, model II */
  {2, {2.735, 3.353}, {{20.0, 20.0}, {259.0, 259.0}}, {{0}},
    {0.018, 0.018}, {1.257}, 1./210, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.2
  },
  /* 9. HPR1982, HCl, model III */
  {2, {0.4, 3.353}, {{20.0, 20.0}, {259.0, 259.0}}, {{0}},
    {0.018, 0.018}, {1.3}, 1./210, LJ_FULL,
    {0.200, -0.200}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.2
  },
  /* 10. PR1982, H2O, model I
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {2.8, 0.4, 0.4}, {{0}},
    { {262.566, 309408} /* O-O */,
      {0, 689.348} /* O-H1 */, {0, 689.348} /* O-H2 */,
      {0, 610.455} /* H1-H1 */, {0, 610.455} /* H1-H2 */, {0, 610.455} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696}, 1./(KBNAC*300), LJ_FULL,
    {-0.866088, 0.433044, 0.433044}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.3
  },
  /* 11. PR1982, H2O, model II (SPC)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.633081624 (109.48 degree) */
  {3, {3.166, 0.4, 0.4}, {{0}},
    { {-625.731, 629624} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */, {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {1.0, 1.0, 1.633081624}, 1./(KBNAC*300), LJ_FULL,
    {-0.82, 0.41, 0.41}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.3
  },
  /* 12. PR1982, H2O, model III (TIPS)
   * atom 0: O, atom 1: H1, atom 2: H2
   * C6/C12 are used instead of sigma/epsilon
   * d(H1, H2) = 1.51369612 (104.5 degree) */
  {3, {3.215, 0.4, 0.4}, {{0}},
    { {-525.000, 580000} /* O-O */,
      {0, 225.180} /* O-H1 */, {0, 225.180} /* O-H2 */,
      {0, 0} /* H1-H1 */, {0, 0} /* H1-H2 */, {0, 0} /* H2-H2 */ },
    {0.03334, 0.03334, 0.03334},
    {0.9572, 0.9572, 1.513696}, 1./(KBNAC*300), LJ_FULL,
    {-0.8, 0.4, 0.4}, KE2C, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.3
  },
  /* 13. SPCE, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/spce */
  {3, {3.1666, 0.4, 0.4},
    {{78.2083543, 78.2083543}, {0, 23.150478}, {0, 23.150478}}, {{0}},
    {0.033314, 0.033314, 0.033314},
    {1.0, 1.0, 1.633}, 1./300, LJ_FULL,
    {-0.8476, 0.4238, 0.4238}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    25, 0.3
  },
  /* 14. TIP3, H2O
   * atom 0: O, atom 1: H1, atom 2: H2
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/tip3 */
  {3, {3.15, 0.4, 0.4},
    { {76.5364, 76.5364}, {0, 23.1509}, {0, 23.1509} }, {{0}},
    {0.033314, 0.033314, 0.033314},
    {0.95719835, 0.95719835, 1.5139}, 1./300, LJ_FULL,
    {-0.834, 0.417, 0.417}, KE2PK, 1.0,
    IE_HNC, 40.96, 2048,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    25, 0.3
  },
  /* 15. HRP1983, solvent + solute (+/- ion pair)
   * Cf. model 7 */
  {4, {3.341, 3.341, 3.341, 3.341},
    { {44.0, 44.0}, {44.0, 44.0}, {44.0, 44.0}, {44.0, 44.0} }, {{0}},
    {0.01867, 0.01867, 0, 0},
    {1.1}, 1./200, LJ_FULL,
    {0.200, -0.200, 1.0, -1.0}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    12, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.2
  },
  /* 16. SPCE, H2O, Na+, Cl-
   * atom 0: O, atom 1: H1, atom 2: H2, atom 3: Na, atom 4: Cl-
   * the following data are copied from /Bossman/Software/3Drism/h2o_lib/spce */
  {5, {3.1666, 0.4, 0.4, 2.35, 4.4},
    { {78.2083543, 78.2083543}, {0, 23.150478}, {0, 23.150478},
      {65.42, 65.42}, {50.32, 50.32} }, {{0}},
    {0.033314, 0.033314, 0.033314, 0, 0},
    {1.0, 1.0, 0, 0, 1.633}, 1./300, LJ_FULL,
    {-0.8476, 0.4238, 0.4238, 1, -1}, KE2PK, 1.0,
    IE_HNC, 20.48, 1024,
    10, 100000, 1e-7,
    SOLVER_LMV,
    0.01,
    10, 0.3
  },
};



#endif /* MODELS_H__ */

