----------------------------------------
Examples and tests for running Tantalus
----------------------------------------

The structure of this folder is as follows

 Examples/tests
   * scripts/   : contains all the individual examples/tests. Note that these
                  are written to function from the *main* directory i.e. 
                      
                          bash scripts/run.$WHATEVER.sh
  
                  works, but direct execution in scripts
                            
                          cd scripts/
                          bash run.$WHATEVER.sh

                  will not.

   * work/      : working directory
   * run.all.sh : execute all tests/examples

 Outputs:
   * out/       : contains the output of all example scripts 
   * reference/ : contains reference output of a previous version, 
                  verified to be correct
 Both of these are subdivided in 
   * STDOUT : output to STDOUT of the code  
   * summary: summary output, BXLFIT files


Individual tests are

  * run.minimal.sh
    - - - - - - - -
    Minimal HF calculation for 16O using SLy4 in a limited box.

    out/STDOUT/Tant.minimal.out

    out/summary/minimal.z008n008num001run001.out

  * run.pairing.sh
    - - - - - - - -

    Example calculation for 44Ca with (i) BCS and (ii) HFB pairing using 
    SLy5s1

    out/STDOUT/Tant.BCS.out
    out/STDOUT/Tant.HFB.out

    out/summary/BCS.z020n022.out
    out/summary/HFB.z020n022.out

  * run.constrained.sh
    - - - - - - - - - - 

    Example calculation for constraints for Ne20 using SLy5s1 in a limited box.
    Three calculations are done: Q20 = 10/20/30 fm^2.

    out/STDOUT/Tant.Q20=10.out
    out/STDOUT/Tant.Q20=20.out
    out/STDOUT/Tant.Q20=30.out

    out/summary/constrained.Q20=10.z010n010.out
    out/summary/constrained.Q20=20.z010n010.out
    out/summary/constrained.Q20=30.z010n010.out

  * run.continuing.sh
    - - - - - - - - - - -

    Example calculation demonstrating the initialisation of a calculation with 
    a wavefunction file. Two calculations are included: a starting one for 
    Ca40 using SLy5s1, followed by one for Ca42.

    out/STDOUT/Tant.Ca40.out
    out/STDOUT/Tant.Ca42.out

    out/summary/continuing.z020n020.out
    out/summary/continuing.z020n022.out

  * run.input.sh
    - - - - - - - - -

    Example calculation showcasing  
    (i) how to add extra spwfs to an existing calculation; and  
    (ii)how to transfer from a time-reversal conserving calculation to 
        a time-reversal breaking one

    out/STDOUT/Tant.O16.out
    out/STDOUT/Tant.O16.extra.out
    out/STDOUT/Tant.O16.T.out
    out/STDOUT/Tant.O16.T.extra.out

  * run.blocking.sh
    - - - - - - - - - -
 
    Example calculation showcasing the various option for blocking
      (1) Converges a false vacuum for Mg25
      (2) Uses the previous run to initialize an EFA calculation for Mg25
      (3) Uses that final run to perform a time-reversal broken calculation, 
         still blocking EFA-style
      (4) Finally, transfer the EFA-blocking to a real blocked calculation
    Note that all of these operate the DIRECT HFB solver, pairingscheme=0.

    out/STDOUT/Tant.Mg25.FV.out
    out/STDOUT/Tant.Mg25.EFA.out
    out/STDOUT/Tant.Mg25.EFA.T.out
    out/STDOUT/Tant.Mg25.block.out

  * run.gradient.sh
    - - - - - - - - - -
 
    Example calculation showcasing the difference between gradient and direct
    HFB solvers. 
       1) Converge a false vacuum for O19
       2) Block the lowest n+ state and try to converge the direct case
       3) Block the lowest n+ state and converge with the gradient solver

    out/STDOUT/Tant.O19.FV.out
    out/STDOUT/Tant.O19.direct.out
    out/STDOUT/Tant.O19.grad.out

  * run.BSkG1.sh
    - - - - - - - -

    Example calculations for the GSk1 and GSk2 parameterizations for Ca48.

    out/STDOUT/Tant.BSkG1.out
    out/STDOUT/Tant.BSkG2.out!

    out/summary/BSkG1.z020n028.out
    out/summary/BSkG2.z020n028.out

  * run.BSkG.T.sh
    - - - - - - - -

    Example calculations for the BSkG1 and BSkG2 parameterizations for Ca48, but
    this time with the time-reversal broken executable, verifying that it gives
    the right results when time-reversal is not explicitly broken.

    out/STDOUT/Tant.BSkG1.T.out
    out/STDOUT/Tant.BSkG2.T.out

    out/summary/BSkG1.T.z020n028.out
    out/summary/BSkG2.T.z020n028.out

  * run.magmoment.sh
    - - - - - - - - -
    This is an example script showing the calculation of the magnetic dipole 
    moment of the 7/2- ground-state in Sc41 with the SLy5s1 EDF.
    
    out/STDOUT/Tant.Sc41.SLy5s1.FV.out
    out/STDOUT/Tant.Sc41.SLy5s1.blocked.out

  * run.cranking.sh
    - - - - - - - - 

    Example cranking calculation using the SLy5s1 parameterisation for Ar36.  
    Note that the example is not so great, as we crank here around the 
    symmetry axis Z, such that we are not really talking about collective
    rotational motion.
        
    out/STDOUT/Tant.Ar36.out
    out/STDOUT/Tant.Ar36.om_z=0.1.out
    out/STDOUT/Tant.Ar36.om_z=0.2.out
    out/STDOUT/Tant.Ar36.om_z=0.3.out
    out/STDOUT/Tant.Ar36.om_z=0.4.out
    out/STDOUT/Tant.Ar36.om_z=0.5.out

  * run.Pb208.sh
    - - - - - - - -

    Example calculation using SLy4 for Pb208 in a more realistic 
    (16x16x16x0.8fm)   box.

    out/STDOUT/Tant.Pb208.out

    out/summary/Pb208.z082n126.out
    
  * run.N2LO.sh
    - - - - - - - -

    Example calculation using the N2LO parameterization SN2LO for Ca48.

    out/STDOUT/Tant.N2LO.out
    
    out/summary/N2LO.z020n028.out
    
  * run.P.sh
   - - - - - - - 
    Example demonstrating the breaking of parity symmetry and then adding some
    points on the mesh.
    
    out/STDOUT/Tant.BSkG1.symmetric.out
    out/STDOUT/Tant.BSkG1.Q30=30.out  
    out/STDOUT/Tant.BSkG1.Q30=30.added.out
