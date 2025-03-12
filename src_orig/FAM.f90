program FAM

 !==============================================================================
 ! ________ _______  _        _ _________ _______  _                 _______
 !(  _____/(  ___  )( (      ) |\__   __/(  ___  )( \      |\     /|(  ____ \
 !| (      | (   ) ||  \    /  |   ) (   | (   ) || (      | )   ( || (    \/
 !| |___   | (___) ||   \  /   |   | |   | (___) || |      | |   | || (_____
 !|  ___)  |  ___  || (\ \/ /) |   | |   |  ___  || |      | |   | |(_____  )
 !| |      | (   ) || | \  / | |   | |   | (   ) || |      | |   | |      ) |
 !| |      | )   ( || )  \/  ( |   | |   | )   ( || (____/\| (___) |/\____) |
 !(_/      |/     \||/        \)   )_(   |/     \|(_______/(_______)\_______)
 !
 !  Copyright W. Ryssens & P. Demol
 !
 !------------------------------------------------------------------------------
 ! A FAM-QRPA implementation to complement MOCCa.
 !==============================================================================

  use compilation
  use IO

  implicit none


  !-----------------------------------------------------------------------------
  ! Read input from STDIN
  ! 
  ! For FAMQRPA, the code should read in addition:
  ! 
  ! -  the type of perturbing operator/external field: E1, E2, M1, M2, ...
  !    and more complicated stuff when targetting beta-decay
  !    Important note: we will need to distinguish
  ! -  the frequency \omega of the perturbing field
  ! -  the 'size' of the perturbation to perform the finite differencing
  ! -  a smearing parameter to avoid discontinuities at the poles of the 
  !    response function
  ! 
  call ReadInput()




  print *, "Reached the end successfully" 

  ! end of one FAM calculation;

end program FAM
