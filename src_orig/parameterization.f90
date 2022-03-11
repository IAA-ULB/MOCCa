module parameterization
 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______ 
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
 !                                                                       
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 
 use iso_fortran_env
 use compilation
 use geninfo
 use pairingcutoffs
 
 implicit none
    
    !---------------------------------------------------------------------------
    ! Automatically generated parameters of the parameterization
$PARAMDECL
    !---------------------------------------------------------------------------
    ! Constants that enter the game, together with the e2 from the Coulomb
    ! module.
    !---------------------------------------------------------------------------
    real(KIND=dp) :: hbm(2)         = 20.73551910_dp 
    real(KIND=dp) :: nucleonmass(2) = (/939.565379_dp , 938.272046_dp /)
    ! Hbar times the speed of light, in units of MeV fm
    ! Only used for calculating moments of inertia
    real(KIND=dp) :: hbarclum       = 197.32697_dp 
    !---------------------------------------------------------------------------  
    ! Treatment of the one-body COM correction.
    ! (0) not included
    ! (1) perturbatively included
    ! (2) completely included
    integer       :: COM1body = 2
    integer       :: COM2body = 0
    !---------------------------------------------------------------------------
    ! Rotational correction (see calcRotationalCorrection() in functional.f90 )
    ! (0) Not included
    ! (1) Included
    integer       :: RotCorr     = 0
    real(KIND=dp) :: RotcorrB    = 0, rotcorrC = 0
    ! Vibrational correction (see calcRotationalCorrection() in functional.f90 )
    real(KIND=dp) :: vibcorrB = 0, vibcorrL=0, vibcorrD = 0    
    ! Whether or not (rotcorr_cut) to use a cutoff in the definition of the 
    ! rotational correction.  
    logical       :: rotcorr_cut = .true.
    ! Parameters of the cutoff for the rotational correction
    ! These are initialised to -1,-1, but get set to the corresponding values
    ! of the pairing cutoff if not read from the .param file.
    real(KIND=dp) :: rotcutwindow(2) = -1.0d0 , rotcutmu(2) = -1.0d0
    !---------------------------------------------------------------------------
    ! Integer, whether or not coulomb is added
    integer       :: coultreatment = 1  
    ! Order of the finite difference discretisation of the laplacian for the 
    ! Coulomb solver
    integer :: CoulOrder = 2
    !---------------------------------------------------------------------------
    ! Stabilisation factors for the stabilised pairing scheme from 
    ! J. Erler et al., EPJA 37, 81-86 (2008).
    real*8 :: Estabp = 0.0, Estabn= 0.0
    !---------------------------------------------------------------------------
    ! Value of the electron charge, squared
    real(KIND=dp) :: e2 =1.43996446_dp 
    !---------------------------------------------------------------------------
    ! Effective size of the proton to take into account. 
    ! The proton density is folded with a Gaussian
    !   G(r = |r1 - r2|) = (r0 * sqrt(pi))**(-1) * exp[-(r/r0)**2]
    ! If r0 = 0, this is a delta-function and no folding is performed.
    real(KIND=dp) :: protonsize(2) = 0.0
    !---------------------------------------------------------------------------
    ! Parameters governing the charge form factor of the neutron.
    ! It is the difference of two Gaussians with different width
    ! See
    !   Chandra et al, PRC 13, 1976, 245
    !   
    ! G(r = |r1 - r2|) = (r_+ * sqrt(pi))**(-1) * exp[-(r/r_+)**2]
    !                  - (r_+ * sqrt(pi))**(-1) * exp[-(r/r_-)**2]
    !
    ! This array stores the r_+ and r_-  (in that order).
    !---------------------------------------------------------------------------
    real(KIND=dp) :: neutronsize(2) = 0.0
    !---------------------------------------------------------------------------
    ! Include (or not) the COM-correction due to the harmonic oscillator basis
    ! in both the neutron and proton form factors.
    logical :: HOCOMForm=.false.
    !---------------------------------------------------------------------------
    ! Whether to treat the inclusion of the nucleon form factors in the Coulomb
    ! module selfconsistently or not.
    logical       :: nucleonsize_selfconsistent=.true.
    !---------------------------------------------------------------------------
    ! Small non-zero value that can be used in a .func file to safeguard against
    ! division by zero in density dependent terms. 
    real(KIND=dp) :: eps =1d-20
    !===========================================================================
    ! ATTENTION: the following parameters are all intended to reenable 
    !            functionality that was found to be "suboptimal". Either these
    !            were true errors or things that did not work as desired. 
    !            Calculations with these flags enables should be handled with
    !            utmost care.
    ! 
    ! NeutronCoulombError:
    ! -------------------- 
    !   If you set this to .true., the charge form factor of the neutron will 
    !   not be taken into account correctly for single-particle hamiltonian. 
    !   This was the way the original GSk1 and GSk2 parameterizations were 
    !   fitted, and this is retained here for reproducing those calculations.
    logical :: neutroncoulomberror = .false.
    !===========================================================================
    
contains
    
   subroutine readparameterization(name_param, func_name) 
    !---------------------------------------------------------------------------
    ! Read the parameterization information from the .param file.
    
    character(len=20) :: name, func_file, toopen
    character(len=*), intent(in) :: name_param, func_name
    integer           :: io
    logical           :: exists
    
    ! predefined options
    namelist /skf/ name, func_file, hbm, e2, COM1body, COM2body, coultreatment,&
    &              protonsize, nucleonsize_selfconsistent, neutronsize,        &
    &              hocomform, rotcorr, rotcorrb, rotcorrc, coulorder,          &
    &              Estabp, Estabn, neutroncoulomberror, cutneutron, cutproton, &
    &              CutType,  rotcorr_cut, rotcutwindow, rotcutmu,              &
    &              vibcorrb, vibcorrl, vibcorrd, eps,                          &
    ! Those generated by Hephaestos
    $READPARAMS                        

    ! Trying some options for the .param file
    ! name_param.param
    ! NAME_PARAM.param (transfer everything to caps)
    ! name_param.param (transfer everything to lowercase)
    ! forces.param
    ! If none of these work, stop.
    toopen = trim(name_param)//'.param'
    inquire(exist=exists, file=toopen)
    if( .not. exists) then   
      toopen = trim(to_upper(name_param))//'.param'
      inquire(exist=exists, file=toopen)
      if(.not.exists) then
        toopen = trim(to_lower(name_param))//'.param'
        inquire(exist=exists, file=toopen)
        if(.not.exists) then
          toopen = trim('forces.param')  
          inquire(exist=exists, file=toopen)
          if(.not.exists) then
            print *, 'No parameterization file found.'
            print *, 'Valid filenames: - {param}.param (as-is)'
            print *, '                 - {PARAM}.param (all uppercase)'
            print *, '                 - {param}.param (all lowercase)'
            print *, '                 - forces.param'
            stop
          endif
        endif
      endif
    endif  
        
    open(unit=12, file=toopen,iostat=io)

    if (io.ne.0) then
      print *, "Problem opening the parameterization file!"
      print *, 'iostat = ', io
      stop
    endif

    do 
      read (unit=12, nml=skf, iostat=io)        ! Read the parameterization info
      if(io .ne. 0 .and. io.ne.iostat_end) then
        ! Something went wrong, but the file isn't at its end yet
        ! Maybe we read a parameter that doesn't exist for this EDF-type?
        ! Maybe something else; in any case, we proceed to the next &SKF/ 
        call resetparameterization()
        cycle
      endif
        
      name = to_upper(name)  ; func_file = to_upper(func_file)
      if(adjustl(to_upper(name_param)) .eq. adjustl(name)) then 
        ! Found the parameterization
        exit
      else
        if(io.eq.iostat_end) then 
          print *, 'Parameterization not found on the .param file.'
          stop
        else
          ! Reset all values
          call resetparameterization()
        endif
      endif
    enddo
    close(unit=12)
    !---------------------------------------------------------------------------
    ! Some sanity checks
    if(adjustl(func_file) .ne. adjustl(func_name)) then
        print *, '============================================================='
        print *, ' .func file used for compilation      = ', adjustl(func_name)
        print *, ' .func file for this parameterization = ', adjustl(func_file)
        print *, '============================================================='
        stop
    endif 
    
    if(adjustl(to_upper(name_param)) .ne. adjustl(name)) then
        print *, '============================================================='
        print *, ' .param asked for = ', adjustl(name_param)
        print *, ' .param read      = ', adjustl(name)
        print *, '============================================================='
        stop
    endif 
      
    !---------------------------------------------------------------------------
    ! Checking that all requested parameters have been read.
    ! Hephaestos generates a list of 'ifs' that check if the parameters are
    ! exactly their initializer
$CHECKPARAMS
    
    
    ! setting the rotcutwindow and rotcorrmu parameters, if not read from
    ! the .param file
    if(any(rotcutwindow .lt. 0.0d0)) then
      rotcutwindow(1) = cutneutron
      rotcutwindow(2) = cutproton
    endif
    if(any(rotcutmu .lt. 0.0d0)) then
      rotcutmu = pairingmu
    endif
    
  end subroutine readparameterization
    
  subroutine resetparameterization()
    !---------------------------------------------------------------------------
    ! Reset all values that can be read from a .param file to defaults.
    !---------------------------------------------------------------------------
    ! Reset physical constants 
    hbm            = 20.73551910_dp 
    nucleonmass(1) = 939.565379_dp 
    nucleonmass(2) = 938.272046_dp 
    hbarclum       = 197.32697_dp 
    e2             = 1.43996446_dp 
    ! COM treatment
    COM1body = 2      ;  COM2body = 0
    ! Coulomb
    coultreatment = 1 ; coulorder = 2
    neutroncoulomberror = .false.
    protonsize = 0.0d0 ; neutronsize= 0.0d0 
    nucleonsize_selfconsistent = .true.
    hocomform  = .false.
    ! Numerical safeguard
    eps        = 1d-20
    ! Rotational and vibrational corrections
    rotcorr    = 0
    rotcorrb   = 0 ; rotcorrc = 0
    vibcorrb   = 0 ; vibcorrl = 0 ; vibcorrd = 0
    rotcorr_cut  = .true.
    rotcutwindow = -1 
    rotcutmu     = -1.0d0
    ! Pairing options
    cutneutron = 5.0 ; cutproton  = 5.0
    cuttype    = 1
    PairingCut = 0   ; PairingMu =0.5
    Estabn = 0.0 ; Estabp = 0.0
    
    ! Reset paramets generated by Hephaestos
$RESETPARAMS 
  end subroutine resetparameterization
  
  subroutine printparameterization(param_name, func_name)
    !---------------------------------------------------------------------------
    ! Printing the available information on the parameterization/functional.
    !
    !---------------------------------------------------------------------------
    character(len=*), intent(in) :: param_name, func_name
    
    1 format (80('-'))
    2 format (' Parameterization name',  5x, 20a)
    3 format (' Functional name      ',  5x, 20a)
    4 format (' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
    5 format (' Parameters' )
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    7  format (' Coulomb Treatment = ', i2)
    70 format ('   No direct, nor exchange energy.')
    71 format ('   Direct and exchange included.')
    72 format ('   Only Direct contribution included.') 
  
    73 format ('   Discretisation of laplacian: ', i3, '-point stencil.')

    8  format ('   Protonsize  = ', 2f8.3, /,                                  &
    &          '   Neutronsize = ', 2f8.3 )
        
    81 format ('   Proton form factor  = One Gaussian')
    82 format ('   Proton form factor  = Difference of two Gaussians')
    83 format ('   Neutron form factor = One Gaussian')
    84 format ('   Neutron form factor = Difference of two Gaussians')
    
    85 format ('   Folding included self-consistently.')
    86 format ('   Folding not included self-consistently.')

    87 format ('   Nucleon radii corrected for HO basis size.')
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    9  format (' Center-of-mass options')
    93 format ('   One-body: self-consistent')
    92 format ('   One-body: perturbative')
    91 format ('   One-body: not included')
    96 format ('   Two-body: self-consistent')
    95 format ('   Two-body: perturbative')
    94 format ('   Two-body: not included')

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    103 format(' Rotational & vibrational correction ')
    104 format('    Not included.')
    105 format('    Included.    ')
   1051 format('     -> Attention: only calculated every ', i3, ' iterations.') 
    106 format('    Rot. Param. b =', f6.3)
    107 format('                c =', f6.3)
   1071 format('    Vib. Param. d =', f6.3)
   1072 format('                l =', f6.3)
   1073 format('                b =', f6.3)
    108 format('    Cutoff :     ACTIVE')
   1081 format('      window      =', 2f6.3, ' (MeV)')
   1082 format('      mu          =', 2f6.3, ' (MeV)')
    109 format('    Cutoff : NOT ACTIVE')
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    100 format( ' Elementary Constants')
    101 format(" e^2            = ", f11.8 ,' (sqrt(MeV fm)) ')
    102 format(" hbar^2/(2*mn)  = ", f11.8 ,' (MeV fm^2)'&
    &           ,/," hbar^2/(2*mp)  = ", f11.8 ,' (MeV fm^2)')
    200 format( ' Numerical details ')
    201 format( '   epsilon = ', es15.8) 

    999 format ('!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!',/&
    &           '!  Attention: neutroncoulomberror = .true.',                 /&
    &           '!  Neutron finite size not correctly propagated to sp-ham.', /&
    &           '!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!-!')  
    
    
    print 1
    print 2, adjustl(param_name)
    print 3, adjustl(func_name)
    print 4
    print 5
    
$PRINTPARAMS
    
    print 4
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    ! All coulomb options
    print 7, coultreatment
    select case(Coultreatment)
    case(0)
      print 70
    case(1)
      print 71
    case(2)
      print 72
    end select 

    print 73, 2*coulorder+1
    
    print 8, protonsize, neutronsize
    if(any(protonsize.ne.0.0) .or. any(neutronsize.ne.0.0)) then
        if(protonsize(2).ne.0.0) then
            print 82
        elseif(protonsize(1).ne.0.0) then
            print 81
        endif
        if(neutronsize(2).ne.0.0) then
            print 84
        elseif(neutronsize(1).ne.0.0) then
            print 83
        endif
      if(nucleonsize_selfconsistent) then
        print 85
      else
        print 86
      endif
      if(hocomform) then
          print 87
      endif
    endif

    if(neutroncoulomberror .and. nucleonsize_selfconsistent) then
      print 999
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    print 4
    print 9
    select case (COM1body)
    case(0)
      print 91
    case(1)
      print 92
    case(2)
      print 93
    end select
    select case (COM2body)
    case(0)
      print 94
    case(1)
      print 95
      print 1051, printiter
    case(2)
      print 96
      print *, 'Self-consistent inclusion of the two-body center of mass',  &
      &        ' correction is not available.'
      stop
    end select
  
    print 4
    print 103
    select case(rotcorr)
    case(0)
      print 104
    case(1)   
      print 105
      print 1051, printiter
      print 106, rotcorrb
      print 107, rotcorrc
      print 1071, vibcorrd
      print 1072, vibcorrl
      print 1073, vibcorrb
      
      if(rotcorr_cut) then
          print 108
          print 1081, rotcutwindow
          print 1082, rotcutmu
      else
          print 109
      endif      
    case DEFAULT
      print *, 'Illegal value of RotCorr = ', Rotcorr
      stop
    end select
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    ! Constants
    print 4
    print 100
    print 101, e2
    print 102, hbm
    
    print 200
    print 201, eps
  end subroutine printparameterization
  
  !=============================================================================
  ! Various functions that might be useful to define coupling constants in 
  ! the .func files.
  !=============================================================================
  
  real(KIND=dp) function Cc(t,x,p,s,it) result(c)
    !---------------------------------------------------------------------------
    ! Michaels convention for coupling coefficients of central Skyrme terms.
    ! 
    !
    ! C^+/-_cST       +              -
    !              t      tx      t    tx
    !       
    !    c00     +3/8     0     +5/8  +1/2 
    !    c01     -1/8   -1/4    +1/8  +1/4
    !    c10     -1/8   +1/4    +1/8  +1/4
    !    c11     -1/8     0     +1/8   0
    !---------------------------------------------------------------------------
    
    integer, intent(in)      :: p,s,it
    real(KIND=dp),intent(in) :: t,x
    
    c = 0
    
    select case(p)
    !---------------------------------------------------------------------------
    case(+1) ! C^+_c
      select case(s)          ! C^+_cS  
        case(0)               ! C^+_c0 
          select case(it)     ! C^+_c0T
          case(0)
            c = +3.0/8.0 * t                     ! C^+_c00
          case(1)
            c = -1.0/8.0 * t  - 1.0/4.0 * t * x  ! C^+_c01 
          end select
        case(1)               ! C^+_c1 
          select case(it)     ! C^+_c1T
          case(0)
            c = -1.0/8.0 * t  + 1.0/4.0 * t * x  ! C^+_c10
          case(+1)
            c = -1.0/8.0 * t                     ! C^+_c11
          end select
      end select
    !---------------------------------------------------------------------------
    case(-1) ! C^-_c
      select case(s)          ! C^-_cS  
        case(0)               ! C^-_c0 
          select case(it)     ! C^-_c0T
          case(0)
            c = +5.0/8.0 * t  + 1.0/2.0 * t * x  ! C^-_c00
          case(1)
            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  ! C^-_c01 
          end select
        case(1)               ! C^-_c1 
          select case(it)     ! C^-_c1T
          case(0)
            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  ! C^-_c10
          case(+1)
            c = +1.0/8.0 * t                     ! C^-_c11
          end select
      end select
    end select
    !---------------------------------------------------------------------------
  end function Cc


end module parameterization
