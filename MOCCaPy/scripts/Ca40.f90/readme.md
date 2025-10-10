tantalus::run_tantalus(input_file='tant.data')
    io::ReadInput(input_file)
        geninfo::read_geninfo()
            reads namelists nucleus and mesh
            call inimesh
    readfunctional(file_number)
    initpairing(file_number)
    ReadEvolution(file_number)
    ReadSCFIteration(file_number)
    ReadWFdata(file_number)
    ReadIOInput(file_number)
    
    if(N_inertia .gt. 4) then
      read_inertia(file_number=file_number)
    else 
      inertia_l = inertia_l_hardcoded
      inertia_m = inertia_m_hardcoded
    endif
    call readmomentdata(file_number)
    call readcranking(file_number)

#if(0 == 1)
    call readfam(file_number)
#endif

    if(present(file_number)) then
      close(unit=file_number)
    endif
