!-------------------------------------------------------------------------------
! Original code by E. Olsen, simplified for testing purposes by W.R.
!-------------------------------------------------------------------------------
! This code is for running Tantalus in parallel and contains an MPI for this purpose.
! It also contains the ability to fit forces, based of off Stephane Goriely's run_mocca.f and chi.f
! Also includes minsq.f for parameter fitting and minimization.
! Reads "datain" to get information on the interaction parameters and for which nuclei to make a mass table
! Each (Z,N) gets sent to a core to have its binding energy and associated ground state properties calculated;
! binding energy is sent back for chi.f to make further calculations

!-------------------------------------------------------------------------------
! 1. Master core determines necessary information for the mass table calculation
!  a. Force parameters (force.param)
!  b. Parameters for each (Z,N) (z000n000_in, where "000" is the respective Z or N)
! 2. Master core sends information to each worker core on which calculation to perform
! 3. All worker cores do their calculations and send their results back to the master core
! 4. Master core calculates \sigma_{rms}
! 5. If \sigma_{rms} is acceptable, the code stops. Elsewise, we start at step 1 again.
!-------------------------------------------------------------------------------
!====================================
!Main program: runs MOCCa with an MPI
!====================================

program run_mocca_parallel
  use iso_fortran_env
  use parallel_material
  use mpi_f08        

  implicit none

  integer(int32) :: my_id, num_workers, an_id, tag=1, ierr   !MPI integers
  integer(int32) :: worker_tag, worker_id
  integer(int32) :: i,j,k     !Iteration variables
  integer(int64) :: isgfit    !Type of calculation (with or without fit)
  integer(int64) :: num_nuc   !Number of nuclei to calculate
  integer(int64) :: file_number, counter, calc_status, file_counter
  integer(int64) :: nuc_left  !Counter for uncalculated nuclei
  integer(int64) :: mass_table_count  !Number of worker cores stopped

  character(11)  :: input_file
  character(11), allocatable :: input_data_files(:)

  integer(int64), allocatable :: isgz(:), isgn(:)

  call MPI_Init(ierr)                                 !Begins parallel execution
  call MPI_Comm_Rank(MPI_Comm_World,my_id,ierr)
  call MPI_Comm_Size(MPI_Comm_World,num_workers,ierr)

  counter=1  !Initialize the counter
  
  !===========
  !Master core
  !===========
  if (my_id == 0) then

    calc_status=0  !There's still work to do
    write(*,*) "Times run:", counter 

    !=======================================
    !Outer loop: all functions of the master
    !=======================================
    !do i=1,2
      mass_table_count=0   !All workers are available
    
      !----------------------------------------------------
      !Calculation not finished, read mass table input file
      !----------------------------------------------------
      if (calc_status == 0) then 
        call create_mass_table_input(counter,isgfit,num_nuc,isgz,isgn,input_data_files)
      end if

      nuc_left=num_nuc   !Total nuclei to calculate

      !=====================================================
      !Inner loop: master sends data to workers to calculate
      !=====================================================
    !do i=1,2
      calc_status=0
      nuc_left=num_nuc 
      do

        !-------------------------------------------------
        !No more nuclei to calculate: tell workers to stop
        !-------------------------------------------------
        if (nuc_left == 0) then
            
          call MPI_Recv(worker_id, 1, MPI_Integer, &
               MPI_Any_Source, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)
          
          calc_status=1   !This mass table calculation is over
          call MPI_Send(calc_status, 1, MPI_Integer, &
               worker_id, tag, MPI_Comm_World, ierr)

          mass_table_count=mass_table_count+1

          !--------------------------------------------
          !All workers have stopped, move to processing
          !--------------------------------------------
          if (mass_table_count >= (num_workers-1)) then
            exit
          else
            cycle
          end if
         
        end if
       
        !---------------------------------
        !Signal from available worker core
        !---------------------------------
        call MPI_Recv(worker_id, 1, MPI_Integer, &
             MPI_Any_Source, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)
      
        !--------------------------------------------------
        !Sends the status of the calculation to worker core
        !  -for calc_status=0, worker core runs something
        !  -for calc_status=1, worker core stops
        !--------------------------------------------------
        call MPI_Send(calc_status, 1, MPI_Integer, &
             worker_id, tag, MPI_Comm_World, ierr)

        !-----------------------------------------
        !Sends the z000n000_in file to this worker
        !-----------------------------------------
        input_file=input_data_files(nuc_left)      
        call MPI_Send(input_file, 11, MPI_Character, &
             worker_id, tag, MPI_Comm_World, ierr)

        nuc_left=nuc_left-1   !One less nucleus to calculate

        !--------------------------------------
        !Receives word when each worker is done
        !--------------------------------------
        !call MPI_Recv(calc_error, 1, MPI_Integer, &
        !     MPI_Any_Source, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)

        !write(*,*) "Core number", dest, "finished with input file", input_file

      end do

    !end do

    !---------------------------------------------
    !Checks to see if every nucleus was calculated
    !---------------------------------------------
    !j=1 
    !do
    !  inquire(file=output_data_files(j), exist=file_exists)

    !  if (file_exists) then
    !    j=j+1
    !    file_counter=file_counter+1

    !    if (file_counter == size(output_data_files)) then
    !      exit
    !    end if

    !  end if
    !end do

    write(*,*) "All the cores finished! Hooray!"
      
    !-----------------------------
    !Calculation is complete, exit 
    !-----------------------------
    !if (calc_status == 1) then
    !  write(*,*) "Everything is finished! Time to go!"
    !  exit
    !end if

    !======================
    !Post-processing begins
    !======================

    !---------------------------------------------
    !Consolidate all mass table data into one file
    !---------------------------------------------
    call mass_table_output(num_nuc,isgz,isgn)

    !-------------------------------------
    !No fit, just a mass table calculation
    !-------------------------------------
    if (isgfit == 0) then
    !if (counter == 4) then
      write(*,*) "Mass table, so we only run once"
      deallocate(input_data_files)
      deallocate(isgz,isgn)
      calc_status=1

      !-------------------------------------------------------------
      !Tells remainng workers to stop
      !  -Note: one worker will already be stopped in the upper loop
      !  (this is why we go to num_cores-2)
      !-------------------------------------------------------------
      !do i=1,num_cores-2
      !  call MPI_Recv(worker_id, 1, MPI_Integer, &
      !       MPI_Any_Source, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)

      !  call MPI_Send(calc_status, 1, MPI_Integer, &
      !       worker_id, tag, MPI_Comm_World, ierr)
      !end do

    end if

    !-------------------------------
    !Calculation of energy rms_error
    !-------------------------------
    !call rms_error_calculator

    !-----------------------------------
    !Goal: energy rms_error of < 0.8 MeV
    !-----------------------------------
    !if (rms_error < 0.8) then
    !   exit
    !else
    !  counter=counter+1
    !   call updated_mass_table_input()
    !end if

    counter=counter+1

      !exit   !Temporary

    !end do
 
  !===========
  !Worker core
  !===========
  else
    worker_id=my_id   !Worker is ready to calculate

    !do j=1,2
      do
        !---------------------------------------
        !Worker tells master core it's available
        !---------------------------------------
        call MPI_Send(worker_id, 1, MPI_Integer, &
             0, tag, MPI_Comm_World, ierr) 
      
        !-----------------------------------------------
        !Receives word from master to keep going or stop
        !-----------------------------------------------
        call MPI_Recv(calc_status, 1, MPI_Integer, &
             0, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)

        !-----------------------------------------
        !No more nuclei for this run, stop working
        !-----------------------------------------
        if (calc_status == 1) then
          write(*,*) "Core", my_id, "is finished!"
          exit

        !-----------------------------------
        !There are still nuclei to calculate
        !-----------------------------------  
        else
          write(*,*) "Core", my_id, "still has work to do!"  
        end if

        !-----------------------------------
        !Receive input file for Tantalus run
        !-----------------------------------
        call MPI_Recv(input_file, 11, MPI_Character, &
             0, tag, MPI_Comm_World, MPI_Status_Ignore, ierr)   !Master says "use this input file"

        !write(*,*) "Core", my_id, "will use file", input_file
        file_number=my_id+100
    
        call Tantalus(file_number, input_file)         !Runs the Tantulus code for this core
        !write(*,*) "Hello from processor", my_id, "which will use file", input_file

        !worker_tag=my_id

        !-------------------------------------------------
        !Calculation complete, sending info back to master
        !-------------------------------------------------
        !call MPI_Send(calc_error, 1, MPI_Integer, &
        !     0, worker_tag, MPI_Comm_World, ierr)
      
      end do
    !end do

  end if

  call MPI_Finalize(ierr)   !Ends the parallel execution

end program run_mocca_parallel
