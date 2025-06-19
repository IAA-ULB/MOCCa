This is a small fortran program to check whether mpi is ok.
[link](https://curc.readthedocs.io/en/latest/programming/MPI-Fortran.html)

Compile it as

    > ml calcua/2024a
    > ml ScaLAPACK
    > ml

    Currently Loaded Modules:
        1) calcua/2024a                   
        2) GCCcore/13.3.0                 
        3) zlib/1.3.1-GCCcore-13.3.0      
        4) binutils/2.42-GCCcore-13.3.0   
        5) GCC/13.3.0                    
        6) numactl/2.0.18-GCCcore-13.3.0       
        7) XZ/5.4.5-GCCcore-13.3.0             
        8) libxml2/2.12.7-GCCcore-13.3.0       
        9) libpciaccess/0.18.1-GCCcore-13.3.0  
        10) hwloc/2.10.0-GCCcore-13.3.0         
        11) OpenSSL/3                        
        12) libevent/2.1.12-GCCcore-13.3.0   
        13) UCX/1.16.0-GCCcore-13.3.0        
        14) libfabric/1.21.0-GCCcore-13.3.0  
        15) PMIx/5.0.2-GCCcore-13.3.0        
        16) PRRTE/3.0.5-GCCcore-13.3.0  
        17) UCC/1.3.0-GCCcore-13.3.0    
        18) OpenMPI/5.0.3-GCC-13.3.0
        19) gompi/2024a
        20) OpenBLAS/0.3.27-GCC-13.3.0
        21) FlexiBLAS/3.4.4-GCC-13.3.0
        22) ScaLAPACK/2.2.0-gompi-2024a-fb
    
    > mpif90 test.f90 -o test

Run it as:

    > mpirun -np 4 ./test.exe 
    Hello World from process:            0 of            4
    Hello World from process:            3 of            4
    Hello World from process:            2 of            4
    Hello World from process:            1 of            4

This works fine. 

When running Tantalus.BXL.exe we got 
*** An error occurred in MPI_Init
*** on a NULL communicator

ik voer script /data/antwerpen/201/vsc20170/tantalus_full/pasta-test/run-pasta-test.sh uit in folder /data/antwerpen/201/vsc20170/tantalus_full/pasta-test  en krijg volgende error An error occurred in MPI_Init on a NULL communicator
het zal wel weer iets doms van mij zijn maar ik heb voorlopig geen idee wat...
Hier is het begin van de output file /user/antwerpen/201/vsc20170/workspace/tantalus_full/pasta-test/slurm-2137249.out
 
    Running /user/antwerpen/201/vsc20170/.bashrc

    The following have been reloaded with a version change:
    1) calcua/all => calcua/2024a


    Currently Loaded Modules:
    1) calcua/2024a                        12) libevent/2.1.12-GCCcore-13.3.0
    2) GCCcore/13.3.0                      13) UCX/1.16.0-GCCcore-13.3.0
    3) zlib/1.3.1-GCCcore-13.3.0           14) libfabric/1.21.0-GCCcore-13.3.0
    4) binutils/2.42-GCCcore-13.3.0        15) PMIx/5.0.2-GCCcore-13.3.0
    5) GCC/13.3.0                          16) PRRTE/3.0.5-GCCcore-13.3.0
    6) numactl/2.0.18-GCCcore-13.3.0       17) UCC/1.3.0-GCCcore-13.3.0
    7) XZ/5.4.5-GCCcore-13.3.0             18) OpenMPI/5.0.3-GCC-13.3.0
    8) libxml2/2.12.7-GCCcore-13.3.0       19) gompi/2024a
    9) libpciaccess/0.18.1-GCCcore-13.3.0  20) OpenBLAS/0.3.27-GCC-13.3.0
    10) hwloc/2.10.0-GCCcore-13.3.0         21) FlexiBLAS/3.4.4-GCC-13.3.0
    11) OpenSSL/3                           22) ScaLAPACK/2.2.0-gompi-2024a-fb
    #---------------------------------------------------------------------------------
    Date:       =   Fri Jun 13 16:45:26 CEST 2025
    Host:       =   r5c05cn2.vaughan
    Directory:  =   /data/antwerpen/201/vsc20170/tantalus_full/pasta-test
    JOB_ID:     =   2137249
    #---------------------------------------------------------------------------------
    Running Tantalus.BXL.exe < data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3 > out_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3
    *** An error occurred in MPI_Init
    *** on a NULL communicator
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    ***    and MPI will try to terminate your MPI job as well)
    [r5c05cn3.vaughan:1732582] Local abort before MPI_INIT completed completed successfully, but am not able to aggregate error messages, and not able to guarantee that all other processes were killed!
    *** An error occurred in MPI_Init
    *** An error occurred in MPI_Init
    *** An error occurred in MPI_Init
    *** An error occurred in MPI_Init
    *** on a NULL communicator
    *** on a NULL communicator
    *** An error occurred in MPI_Init
    *** on a NULL communicator
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    ***    and MPI will try to terminate your MPI job as well)
    [r5c05cn3.vaughan:1732585] Local abort before MPI_INIT completed completed successfully, but am not able to aggregate error messages, and not able to guarantee that all other processes were killed!
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    *** on a NULL communicator
    *** An error occurred in MPI_Init
    *** on a NULL communicator
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    ***    and MPI will try to terminate your MPI job as well)
    [r5c06cn1.vaughan:1729295] Local abort before MPI_INIT completed completed successfully, but am not able to aggregate error messages, and not able to guarantee that all other processes were killed!
    ***    and MPI will try to terminate your MPI job as well)
    *** MPI_ERRORS_ARE_FATAL (processes in this communicator will now abort,
    ***    and MPI will try to terminate your MPI job as well)
    [r5c06cn2.vaughan:1705323] Local abort before MPI_INIT completed completed successfully, but am not able to aggregate error messages, and not able to guarantee that all other processes were killed!
    ...


This is a small fortran to check whether mpi is ok