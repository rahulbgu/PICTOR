program MAIN
     use parameters 
     use vars
     use initialise
#ifdef cyl
     use cyl_vars
	 use cyl_initialise 
#endif 	 
     use movdep
     use particles
	 use em_update
     use fields
     use savedata
     use comm_fldprtl
     use memory
	 use reload
	 use loadbalance	 
#ifdef cyl 
     use cyl_bc 
#else
     use bc 
#endif   	 
#ifdef gpu 
	 use var_gpu
	 use initialise_gpu
	 use particles_gpu 
	 use comm_fld_gpu
	 use em_update_gpu
	 use movdep_gpu 
	 use fields_gpu
#endif  
     use setup

     implicit none

call CheckRestart
call StartMPI 
if(restart) then 
	call RestartSteps
else  
    call Startups
endif 
         call StartTimer(1)
		 
do t=tstart,tfinish
          call StartTimer(31) !total time taken is stored in 31
		  call StartTimer(32) !this execution time is used in balancing the load
	 	  
!----------------------------------------------
! Update magnetic field by Half a time step 
!----------------------------------------------           
		   call StartTimer(2)
#ifdef CPU		   
     call UpdateBfieldHalf
#endif	 
#ifdef gpu
     call UpdateBfieldHalfGPU
#endif
           call StopTimer(2)
		   
!----------------------------------------------
! Mover particles and test particles and deposit current 
!----------------------------------------------   		   
           call StartTimer(4)
#ifdef CPU 
	     call MoveDepositPrtl
	  	 call MoveTestPrtl
#endif		 		   
#ifdef gpu 
         call MoveDepositPrtlGPU
         call MoveTestPrtlGPU	
#endif 	 

           call StopTimer(4)
		   
!----------------------------------------------
! Some operations performed right after moving the particles 
!----------------------------------------------   		   
!GPU version of BC routines should be defiend and called from these subroutines if needed 		  
     call EnforceBC_PostMovDep 
	 !This subroutine is defined in setup_*.F90
     call PostMovDep ! It allows for some setup dependent operations right after moving the particles and depositing their current   
 !----------------------------------------------
 !Update magnetic field by another half step and electric field by full step 
 !----------------------------------------------  	  
           call StartTimer(5)
#ifdef CPU		   
     call UpdateBfieldHalf 
	 call UpdateEfield
#endif	 
#ifdef gpu
     call UpdateBfieldHalfGPU
	 call UpdateEfieldGPU
#endif 	 
          call StopTimer(5)
		   
!----------------------------------------------
! Exchange Particles [Half of This step should go above Fld Update][Consider overlapping stuff here for better performance]
!----------------------------------------------  		      

            call StartTimer(12)
#ifdef CPU			
     call ExchangePrtl !Test particles are also exchanged 
#endif	 

#ifdef gpu
     call ExchangePrtlGPU
#endif 
           call StopTimer(12)
		   
		   

           call StartTimer(13)
#ifdef CPU		   
     call AppendParticles
	 call AppendTestParticles
#endif	 

           call StopTimer(13) 


!----------------------------------------------
! Exchange Current,Consider overlapping communication
!----------------------------------------------  	
           call StartTimer(14)
#ifdef CPU		   
     call UpdateCurrentsAllEdges
#endif 	 

#ifdef gpu
     call UpdateCurrentsAllEdgesGPU
#endif   
           call StopTimer(14)

!----------------------------------------------
! Some user mandated operations on the current matrix 
!----------------------------------------------  	
	 !This subroutine is defined in setup_*.F90	   NOTE: if needed, the GPU routines should be defined and called from within
	 call PreAddCurrent! Allows for some setup dependent operations before '-current' is finally added to E. For example filtering	   
     
	 call AddExternalCurrent
 !----------------------------------------------
 ! Filter Current
 !----------------------------------------------  
            call StartTimer(17)
#ifdef CPU
     call smoothen_current !default moving average filtering of the current
#endif

#ifdef gpu
     call smoothen_current_gpu
#endif
 	 
            call StopTimer(17)
			
!----------------------------------------------
! Add -Current to the electric field 
!---------------------------------------------- 
#ifdef CPU 
     call AddCurrent
#endif 	 
#ifdef gpu 
	 call AddCurrentGPU
#endif

!----------------------------------------------
! Exchange the updated value of E,B at the edges 
!---------------------------------------------- 	  
	  
          call StartTimer(18)
#ifdef CPU	  
     call ExchangeYZEdgeField(Ex,Ey,Ez)
     call ExchangeYZEdgeField(Bx,By,Bz) 	 
     call ExchangeZXEdgeField(Ex,Ey,Ez)
     call ExchangeZXEdgeField(Bx,By,Bz)
#ifndef twoD      
     call ExchangeXYEdgeField(Ex,Ey,Ez)
     call ExchangeXYEdgeField(Bx,By,Bz)
#endif 
#endif    

#ifdef gpu		   
		   call EMfldExchangeGPU
#endif

          call StopTimer(18)


!----------------------------------------------
!Filter the Electric field for moving test particles
!---------------------------------------------- 
#ifdef CPU		
           call StartTimer(20) 
    call SetFilteredEfield      
           call StopTimer(20)
#endif
!----------------------------------------------
!Enforce the boundary conditions
!---------------------------------------------- 

    call EnforceBC_Final ! calls to gpu routines are nested within	

!----------------------------------------------
!Setup specific steps
!---------------------------------------------- 			   
    call Finalsubroutines ! case dependent and is supplied by the setup module
	
!----------------------------------------------
!Balance the load across all MPI ranks, if needed
!---------------------------------------------- 
    call BalanceLoad 
		
!----------------------------------------------
!save the output 
!---------------------------------------------- 
      !call SaveOutputMaster !all data is transferred to and saved by Master
	      call StartTimer(21)
#ifdef gpu 
    call CopyGPUDataToHostForOutput
#endif 
   		  
	call SaveOutputCollective !CPU subroutine 
	call SaveRestartData    

#ifdef gpu 
    call FinishSaveDataGPU 
#endif 	
	      call StopTimer(21)

!----------------------------------------------
!Reset Current Particles
!---------------------------------------------- 		  
#ifdef CPU 		   
    call ResetCurrent
#endif 	
#ifdef gpu 	   
	call ResetCurrentGPU
#endif


!----------------------------------------------
!Reorder Particles in the memory
!---------------------------------------------- 	
	
          call StartTimer(19)
#ifdef CPU     
	    call ReorderPrtl
#endif 	
#ifdef gpu 	   
		call ReorderPrtlGPU
#endif	
	
          call StopTimer(19)
		  call StopTimer(32) !loadbalancing time 
	
	
!----------------------------------------------
!Some final subroutines 
!---------------------------------------------- 		
		   
    !call RegulatePrtlMemory
  

	
         call StopTimer(31) 
    call StepInfo !output the profiling info 
    call SavePerformanceData !save the execution time to analyse performance of the code 
	
          
end do !!!end of MAIN loop!!! 

    call StopTimer(1)
if(proc.eq.0) then
     print*,'Simulation Finished!'
     print*,'Total Running Time:',exec_time(1),'sec'
end if

call mpi_finalize(ierr)


contains 

!subroutines to run once at the begining of the simulation      
     subroutine Startups
#ifndef cyl 		 
          call CheckNum_MPI_Task
#endif		  
          call StartupMessage
		  call Initialisations	
  		    
          if(proc.eq.0) then
			   call SaveParameters
		       call SaveSetupParameters
		  end if 
		  
          call ReorderPrtl
		  
          if(.not.restart) then 
               t=0			   
               call SaveOutputCollective ! It saves the intial condition with extension _0
			   call ResetCurrent
#ifdef gpu
               call ReorderPrtlGPU
               call ResetCurrentGPU
#endif			   
          end if
          if(nMoverEMFilter.gt.0) call SetFilteredEfield ! in Fields.F90      		       
     end subroutine Startups
    subroutine CheckNum_MPI_Task
#ifdef twoD          
          if(nproc.ne.nSubDomainsX*nSubDomainsY) then
                    STOP 'Error :: Number of MPI tasks does not match (nSubDomainsX*nSubDomainsY) the number declared in parameter.F90 file'
          end if
#else 
          if(nproc.ne.nSubDomainsX*nSubDomainsY*nSubDomainsZ) then
                    STOP 'Error :: Number of MPI tasks does not match (nSubDomainsX*nSubDomainsY*nSubDomainsZ) the number declared in parameter.F90 file'
          end if
#endif     
     end subroutine CheckNum_MPI_Task
     
!Subroutines to print info about the simulation      
     subroutine StartupMessage
        
        if(proc.eq.0) then
#ifdef twoD
            print*,'            Started a 2D PIC sumulation with ',nproc,' MPI tasks'
#else
            print*,'            Started a 3D PIC sumulation with',nproc,' MPI tasks'            
#endif
    print*,'---------------------------------------------------------------------------'
    print*,'---------------------------------------------------------------------------'     
        
        end if
          
     end subroutine StartupMessage
	 
	 subroutine Initialisations
		  call InitRandomNumberSeed
		  call InitScale
		  call InitOpenMP
#ifdef cyl
          call InitParam_cyl
#endif		  

	      call Initboundaries
#ifdef cyl			  
		  call Initboundaries_cyl
#endif			   

		   
		  call InitDomainSkelton
#ifdef cyl 
		  call UpdateDomainSkelton_cyl
#endif
          
		   call InitDomainSize
#ifdef cyl			  
		   call IncAxisDomainSize
#endif			  

		  call InitPrtlBoundaries
		  call InitBCpos
		  
#ifdef cyl		  
		  call InitPrtlArrSize_cyl
#else
		  call InitPrtlArrSize
#endif 
          
		  call InitSaveOutput		  	  		  
          call InitAll    
#ifdef cyl
          call InitAll_cyl
#endif			  
		  call InitMoveDeposit 		

          call InitSetup ! Note: InitUser in the setup files must be designed to be consistent with this call in a restart situation 
		  call BoundaryConditions

          call InitOverride
#ifdef cyl
          call ComplyBC_cyl 
#else 
          call ComplyBC		 
#endif 
		  call InitPrtlScan
		  call InitPrtlStatVars

		   
		   call InitPrtlTransferInOutArr
#ifdef cyl 
           call InitPrtlTransferInOutArr_cyl 
#endif		   
		   
#ifdef gpu
           call InitAll_gpu
#ifdef cyl
           call InitAll_cyl_gpu 
#endif		   
#endif 
		  
	 end subroutine Initialisations
	 
	 subroutine RestartSteps
		 call InitRandomNumberSeed
		 call InitScale
		 call InitOpenMP
#ifdef cyl
         call InitParam_cyl
#endif
         call InitRestart
		 call RestartMPI
		 call RestartInitParam
		 call RestartAllVars
		 
		 call InitPrtlBoundaries
		 
		 call AllocateFldVars
		 call RestartAllocatePrtlMem 
		 
		 call RestartInitPrtl
		 call RestartInitFlds
		 call InitCurrent
		 call InitTransferVars
		 call InitPrtlStatVars
		 
		 call BoundaryConditions
		 call InitSaveOutput
		 call InitMoveDeposit 
		 call InitPrtlTransferInOutArr	
#ifdef cyl 
         call InitPrtlTransferInOutArr_cyl 
		 call InitAll_cyl
		 call SetRborders
#endif	
	 
#ifdef gpu
         call InitAll_gpu
#ifdef cyl
         call InitAll_cyl_gpu 
#endif		   
#endif 		
		 
	 end subroutine RestartSteps
	 
     subroutine StepInfo
          if(proc.eq.0) then
#ifdef gpu	  
               print*,'Step No. :',t,'Total Time: ',real(exec_time(31)) , 'No. Prtl :',np,'(on GPU ',np_gpu,') Prtl Arr Size:',prtl_arr_size
#else 
               print*,'Step No. :',t,'Total Time: ',real(exec_time(31)) , 'No. Prtl :',np,'Prtl Arr Size:',prtl_arr_size
#endif 			   
               print*,'Flds Update :',real(exec_time(2)+exec_time(5)),'Move Particles :',real(exec_time(4)),'Prtl Exchange:',real(exec_time(12))
               print*,'Smoothening :', real(exec_time(17)),'Reordering',real(exec_time(19)),'Append Prtl',real(exec_time(13))
               print*,'Current Exchange',real(exec_time(14)),'Fld Exchange :',real(exec_time(18)),'Load Outliers :',real(exec_time(9)),'Output:',real(exec_time(21))
               print*,'----------------------------------------------------------------'         
          end if
     end subroutine StepInfo

end program MAIN 