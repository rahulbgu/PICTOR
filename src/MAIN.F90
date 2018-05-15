program MAIN
     use parameters 
     use vars
     use initialise
     use movdep
     use particles
     use fields
     use savedata
     use comm_fldprtl
     use memory
	 use reload
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
     !use cudafor
	 use var_gpu
	 use initialise_gpu
	 use memory_gpu 
	 use comm_fld_gpu
	 use movdep_gpu 
	 use fields_gpu
#endif  

     implicit none

call StartMPI 
if(restart) call RestartMPI
call Startups
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
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
     if(use_device) call UpdateBfieldHalfGPU
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
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
    if(use_device) call MoveDepositPrtlGPU
    if(use_device) call MoveTestPrtlGPU	
#endif 	 

           call StopTimer(4)

#ifdef GPU_EXCLUSIVE
     call ExchangePrtlGPU_Exclusive
#endif 
		   
!----------------------------------------------
! Some operations performed right after moving the particles 
!----------------------------------------------   		   
!Not available in the GPU version
#ifdef CPU		  
	 !This subroutine is defined in setup_*.F90
     call PostMovDep ! It allows for some setup dependent operations right after moving the particles and depositing their current  
     
#endif 	 
 !----------------------------------------------
 !Update magnetic field by another half step and electric field by full step 
 !----------------------------------------------  	  
           call StartTimer(5)
#ifdef CPU		   
     call UpdateBfieldHalf 
#endif	 
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
     if(use_device) call UpdateBfieldHalfGPU
#endif 	 
          call StopTimer(5)
 
           call StartTimer(7)
#ifdef CPU 		   
     call UpdateEfield
#endif	 
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
  if(use_device) call UpdateEfieldGPU
#endif
           call StopTimer(7)
		   
		   
!----------------------------------------------
! Exchange Particles [Half of This step should go above Fld Update][Consider overlapping stuffs here]
!----------------------------------------------  		      

            call StartTimer(12)
#ifdef CPU			
     call ExchangePrtl !Test particles are also exchanged 
#endif	 
#ifdef gpu  
	  if(use_device) then 
		  call ExchangePrtlGPU
		  call ExchangeTestPrtlGPU
	  end if 
#endif 
           call StopTimer(12)
		   
		   

           call StartTimer(13)
#ifdef CPU		   
     call AppendParticles
	 call AppendTestParticles
#endif	 
#ifdef gpu
     if(use_device) then 		
		call AppendHostPrtlArr
		call AppendGPUPrtlArr
		call AppendHostTestPrtlArr
		call AppendGPUTestPrtlArr	
	 end if 
#endif 
#ifdef GPU_EXCLUSIVE
        call AppendGPUPrtlArr
		call AppendGPUTestPrtlArr
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
     if(use_device) call UpdateCurrentsAllEdgesGPU
#endif 	
#ifdef GPU_EXCLUSIVE
    call UpdateCurrentsAllEdgesGPU_Exclusive
#endif   
           call StopTimer(14)

!----------------------------------------------
! Some user mandated operations on the current matrix 
!----------------------------------------------  	
#ifdef CPU	
	 !This subroutine is defined in setup_*.F90	   NOTE: Not avaialble in GPU version
	 call PreAddCurrent! Allows for some setup dependent operations before '-current' is finally added to E. For example filtering	   
#endif 

 !----------------------------------------------
 ! Filter Current
 !----------------------------------------------  
            call StartTimer(17)
#ifdef CPU
     call smoothen_current !default moving average filtering of the current
#endif
#ifdef gpu
     if(use_device) then 
		  call smoothen_current_gpu
	 else 
		  call smoothen_current
	 end if 
#endif
#ifdef GPU_EXCLUSIVE
      call smoothen_current_gpu_exclusive
#endif
 	 
            call StopTimer(17)
			
!----------------------------------------------
! Add -Current to the electric field 
!---------------------------------------------- 
#ifdef CPU 
     call AddCurrent
#endif 	 
#ifdef gpu 
	 if(use_device) call AddCurrentGPU
#endif
#ifdef GPU_EXCLUSIVE
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
		   if(use_device) call ExchangeEMfldGPU
#endif

#ifdef GPU_EXCLUSIVE
      call EMfldExchangeGPU_Exclusive
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
!save the output 
!---------------------------------------------- 
      !call SaveOutputMaster !all data is transferred to and saved by Master
	      call StartTimer(21)
#ifdef gpu 
    if(use_device) call CopyGPUDataToHostForOutput
#endif 
#ifdef GPU_EXCLUSIVE
    call CopyGPUDataToHostForOutput
#endif 		   
		  
	call SaveOutputCollective !CPU subroutine 

#ifdef gpu 
    if(use_device) call FinishSaveDataGPU 
#endif 	
#ifdef GPU_EXCLUSIVE
     call FinishSaveDataGPU 
#endif 
	      call StopTimer(21)

!----------------------------------------------
!Reset Current, Reorder Particles
!---------------------------------------------- 		  
#ifdef CPU 		   
    call ResetCurrent
#endif 	
#ifdef gpu 	   
		if(use_device) call ResetCurrentGPU
#endif
#ifdef GPU_EXCLUSIVE
        call ResetCurrentGPU
#endif 

	
          call StartTimer(19)
#ifdef CPU     
	call ReorderPrtl
#endif 	
#ifdef gpu 	   
		if(use_device) call ReorderPrtlGPU
		if(use_device) call ReorderTestPrtlGPU
#endif	
#ifdef GPU_EXCLUSIVE
        call ReorderPrtlGPU
		call ReorderTestPrtlGPU
#endif	
          call StopTimer(19)
		  call StopTimer(32) !loadbalancing time 
	
	
!----------------------------------------------
!Some final subroutines 
!---------------------------------------------- 		
! #ifdef gpu
!     call AdjustDomainGPU
! #endif

		  
    call SaveRestartData     
    !call RegulatePrtlMemory
    call Finalsubroutines ! case dependent and is supplied by setup module
	
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
          call StartupChecks
          call StartupMessage
          call InitSaveOutput		  
          call InitAll    
		  
		  call InitMoveDeposit 		
#if defined(gpu) || defined(GPU_EXCLUSIVE) 
          call InitAll_gpu
#endif 		  
		    
          if(proc.eq.0) call SaveParameters
          call ReorderPrtl
		  		  
          if(.not.restart) then 
               t=0			   
               call SaveOutputCollective ! It saves the intial condition with extension _0
			   call ResetCurrent
#if defined(gpu) || defined(GPU_EXCLUSIVE)
               call ReorderPrtlGPU
               call ReorderTestPrtlGPU
               call ResetCurrentGPU
#endif			   
          end if
          if(nMoverEMFilter.gt.0) call SetFilteredEfield ! in Fields.F90           
     end subroutine Startups
    subroutine StartupChecks
#ifdef twoD          
          if(nproc.ne.nSubDomainsX*nSubDomainsY) then
                    STOP 'Error :: Number of cores assgined does not match (nSubDomainsX*nSubDomainsY) the number declared in parameter.F90 file'
          end if
#else 
          if(nproc.ne.nSubDomainsX*nSubDomainsY*nSubDomainsZ) then
                    STOP 'Error :: Number of cores assgined does not match (nSubDomainsX*nSubDomainsY*nSubDomainsZ) the number declared in parameter.F90 file'
          end if
#endif     
               
     end subroutine StartupChecks
     
!Subroutines to print info about the simulation      
     subroutine StartupMessage
        
        if(proc.eq.0) then
#ifdef twoD
            print*,'            Started a 2D PIC sumulation with ',nproc,' MPI jobs'
#else
            print*,'            Started a 3D PIC sumulation with',nproc,' MPI jobs'            
#endif
    print*,'---------------------------------------------------------------------------'
    print*,'---------------------------------------------------------------------------'     
        
        end if
          
     end subroutine StartupMessage
     subroutine StepInfo
          if(proc.eq.0) then
#if defined(gpu) || defined(GPU_EXCLUSIVE) 		  
               print*,'Step No. :',t,'Total Time: ',real(exec_time(31)) , 'No. Prtl :',np,'(on GPU ',np_gpu,') Prtl Arr Size:',prtl_arr_size
#else 
               print*,'Step No. :',t,'Total Time: ',real(exec_time(31)) , 'No. Prtl :',np,'Prtl Arr Size:',prtl_arr_size
#endif 			   
               print*,'Flds Update :',real(exec_time(2)+exec_time(5)+exec_time(7)),'Move Particles :',real(exec_time(4)),'Prtl Exchange:',real(exec_time(12))
               print*,'Smoothening :', real(exec_time(17)),'Reordering',real(exec_time(19)),'Append Prtl',real(exec_time(13))
               print*,'Current Exchange',real(exec_time(14)),'Fld Exchange :',real(exec_time(18)),'Load Outliers :',real(exec_time(9)),'Output:',real(exec_time(21))
               print*,'----------------------------------------------------------------'         
          end if
     end subroutine StepInfo

end program MAIN 