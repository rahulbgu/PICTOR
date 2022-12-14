module interp_gpu
	use parameters 
	use vars 
	use var_gpu 
	use cudafor
contains 
	
	subroutine SetTexFlds
		call InitTexExGPUKernel<<<grid,tBlock>>>(TexEx_gpu,Ex_gpu,mx,my,mz)
		call InitTexEyGPUKernel<<<grid,tBlock>>>(TexEy_gpu,Ey_gpu,mx,my,mz)
		call InitTexEzGPUKernel<<<grid,tBlock>>>(TexEz_gpu,Ez_gpu,mx,my,mz)
		call InitTexBxGPUKernel<<<grid,tBlock>>>(TexBx_gpu,Bx_gpu,mx,my,mz,Bx_ext0)
		call InitTexByGPUKernel<<<grid,tBlock>>>(TexBy_gpu,By_gpu,mx,my,mz,By_ext0)
		call InitTexBzGPUKernel<<<grid,tBlock>>>(TexBz_gpu,Bz_gpu,mx,my,mz,Bz_ext0)
	end subroutine SetTexFlds
	
	attributes(global) subroutine InitTexExGPUKernel(TexFld,Fld,mx,my,mz)
 		integer, value :: mx,my,mz
 		real, dimension(mx,my,mz) :: TexFld,Fld 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +1
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 		
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.5*(Fld(i-1,j,k)+Fld(i,j,k))
	end subroutine InitTexExGPUKernel
	
	attributes(global) subroutine InitTexEyGPUKernel(TexFld,Fld,mx,my,mz)
 		integer, value :: mx,my,mz	
 		real, dimension(mx,my,mz) :: TexFld,Fld 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +1
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 		
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.5*(Fld(i,j-1,k)+Fld(i,j,k))
	end subroutine InitTexEyGPUKernel
	
	attributes(global) subroutine InitTexEzGPUKernel(TexFld,Fld,mx,my,mz)
 		integer, value :: mx,my,mz		
 		real, dimension(mx,my,mz) :: TexFld,Fld		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +1
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.5*(Fld(i,j,k-1)+Fld(i,j,k))
#else
        k=1
		if((i.le.mx).and.(j.le.my)) TexFld(i,j,k)=Fld(i,j,k)
#endif 		
	end subroutine InitTexEzGPUKernel
	
	attributes(global) subroutine InitTexBxGPUKernel(TexFld,Fld,mx,my,mz,extB)
 		integer, value :: mx,my,mz
 		real, dimension(mx,my,mz) :: TexFld,Fld		
        real, value :: extB 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x 
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y+1 
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z+1 
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.25*(Fld(i,j-1,k-1)+Fld(i,j,k-1)+Fld(i,j-1,k)+Fld(i,j,k)) + extB
#else
        k=1
		if((i.le.mx).and.(j.le.my)) TexFld(i,j,k)=0.5*(Fld(i,j-1,k)+Fld(i,j,k)) +extB
#endif 		
	end subroutine InitTexBxGPUKernel
	attributes(global) subroutine InitTexByGPUKernel(TexFld,Fld,mx,my,mz,extB)
 		integer, value :: mx,my,mz	
 		real, dimension(mx,my,mz) :: TexFld,Fld
        real, value :: extB 		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +1
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z +1
		if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.25*(Fld(i-1,j,k-1)+Fld(i-1,j,k)+Fld(i,j,k-1)+Fld(i,j,k)) +extB
#else
        k=1
		if((i.le.mx).and.(j.le.my)) TexFld(i,j,k)=0.5*(Fld(i-1,j,k)+Fld(i,j,k)) +extB
#endif 		
	end subroutine InitTexByGPUKernel
	attributes(global) subroutine InitTexBzGPUKernel(TexFld,Fld,mx,my,mz,extB)
 		integer, value :: mx,my,mz
 		real, dimension(mx,my,mz) :: TexFld,Fld 
        real, value :: extB		
 		integer :: i,j,k
		
 		i = (blockIdx%x-1)*blockDim%x + threadIdx%x +1
 		j = (blockIdx%y-1)*blockDim%y + threadIdx%y +1
#ifndef twoD 		
 		k = (blockIdx%z-1)*blockDim%z + threadIdx%z 
#else
        k=1
#endif 		
        if((i.le.mx).and.(j.le.my).and.(k.le.mz)) TexFld(i,j,k)=0.25*(Fld(i-1,j-1,k)+Fld(i-1,j,k)+Fld(i,j-1,k)+Fld(i,j,k)) +extB
	end subroutine InitTexBzGPUKernel
end module interp_gpu