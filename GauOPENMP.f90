
program Gaussian_par
    USE OMP_LIB
  
    implicit none
    integer :: narg, i
    CHARACTER(len=120), DIMENSION(:), allocatable :: files
    CHARACTER(len=120) :: fin, fout, cmd
    narg = iargc()
  
    if(narg==0) THEN
      WRITE(*,*) 'This is a program to run Gaussian jobs on parralel.'
      WRITE(*,*) 'Please define the Gauexe before: g16 or g09'
      WRITE(*,*) 'Usage: gjf1 gjf2 ... gjfn'
  
    else
      allocate(files(narg))
      !get all the gjf names
      do i=1,narg
        call getarg(i, files(i))
      enddo
  
      !$OMP PARALLEL PRIVATE(fin,fout,cmd)
      !$OMP DO SCHEDULE(RUNTIME)
      do i=1,narg
        fin=files(i)
        fout=fin(1:len_trim(fin)-4)//'.log'
        cmd='$Gauexe < '//trim(fin)//' > '//trim(fout)
        WRITE(*,*) cmd 
        call system(trim(cmd))
        WRITE(*,*) 'Computation of '//trim(fin)//'is finished!'
      enddo
      !$OMP END DO
      !$OMP END PARALLEL
  
  
    endif
  
  end program Gaussian_par