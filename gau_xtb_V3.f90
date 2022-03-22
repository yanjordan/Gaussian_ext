!this external program was written by Dr. YAN Zeyin at 2019-04-17, Southern University of Science and Technology
!contact: zeyin.yan@outlook.com

program gau2xtb
    implicit none
    DOUBLEPRECISION, PARAMETER :: b2a=0.529177249D0
    CHARACTER(len=2), PARAMETER :: element(0:118)=(/'Bq','H ','He',& !0 is ghost atom, 1~2
        'Li','Be','B ','C ','N ','O ','F ','Ne',& !3~10
        'Na','Mg','Al','Si','P ','S ','Cl','Ar',& !11~18
        'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',& !19~36
        'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',& !37~54
        'Cc','Ba',& !55~56
        'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',& !La, 57~71
        'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',& !72~86
        'Fr','Ra',& !87~88
        'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',& !Ac, 89~103
        'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/) !104~118
    integer :: narg, i, j, atomic
    integer :: atoms, deriva, charge, spin
    character(len=20) :: chrg, uhf
    character(len=120) :: filein, fileout, filemsg, filecon
    character(len=120) :: commandline
    character(len=80) :: c1, c2, c3, c4, c5, c6
    doubleprecision :: x, y, z, ene
    doubleprecision, allocatable :: ddip(:), hess(:,:)
    logical :: alive

    !V3 add time info
    integer :: int1, int2 

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

    call system_clock(int1)
    narg = iargc()
    !no extra option
    if(narg == 6) then
        call getarg(2, filein) 
        call getarg(3, fileout)
        call getarg(4, filemsg)
        !read atoms  derivatives-requested  charge  spin
        open(10, file = trim(filemsg), status='replace')
        write(10,*) 'begin read file: ',trim(filein) 
        write(10,*) 'and write the xyz file'

        open(11,file = trim(filein), action='read')
        open(12, file = 'gau_xtb.xyz', status='replace')

        read(11,1001) atoms, deriva, charge, spin
        write(12, '(I5)') atoms 
        write(12, '(A)') 'Gau_xtb'
        do i = 1, atoms
            read(11, 1002) atomic, x, y, z
            write(12, "(A2,1X,3F20.12)") element(atomic), x*b2a, y*b2a, z*b2a
        enddo
        close(11)
        close(12)
        write(10,*) 'xyz file is generated as gau_xtb.xyz'
 
        !lanch xtb computation
        call system('rm -f charges energy xtbrestart gradient hessian')
        call system('rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord')
        write(chrg, *) charge 
        write(uhf, *) spin-1
     
        if (deriva == 0) then
            commandline='xtb gau_xtb.xyz --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --sp > gau_xtb.out' 
        elseif (deriva == 1) then
            commandline='xtb gau_xtb.xyz --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --grad > gau_xtb.out'          
        elseif (deriva == 2) then
            commandline='xtb gau_xtb.xyz --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --hess --grad > gau_xtb.out'
        endif

        write(10,*) 'xtb lanch: '//trim(commandline)
        call system(trim(commandline))
        write(10,*) 'xtb job finished!'
       
        call system("grep 'TOTAL ENERGY' gau_xtb.out | awk '{print $4}' > energy.tmp")
        open (9, file='energy.tmp', action='read')
            read(9,*) ene
        close(9)
        call system("rm -f energy.tmp")

        open(14,file=trim(fileout),status='replace')
        write(14,"(4D20.12)") ene,0D0,0D0,0D0

        write(10,*) 'Energy has been extracted'
        
        if (deriva == 1 .or. deriva==2) then
            inquire(file='gradient',exist=alive)
            if (.not. alive) then
                write(*,*) 'Error: gradient is not existed in current folder!'
                stop
            endif

            open(13, file='gradient',action='read')          
            read(13,*)
            read(13,*)            
            do i = 1, atoms
                read(13,*) 
            enddo
            do i = 1, atoms
                read(13,*) x, y, z
                write(14,"(3D20.12)") x, y, z
            enddo
            close(13)
            write(10,*) 'XTB energy: ', ene, ' Eh'
            write(10,*) 'Gradient has been extracted'
        endif


        if (deriva==2) then
            !polarizability
            write(14, "(3D20.12)") 0.0, 0.0, 0.0            
            write(14, "(3D20.12)") 0.0, 0.0, 0.0
            !dipole derivatives
            allocate(ddip(9*atoms))
            ddip=0.0
            write(14,"(3D20.12)") ddip
            deallocate(ddip)
            inquire(file='hessian',exist=alive)
            if(.not. alive) then
                write(*,*) "Error: hessian is not existed in current folder!"
                stop
            endif
            allocate(hess(3*atoms,3*atoms))
            open(15,file='hessian',action='read')
            read(15,*)
            read(15,*) ((hess(i,j),j=1,3*atoms),i=1,3*atoms)
            close(15)
            write(14,"(3D20.12)") ((hess(i,j),j=1,i),i=1,3*atoms)
            write(10,*) 'Hessian has been extracted'
        endif

        close(14)

        call system_clock(int2)

        write(10,*) 'System time: ', (int2-int1)/1000.0
        close(10)

        call system('rm -f charges energy xtbrestart gradient hessian')
        call system('rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord')
    !add control file to add more option
    elseif(narg == 7) then
        call getarg(1, filecon)
        call getarg(3, filein) 
        call getarg(4, fileout)
        call getarg(5, filemsg)
        !read atoms  derivatives-requested  charge  spin
        open(10, file = trim(filemsg), status='replace')

        inquire(file=trim(filecon),exist=alive)
        if(.not. alive) then
            write(10,*) "File "//trim(filecon)//" does not exist!"
            stop 
        endif

        write(10,*) 'begin read file: ',trim(filein) 
        write(10,*) 'and write the xyz file'

        open(11,file = trim(filein), action='read')
        open(12, file = 'gau_xtb.xyz', status='replace')

        read(11,1001) atoms, deriva, charge, spin
        write(12, *) atoms 
        write(12, *) 'Gaussian external xyz file for xtb'
        do i = 1, atoms
            read(11, 1002) atomic, x, y, z
            write(12, "(A2,1X,3F20.12)") element(atomic), x*b2a, y*b2a, z*b2a
        enddo
        close(11)
        close(12)
        write(10,*) 'xyz file is generated as gau_xtb.xyz'
    
        !lanch xtb computation
        call system('rm -f charges energy xtbrestart gradient hessian')
        call system('rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord')
        write(chrg, *) charge 
        write(uhf, *) spin-1
        
        if (deriva == 0) then
            commandline='xtb gau_xtb.xyz '//'-I '//trim(filecon)//' --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --sp > gau_xtb.out' 
        elseif (deriva == 1) then
            commandline='xtb gau_xtb.xyz '//'-I '//trim(filecon)//' --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --grad > gau_xtb.out'          
        elseif (deriva == 2) then
            commandline='xtb gau_xtb.xyz '//'-I '//trim(filecon)//' --chrg '//trim(adjustl(chrg))//' --uhf '//trim(adjustl(uhf))//' --hess --grad > gau_xtb.out'
        endif

        write(10,*) 'xtb lanch: '//trim(commandline)
        call system(trim(commandline))
        write(10,*) 'xtb job finished!'
        
        call system("grep 'TOTAL ENERGY' gau_xtb.out | awk '{print $4}' > energy.tmp")
        open (9, file='energy.tmp', action='read')
            read(9,*) ene
        close(9)
        call system("rm -f energy.tmp")

        open(14,file=trim(fileout),status='replace')
        write(14,"(4D20.12)") ene,0D0,0D0,0D0

        write(10,*) 'Energy has been extracted'
        
        if (deriva == 1 .or. deriva==2) then
            inquire(file='gradient',exist=alive)
            if (.not. alive) then
                write(*,*) 'Error: gradient is not existed in current folder!'
                stop
            endif

            open(13, file='gradient',action='read')          
            read(13,*)
            read(13,*)            
            do i = 1, atoms
                read(13,*) 
            enddo
            do i = 1, atoms
                read(13,*) x, y, z
                write(14,"(3D20.12)") x, y, z
            enddo
            close(13)
            write(10,*) 'XTB energy: ', ene, ' Eh'
            write(10,*) 'Gradient has been extracted'
        endif


        if (deriva==2) then
            !polarizability
            write(14, "(3D20.12)") 0.0, 0.0, 0.0            
            write(14, "(3D20.12)") 0.0, 0.0, 0.0
            !dipole derivatives
            allocate(ddip(9*atoms))
            ddip=0.0
            write(14,"(3D20.12)") ddip
            deallocate(ddip)
            inquire(file='hessian',exist=alive)
            if(.not. alive) then
                write(*,*) "Error: hessian is not existed in current folder!"
                stop
            endif
            allocate(hess(3*atoms,3*atoms))
            open(15,file='hessian',action='read')
            read(15,*)
            read(15,*) ((hess(i,j),j=1,3*atoms),i=1,3*atoms)
            close(15)
            write(14,"(3D20.12)") ((hess(i,j),j=1,i),i=1,3*atoms)
            write(10,*) 'Hessian has been extracted'
        endif

        close(14)

        call system_clock(int2)
   
        write(10,*) 'System time: ', (int2-int1)/1000.0
        close(10)

        
        call system('rm -f charges energy xtbrestart gradient hessian')
        call system('rm -f hessian xtb_normalmodes g98_canmode.out g98.out wbo xtbhess.coord')
    else
        write(*,*) 'No extra input, please see keyword External'
        write(*,*) 'in Gaussian xx website for more information.'
    endif
    
endprogram
