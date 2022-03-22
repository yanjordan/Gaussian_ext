program gau_srcf
    implicit none
    DOUBLEPRECISION, PARAMETER :: b2a=0.529177249D0
    !CHARACTER(len=2), PARAMETER :: element(0:118)=(/'Bq','H ','He',& !0 is ghost atom, 1~2
    !'Li','Be','B ','C ','N ','O ','F ','Ne',& !3~10
    !'Na','Mg','Al','Si','P ','S ','Cl','Ar',& !11~18
    !'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr',& !19~36
    !'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe',& !37~54
    !'Cc','Ba',& !55~56
    !'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu',& !La, 57~71
    !'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn',& !72~86
    !'Fr','Ra',& !87~88
    !'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr',& !Ac, 89~103
    !'Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/) !104~118
    integer :: narg, i, j
    integer :: atoms, deriva, charge, spin
    integer :: ncol, nblock, nb, ne
    character(len=20) :: chrg, uhf, solvation
    character(len=120) :: filein, fileout, filemsg, filehead, filetail
    character(len=120) :: commandline
    character(len=80) :: c1, c2, c3, c4, c5, c6
    doubleprecision :: x, y, z, ene
    !doubleprecision, dimension(3) :: dip
    !doubleprecision, allocatable :: ddip(:), hess(:,:)
    logical :: alive, isgrad, ishess

    integer, dimension(:), allocatable :: atomic
    doubleprecision, dimension(:,:), allocatable :: atomiccoord, grad

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

    narg = iargc()
    if(narg == 8) then
        call getarg(1, filehead)
        call getarg(2, filetail)
        call getarg(4, filein) 
        call getarg(5, fileout)
        call getarg(6, filemsg)

            !read gaussian input
        open(11,file = trim(filein), action='read')
        open(12,file = trim(filemsg), status='replace')
        read(11,1001) atoms, deriva, charge, spin

        write(12, *) atoms, deriva, charge, spin

        allocate(atomic(atoms), atomiccoord(3,atoms), grad(3,atoms))
        !add coordinates
        do i = 1, atoms
            read(11, 1002) atomic(i), atomiccoord(:,i) 
        enddo
        close(11)

        call rungauss(filehead, filetail, atoms, atomic, atomiccoord, ene, grad)

        open(14,file = trim(fileout), status='replace')
        !energy, dipole-moment (xyz)
        write(14, "(4D20.12)") ene, 0D0,0D0,0D0
        do i=1, atoms
            write(14,"(3D20.12)") grad(:,i)
        enddo
        close(14)
        deallocate(atomic,atomiccoord,grad)
     
        close(12)
    else
        write(*,*) "8 inputs: headf, tailf, Gaussian external 6 input"
        write(*,*) "headf for Gaussian gjf including force punch=(coord,derivatives) "
        write(*,*) "tailf for Gaussian gjf for extral bs or info"
        write(*,*) "Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile"
        write(*,*) "More infomation see Gaussian website"
    endif

endprogram

!run gaussian job using head tail and atomic coordnations
subroutine rungauss(hfile, tailf, numatom, atomic, atomcoord, ene, grad)
    implicit none

    DOUBLEPRECISION, PARAMETER :: b2a=0.52917721067121D0
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


    character(len=*), intent(in) :: hfile, tailf
    integer, intent(in) :: numatom
    integer, dimension(numatom), intent(in) :: atomic
    doubleprecision, dimension(3,numatom), intent(in) :: atomcoord
    doubleprecision, intent(out) :: ene 
    doubleprecision, dimension(3,numatom), intent(out) :: grad

    integer :: i, j, isend
    logical :: alive, b
    character(len=120) :: line, c1, c2, c3

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)
    !prepare the gaussian gjf file
    inquire(file=trim(hfile), exist=alive)
    if (.not. alive) then
        write(*,*) trim(hfile)//' does not exist'
        stop 
    endif 
        !read the head  file
    open(30,file=trim(hfile),action='read')
        !write gjf of fragment 
    open(40,file='Center_temp.gjf',status='replace')

    b = .false.
    do while(.not. b)
        read(30,'(A)') line
        write(40,*) trim(line)
        b = (trim(line) == '')
    enddo

    b = .false.
    do while(.not. b)
        read(30,'(A)') line
        write(40,*) trim(line)
        b = (trim(line) == '')
    enddo

    read(30,'(A)') line
    write(40,*) trim(line)
    close(30)

    !write atom block
    do i=1,numatom 
        write(40,'(1X,A2,13X,3F14.8)') element(atomic(i)), atomcoord(:,i)*b2a
    enddo

    !read the tail
    inquire(file=trim(tailf), exist=b)
                    
    !withno extral information
    if (.not. b) then
        write(40,*) ' '
        write(40,*) ' '
        write(40,*) ' '
    else
        open(31,file=trim(tailf),action='read')
        do while(.true.)
            read(31, '(A)', IOSTAT = isend) line
            if (isend < 0) then
                EXIT
            else
                write(40,*) trim(line)
            endif
        enddo
        close(31)
    endif
    write(40,*) ' '
    close(40)

    !run Gaussian job
    call system('cat Center_temp.gjf >> center_gjf.bak')
    call system('echo >> center_gjf.bak')
    !run gaussian job
    call system('$Gauexe < Center_temp.gjf > Center_temp.out')
    !extract energy and force 
    call system("grep 'extrapolated energy' Center_temp.out | awk '{print $5}' > result.temp ")
    open(16, file='result.temp',action='read')
    if(EOF(16)) then
        close(16)
        call system("grep 'SCF Done:' Center_temp.out | awk '{print $5}' > result.temp ")
    else
        close(16)
    endif

    !get energy
    open(20, file='result.temp', action='read')
    read(20,*) ene
    close(20)

    open(21,file='fort.7', action='read')

    do i= 1, numatom
        read(21,*)
    enddo

    do i= 1, numatom
        read(21,*) grad(:,i)
    enddo
    close(21)

endsubroutine


