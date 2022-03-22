program g09_orca
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
    integer :: ncol, nblock, nb, ne
    character(len=20) :: chrg, uhf, solvation
    character(len=120) :: filein, fileout, filemsg, filehead
    character(len=120) :: commandline
    character(len=80) :: c1, c2, c3, c4, c5, c6
    doubleprecision :: x, y, z, ene
    doubleprecision, dimension(3) :: dip
    doubleprecision, allocatable :: ddip(:), hess(:,:)
    logical :: alive, isgrad, ishess

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

    narg = iargc()
    if(narg == 7) then
        call getarg(1, filehead)
        call getarg(3, filein) 
        call getarg(4, fileout)
        call getarg(5, filemsg)

        !call system("cp "//trim(filein)//" .")

        !read gaussian input
        open(11,file = trim(filein), action='read')
        open(12,file = trim(filemsg), status='replace')
        read(11,1001) atoms, deriva, charge, spin

        write(12, *) atoms, deriva, charge, spin

        allocate(ddip(9*atoms), hess(3*atoms,3*atoms))
        !check the keywords in gau_orca.inp
        call system("grep -i engrad "//trim(filehead)//" > isengrad.tmp")
        open(20, file = 'isengrad.tmp')
        isgrad = EOF(20)
        close(20)
        call system("rm -f isengrad.tmp")

        call system("grep -i freq "//trim(filehead)//" > ishess.tmp")
        open(20, file = 'ishess.tmp')
        ishess = EOF(20)
        close(20)
        call system("rm -f ishess.tmp")

        !check computational type with keywords in gau_orca.inp
        if (deriva == 0) then
            write(12, *) 'Run ORCA energy computation'
        elseif (deriva == 1) then
            if (isgrad) then
                write(12, *) 'There is no engrad keyword in '//trim(filehead)
                stop
            endif
            write(12, *) 'Run ORCA engrad computation'
        elseif (deriva == 2) then
            if (isgrad .or. ishess) then
                write(12, *) 'There is no engrad and freq keywords in '//trim(filehead)
                stop
            endif
            write(12, *) 'Run ORCA engrad freq computation'
        endif

        !creat ORCA input   
        inquire(file=trim(filehead),exist=alive)
        if (alive) then
            call system("cp "//trim(filehead)//" gau_orca.inp")
        else
            !default calculation method and ncps, mem
            write(12, *) 'file '//trim(filehead)//'does not exist!'
            stop
        endif

        !add coordinates
        open(13, file='gau_orca.inp', access='append')
        write(13, *) 
        write(13, "(A5,2I4)") "* xyz", charge, spin
        do i = 1, atoms
            read(11, 1002) atomic, x, y, z
            write(13, "(A2,1X,3F20.12)") element(atomic), x*b2a, y*b2a, z*b2a
        enddo
        write(13, "(A1)") "*"
        write(13, *) 
        close(13)
        close(11)

        !run orca computation
        
        call system("$ORCA_BIN/orca gau_orca.inp > gau_orca.log")

        call system("rm -f *.tmp*")
        !get energy, force and hess
        call system("grep 'Total Energy       :' gau_orca.log | awk '{print $4}' > result_orca.tmp")
        call system("grep 'Total Dipole Moment    :' gau_orca.log | awk '{print $5,$6,$7}' >> result_orca.tmp")

        open(14,file = trim(fileout), status='replace')
        open(15,file = 'result_orca.tmp', action='read')
        read(15, *) ene 
        read(15, *) dip(1:3)
        close(15)
        call system("rm -f result_orca.tmp")
        !energy, dipole-moment (xyz)
        write(14, "(4D20.12)") ene, dip(1:3)

        

        if (deriva == 1 .or. deriva == 2) then

            

            open(16, file='gau_orca.engrad', action='read')
            do i=1,11
                read(16,*)
            enddo
            !gradient on atom (xyz)
            do i = 1,atoms 
                read(16,*) x
                read(16,*) y
                read(16,*) z
                write(14, "(3D20.12)")  x,y,z
            enddo
            close(16)
            
        endif
        write(12, *) 2
        if (deriva == 2) then
            
            !polarizability
            write(14, "(3D20.12)")  0.0, 0.0, 0.0
            write(14, "(3D20.12)")  0.0, 0.0, 0.0

            !read force constants
            open(17, file='gau_orca.hess', action='read')
            do while(.true.)
                read(17,"(A)") c1
                if (trim(c1)=='$hessian') then
                    exit
                endif
            enddo
            read(17, *)
           
            ncol=5
            nblock=ceiling(3*atoms/5D0)

            do i = 1,nblock 
                nb=(i-1)*ncol + 1
                if (i /= nblock) then
                    ne = (i-1)*ncol+ncol
                else
                    ne = 3*atoms
                endif

                read(17,*)

                do j = 1, 3*atoms
                    read(17, *) c1, hess(j, nb:ne) 
                enddo
            enddo
            
            !read dipole derivatives
            do while(.true.)
                read(17,"(A)") c1
                if (trim(c1)=='$dipole_derivatives') then
                    exit
                endif
            enddo
            read(17,*)
            do i = 1, 3*atoms 
                read(17, *) ddip((i-1)*3+1: 3*i)
            enddo
            close(17)
         

            !write dipole derivatives
            write(14, "(3D20.12)") ddip
            !write force constants
            write(14, "(3D20.12)") ((hess(i,j),j=1,i),i=1,3*atoms)

     
        endif
        deallocate(ddip, hess) 
        close(14)
       
        close(12)
    else
        write(*,*) "7 inputs: orca_input_head, Gaussian external 6 input"
        write(*,*) "orca_input_head with EnGrad for force computation "
        write(*,*) "orca_input_head with EnGrad freq for force and hessian computation "
        write(*,*) "Gau_External layer InputFile OutputFile MsgFile FChkFile MatElFile"
        write(*,*) "More infomation see Gaussian website"
    endif

endprogram