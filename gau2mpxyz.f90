!this program is to call the molpro WF in DFT function for high level computation in ONIOM
program gau2mpxyz
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
    character(len=20) :: chrg, uhf, solvation, atomstr,numstr
    character(len=120) :: filein, fileout, filemsg
    character(len=120) :: commandline
    character(len=80) :: c1, c2, c3, c4, c5, c6
    doubleprecision :: x, y, z, ene
    doubleprecision, allocatable :: ddip(:), hess(:,:)
    logical :: alive

    1001 FORMAT(4I10)
    1002 FORMAT(I10, 3F20.12)

    narg = iargc()
    !no extra option
    if (narg == 2) then
        call getarg(1, filein)
        call getarg(2, fileout)
        !read atoms derivatives-requested  charge  spin
        open(11,file = trim(filein), action='read')
        open(12, file = trim(fileout), status='replace')

        read(11,1001) atoms, deriva, charge, spin
        write(numstr, *) atoms 
        write(12, *) trim(adjustl(numstr))
        write(12, *) 'Gaussian external xyz file for Molpro'
        do i = 1, atoms
            read(11, 1002) atomic, x, y, z
            write(numstr,*) i 
            atomstr=trim(element(atomic))//trim(adjustl(numstr))
            !write(*,*) atomstr
            write(12, "(A,1X,3F20.12)") trim(atomstr), x*b2a, y*b2a, z*b2a
        enddo
        close(11)
        close(12)
        write(*,'(A12, 3A8)') 'atom number','deriva', 'charge', 'spin'
        write(*,'(I12, 3I8)') atoms, deriva, charge, spin
        write(*,*) 'Molpro xyz has been written as ', trim(fileout)
    else 
        write(*,*) 'input: inputfile, outputfile'
    endif


endprogram