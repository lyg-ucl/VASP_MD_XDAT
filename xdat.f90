!==========================================================================
! Program: XDAT                                              
! Version: xdat.x_1.0
! Copyright (C) Yunguo Li @ UCL                              
! Email :: yunguo.li@ucl.ac.uk                               
! Function :: extract lattice parameters, coordination, RMSD and .xyz
! Usage ::  follow prompts                              
! Input ::  XDATCAR (vasp5 format)                          
! Output::  rmsd.dat, xdat.dat, xdat.xyz, coord_num_.dat
!==========================================================================
PROGRAM  XDAT
!==========================================================================
!-----------------------------------------------------------------------
! announcement
!-----------------------------------------------------------------------
IMPLICIT  NONE

INTEGER, PARAMETER :: ap = selected_real_kind(15,300)   
INTEGER :: ntype,nr,maxatom,filein,i,k,m,IOstatus=0,p,q,j,n,s,t,bondtype
REAL(ap), DIMENSION(:), ALLOCATABLE :: xarray, yarray, zarray, disp
REAL(ap), DIMENSION(:,:), ALLOCATABLE :: xarrayimag,yarrayimag,zarrayimag
REAL(ap), DIMENSION(:,:,:), ALLOCATABLE :: cutoff
REAL(ap) :: a, b, c, scale, vol, alpha, beta,gamma,x,y,z,xt,yt,zt,pi = 4.*atan(1.), bondlength
REAL(ap), DIMENSION(3) :: avect, bvect, cvect
INTEGER, DIMENSION(:), ALLOCATABLE :: atnr
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: coordnum
CHARACTER(len=255) :: line, title1,title2
CHARACTER(len=8) :: fnm
CHARACTER(len=2),dimension(:),ALLOCATABLE :: elemnm
LOGICAL  :: exist

!-----------------------------------------------------------------------
!remove old files
!-----------------------------------------------------------------------
INQUIRE(file="xdat.xyz",EXIST=exist) 
IF(exist)  CALL SYSTEM("rm xdat.xyz")
INQUIRE(file="xdat.dat",EXIST=exist) 
IF(exist)  CALL SYSTEM("rm xdat.dat")
INQUIRE(file="rmsd.dat",EXIST=exist) 
IF(exist)  CALL SYSTEM("rm rmsd.dat")
INQUIRE(file="defectfrac.dat",EXIST=exist)
IF(exist)  CALL SYSTEM("rm defectfrac.dat")
INQUIRE(file="coord.dat",EXIST=exist)
IF(exist)  CALL SYSTEM("rm coord.dat")

!-----------------------------------------------------------------------
! open XDATCAR
!-----------------------------------------------------------------------
filein = 12
OPEN(filein, file='XDATCAR', status='old')
OPEN(16,file='xdat.dat', status='new')
WRITE(16,"(A6,7A16)") "#Index","a","b","c","alpha","beta","gamma","volume"  
OPEN(26,file='xdat.xyz', status='new')
OPEN(36,file='rmsd.dat', status='new')

!-----------------------------------------------------------------------
! start main loop
!-----------------------------------------------------------------------
i=0
DO WHILE (IOstatus .eq. 0)
  i=i+1
  IF(i==1) THEN
    READ(filein,*,IOSTAT=IOstatus) title1 ! read title
    title2=title1
  END IF

  IF (IOstatus .ne. 0) EXIT
  ! read XDATCAR format  
  IF(i>1) THEN
    READ(filein,*,IOSTAT=IOstatus) title2
    IF (IOstatus .ne. 0) EXIT              ! reading EOF, exit
    IF ((title2 .ne. title1) .and. i==3) WRITE(*,*) "VASP NVT/NVE FORMAT"
    IF (title2 .ne. title1)  GOTO 100      
    IF ((title2 .eq. title1) .and. i==3) WRITE(*,*) "VASP NPT FORMAT"    
  END IF
  ! read lattice
  READ(filein,*) scale
  READ(filein,*) avect(1), avect(2), avect(3)
  a = scale*sqrt(avect(1)*avect(1) + avect(2)*avect(2) + avect(3)*avect(3))
  READ(filein,*) bvect(1), bvect(2), bvect(3)
  b = scale*sqrt(bvect(1)*bvect(1) + bvect(2)*bvect(2) + bvect(3)*bvect(3))
  READ(filein,*) cvect(1), cvect(2), cvect(3)
  c = scale*sqrt(cvect(1)*cvect(1) + cvect(2)*cvect(2) + cvect(3)*cvect(3))
  alpha = ACOS((bvect(1)*cvect(1) + bvect(2)*cvect(2) + bvect(3)*cvect(3))/(b*c))
  beta = ACOS((avect(1)*cvect(1) + avect(2)*cvect(2) + avect(3)*cvect(3))/(a*c))
  gamma = ACOS((avect(1)*bvect(1) + avect(2)*bvect(2) + avect(3)*bvect(3))/(a*b))
  vol = a*b*c*   & 
      sqrt(1 - COS(alpha)*COS(alpha) - COS(beta)*COS(beta) - COS(gamma)*COS(gamma) - 2*COS(alpha)*COS(beta)*COS(gamma))

  READ(filein,'(A255)') line
  ! subroutine to get element name and number of kinds of element
  CALL getelemtnm(line,elemnm,ntype) 

  ! subrooutine to transfer string to number
  READ(filein,'(A255)') line
  WRITE(16,"(I6,7F16.8)") i, a, b, c, alpha/pi*180, beta/pi*180, gamma/pi*180, vol
  CALL getatmnr(line, maxatom, atnr) 

  READ(filein,*)

  !read coordinates and do manipulation
  100   WRITE(26,"(I6)") maxatom
  WRITE(26,*) "Conifiguration", i
  nr = 0 
  ALLOCATE(disp(ntype))
 
  DO k = 1, ntype
    disp(k) = 0
    DO m = 1, atnr(k)
      read(filein,*) x, y, z
      if(x<0) x=x+1; if(x>1) x=x-1
      if(y<0) y=y+1; if(y>1) y=y-1
      if(z<0) z=z+1; if(z>1) z=z-1
 !     xt = x*avect(1) + y*bvect(1) + z*cvect(1)
 !     yt = x*avect(2) + y*bvect(2) + z*cvect(2)
 !     zt = x*avect(3) + y*bvect(3) + z*cvect(3)
 !     x = xt * scale
 !     y = yt * scale
 !     z = zt * scale

      WRITE(26,"(A10,3F16.8)") elemnm(k),(x*avect(1)+y*bvect(1)+z*cvect(1))*scale, &
      (x*avect(2)+y*bvect(2)+z*cvect(2))*scale,(x*avect(3)+y*bvect(3)+z*cvect(3))*scale

      nr = nr + 1
      IF(i==1) THEN   ! store initial positions
        IF(nr==1) THEN
          ALLOCATE(xarray(maxatom))
          ALLOCATE(yarray(maxatom))
          ALLOCATE(zarray(maxatom))
        END IF
        xarray(nr) = x
        yarray(nr) = y
        zarray(nr) = z
      ELSE            ! calc deviations
        IF((x - xarray(nr)) > 0.5) Then
           x = x - xarray(nr) - 1
        ELSE IF((x - xarray(nr)) < -0.5) THEN
           x = (x - xarray(nr)) + 1 
        ELSE
           x = x - xarray(nr)
        END IF
        IF((y - yarray(nr)) > 0.5) then
           y =y - yarray(nr)-1
        ELSE IF((y - yarray(nr)) < -0.5) THEN
           y = y - yarray(nr) +1
        ELSE
           y = y - yarray(nr)
        END IF
        IF((z - zarray(nr)) > 0.5) then
           z = z - zarray(nr)- 1
        ELSE IF((z - zarray(nr)) < -0.5) THEN
           z = z - zarray(nr) +1
        ELSE
           z = z - zarray(nr)
        END IF

       disp(k) = disp(k) + (x*avect(1)+y*bvect(1)+z*cvect(1))**2 + &
                 (x*avect(2)+y*bvect(2)+z*cvect(2))**2 + &
                 (x*avect(3)+y*bvect(3)+z*cvect(3))**2
      END IF

    END DO

    disp(k) = sqrt(disp(k) / atnr(k))

  END DO
  WRITE(36,*) i, disp
  DEALLOCATE(disp)

END DO

PRINT *, "RMSD done!"
PRINT *, "Movie done!"
PRINT *, "Lattice done!"

DEALLOCATE(xarray)
DEALLOCATE(yarray)
DEALLOCATE(zarray)
CLOSE(filein)
CLOSE(16)
CLOSE(26)
CLOSE(36)

!-----------------------------------------------------------------------
! Calculate coord number 
!-----------------------------------------------------------------------
WRITE(*,*) "Analysis? y/n"
READ(*,*) line
IF(line=='y' .or. line=='yes') THEN
  bondtype = CHOOSE(ntype,2,IOstatus) + ntype   ! num of type of bond
  IF(IOstatus .ne. 0) STOP
  ALLOCATE(xarray(maxatom))
  ALLOCATE(yarray(maxatom))
  ALLOCATE(zarray(maxatom))
  ALLOCATE(xarrayimag(27,maxatom))
  ALLOCATE(yarrayimag(27,maxatom))
  ALLOCATE(zarrayimag(27,maxatom))
  ALLOCATE(cutoff(ntype,ntype,2))
  ALLOCATE(coordnum(maxatom,ntype,ntype))
  n = 0      
  DO k=1, ntype    ! read bond length and tolerance criteria
    DO m=k, ntype
      n=n+1
      WRITE(*,*) "Input the bond length and tolerance for: ",elemnm(k),"--",elemnm(m)
      READ(*,*) cutoff(k,m,:)   ! the correct tolerance should be the RMSD value
      IF(k.ne.m) cutoff(m,k,:)=cutoff(k,m,:)
      WRITE(fnm,'(A2,A2,A4)') elemnm(k),elemnm(m),".dat"
      OPEN(45+n,file=fnm, status='new')
    END DO
  END DO
ELSE 
  STOP
END IF

OPEN(26,file='xdat.xyz', status='old')
! loop over all (i-1) steps
DO m=1, i-1
  READ(26,*)
  READ(26,*)
  n = 0
  ! read coordinates and creat image coordinates
  DO k=1, ntype
    DO j=1, atnr(k)
      n = n + 1
      READ(26,*) line, xarray(n), yarray(n), zarray(n)
     !a-axial positive translated image
      xarrayimag(1,n) = xarray(n) + avect(1) * scale
      yarrayimag(1,n) = yarray(n) + avect(2) * scale
      zarrayimag(1,n) = zarray(n) + avect(3) * scale
     !a-axial negative translated image
      xarrayimag(2,n) = xarray(n) - avect(1) * scale
      yarrayimag(2,n) = yarray(n) - avect(2) * scale
      zarrayimag(2,n) = zarray(n) - avect(3) * scale
     !b-axial positive translated image
      xarrayimag(3,n) = xarray(n) + bvect(1) * scale
      yarrayimag(3,n) = yarray(n) + bvect(2) * scale
      zarrayimag(3,n) = zarray(n) + bvect(3) * scale
     !b-axial negative translated image
      xarrayimag(4,n) = xarray(n) - bvect(1) * scale
      yarrayimag(4,n) = yarray(n) - bvect(2) * scale
      zarrayimag(4,n) = zarray(n) - bvect(3) * scale
     !c-axial positive translated image
      xarrayimag(5,n) = xarray(n) + cvect(1) * scale
      yarrayimag(5,n) = yarray(n) + cvect(2) * scale
      zarrayimag(5,n) = zarray(n) + cvect(3) * scale
     !c-axial negative translated image
      xarrayimag(6,n) = xarray(n) - cvect(1) * scale
      yarrayimag(6,n) = yarray(n) - cvect(2) * scale
      zarrayimag(6,n) = zarray(n) - cvect(3) * scale
     ! a+b
      xarrayimag(7,n) = xarray(n) + ( avect(1) + bvect(1) ) * scale
      yarrayimag(7,n) = yarray(n) + ( avect(2) + bvect(2) ) * scale
      zarrayimag(7,n) = zarray(n) + ( avect(3) + bvect(3) ) * scale
      xarrayimag(8,n) = xarray(n) + ( avect(1) - bvect(1) ) * scale
      yarrayimag(8,n) = yarray(n) + ( avect(2) - bvect(2) ) * scale
      zarrayimag(8,n) = zarray(n) + ( avect(3) - bvect(3) ) * scale
      xarrayimag(9,n) = xarray(n) - ( avect(1) + bvect(1) ) * scale
      yarrayimag(9,n) = yarray(n) - ( avect(2) + bvect(2) ) * scale
      zarrayimag(9,n) = zarray(n) - ( avect(3) + bvect(3) ) * scale
      xarrayimag(10,n) = xarray(n) - ( avect(1) - bvect(1) ) * scale
      yarrayimag(10,n) = yarray(n) - ( avect(2) - bvect(2) ) * scale
      zarrayimag(10,n) = zarray(n) - ( avect(3) - bvect(3) ) * scale
     ! a+c
      xarrayimag(11,n) = xarray(n) + ( avect(1) + cvect(1) ) * scale 
      yarrayimag(11,n) = yarray(n) + ( avect(2) + cvect(2) ) * scale
      zarrayimag(11,n) = zarray(n) + ( avect(3) + cvect(3) ) * scale
      xarrayimag(12,n) = xarray(n) + ( avect(1) - cvect(1) ) * scale
      yarrayimag(12,n) = yarray(n) + ( avect(2) - cvect(2) ) * scale
      zarrayimag(12,n) = zarray(n) + ( avect(3) - cvect(3) ) * scale
      xarrayimag(13,n) = xarray(n) - ( avect(1) + cvect(1) ) * scale
      yarrayimag(13,n) = yarray(n) - ( avect(2) + cvect(2) ) * scale
      zarrayimag(13,n) = zarray(n) - ( avect(3) + cvect(3) ) * scale
      xarrayimag(14,n) = xarray(n) - ( avect(1) - cvect(1) ) * scale
      yarrayimag(14,n) = yarray(n) - ( avect(2) - cvect(2) ) * scale
      zarrayimag(14,n) = zarray(n) - ( avect(3) - cvect(3) ) * scale
     ! b+c
      xarrayimag(15,n) = xarray(n) + ( bvect(1) + cvect(1) ) * scale
      yarrayimag(15,n) = yarray(n) + ( bvect(2) + cvect(2) ) * scale
      zarrayimag(15,n) = zarray(n) + ( bvect(3) + cvect(3) ) * scale
      xarrayimag(16,n) = xarray(n) + ( bvect(1) - cvect(1) ) * scale
      yarrayimag(16,n) = yarray(n) + ( bvect(2) - cvect(2) ) * scale
      zarrayimag(16,n) = zarray(n) + ( bvect(3) - cvect(3) ) * scale
      xarrayimag(17,n) = xarray(n) - ( bvect(1) + cvect(1) ) * scale
      yarrayimag(17,n) = yarray(n) - ( bvect(2) + cvect(2) ) * scale
      zarrayimag(17,n) = zarray(n) - ( bvect(3) + cvect(3) ) * scale
      xarrayimag(18,n) = xarray(n) - ( bvect(1) - cvect(1) ) * scale
      yarrayimag(18,n) = yarray(n) - ( bvect(2) - cvect(2) ) * scale
      zarrayimag(18,n) = zarray(n) - ( bvect(3) - cvect(3) ) * scale
     ! a+b+c
      xarrayimag(19,n) = xarray(n) + ( avect(1) + bvect(1) + cvect(1) ) * scale
      yarrayimag(19,n) = yarray(n) + ( avect(2) + bvect(2) + cvect(2) ) * scale
      zarrayimag(19,n) = zarray(n) + ( avect(3) + bvect(3) + cvect(3) ) * scale
      xarrayimag(20,n) = xarray(n) + ( avect(1) + bvect(1) - cvect(1) ) * scale
      yarrayimag(20,n) = yarray(n) + ( avect(2) + bvect(2) - cvect(2) ) * scale
      zarrayimag(20,n) = zarray(n) + ( avect(3) + bvect(3) - cvect(3) ) * scale
      xarrayimag(21,n) = xarray(n) + ( avect(1) - bvect(1) + cvect(1) ) * scale
      yarrayimag(21,n) = yarray(n) + ( avect(2) - bvect(2) + cvect(2) ) * scale
      zarrayimag(21,n) = zarray(n) + ( avect(3) - bvect(3) + cvect(3) ) * scale
      xarrayimag(22,n) = xarray(n) + ( avect(1) - bvect(1) - cvect(1) ) * scale
      yarrayimag(22,n) = yarray(n) + ( avect(2) - bvect(2) - cvect(2) ) * scale
      zarrayimag(22,n) = zarray(n) + ( avect(3) - bvect(3) - cvect(3) ) * scale
      xarrayimag(23,n) = xarray(n) - ( avect(1) + bvect(1) + cvect(1) ) * scale
      yarrayimag(23,n) = yarray(n) - ( avect(2) + bvect(2) + cvect(2) ) * scale
      zarrayimag(23,n) = zarray(n) - ( avect(3) + bvect(3) + cvect(3) ) * scale
      xarrayimag(24,n) = xarray(n) - ( avect(1) + bvect(1) - cvect(1) ) * scale
      yarrayimag(24,n) = yarray(n) - ( avect(2) + bvect(2) - cvect(2) ) * scale
      zarrayimag(24,n) = zarray(n) - ( avect(3) + bvect(3) - cvect(3) ) * scale
      xarrayimag(25,n) = xarray(n) - ( avect(1) - bvect(1) + cvect(1) ) * scale
      yarrayimag(25,n) = yarray(n) - ( avect(2) - bvect(2) + cvect(2) ) * scale
      zarrayimag(25,n) = zarray(n) - ( avect(3) - bvect(3) + cvect(3) ) * scale
      xarrayimag(26,n) = xarray(n) - ( avect(1) - bvect(1) - cvect(1) ) * scale
      yarrayimag(26,n) = yarray(n) - ( avect(2) - bvect(2) - cvect(2) ) * scale
      zarrayimag(26,n) = zarray(n) - ( avect(3) - bvect(3) - cvect(3) ) * scale

      xarrayimag(27,n) = xarray(n)
      yarrayimag(27,n) = yarray(n)
      zarrayimag(27,n) = zarray(n)
    END DO
  END DO
  ! calculate bonding length for each pair, calcullate coord num for each kind of bond for each atom
  coordnum = 0
  n = 0   ! lable atom in lattice
  DO k=1, ntype   ! loop over elements 
    DO j=1, atnr(k)  ! loop over each element
      n = n + 1
      DO p=1, 27   ! loop over lattice images
        t = 0 ! label atom in image lattice
        DO q=1,ntype  ! loop over elements in image
          DO s=1, atnr(k) ! loop over each element
            t = t + 1
            bondlength = 0
            bondlength = (xarray(n)-xarrayimag(p,t))**2 + (yarray(n)-yarrayimag(p,t))**2 + (zarray(n)-zarrayimag(p,t))**2
            bondlength = SQRT(bondlength)
            IF( ABS(bondlength - cutoff(k,q,1)) .le. cutoff(k,q,2) )  coordnum(n,k,q) = coordnum(n,k,q) + 1
          END DO
        END DO     
      END DO
    END DO
  END DO

  ! write into filec
  n = 0   ! label bondtype
  DO k=1, ntype  
    DO j=k, ntype
      n=n+1
      WRITE(45+n,*) coordnum(:,k,j)
    END DO
  END DO

! end of loop of steps
END DO

DEALLOCATE(xarray)
DEALLOCATE(yarray)
DEALLOCATE(zarray)
DEALLOCATE(xarrayimag)
DEALLOCATE(yarrayimag)
DEALLOCATE(zarrayimag)
DEALLOCATE(cutoff)
DEALLOCATE(coordnum)
DO i=1, bondtype
  CLOSE(45+i)
END DO

PRINT *, "Analysis done!"
!-----------------------------------------------------------------------
! End of analysis  part

STOP

!-----------------------------------------------------------------------
! subroutines
!-----------------------------------------------------------------------
CONTAINS
 ! get the number of atoms
 SUBROUTINE getatmnr(line,maxatom,atnr)
 CHARACTER(len=255),intent(IN):: line
 INTEGER,intent(OUT) ::maxatom
 INTEGER :: j, ios,nratoms,atnrpos
 INTEGER,DIMENSION(:),ALLOCATABLE,intent(OUT) :: atnr
 INTEGER,DIMENSION(:),ALLOCATABLE :: atnrtmp
 CHARACTER(len=255):: line2,atl

 ios=0
 nratoms=0
 line2=line
 j = 0
 ALLOCATE(atnrtmp(20))

 DO WHILE (ios==0)
    j = j + 1
    read(line2,*,IOSTAT=ios) atnrtmp(j)
    IF (ios==0) then
       nratoms=nratoms+atnrtmp(j)
       write(atl,*) atnrtmp(j)
       atl=ADJUSTL(atl)
       atnrpos=Index(line2,trim(atl))+len_trim(atl)
! Index finds the location of a substring in another string,returns 0 if not found
! len_trim counts the number of characters including blanks
! this sentence locates the point in the line
       line2=""
       line2(atnrpos:255)=line(atnrpos:255)
       IF (len_trim(line2)<=0) then
          ios=10
       END IF
    END IF
 END DO
 maxatom=nratoms
 ALLOCATE(atnr(j))
 DO WHILE (j>0)
   atnr(j) = atnrtmp(j)   
   j = j - 1
 END DO
 
 END SUBROUTINE getatmnr

 ! get the element name
 SUBROUTINE getelemtnm(line,elemnm,ntype)
 CHARACTER(len=255),intent(IN):: line
 CHARACTER(len=255):: line2,line3
 CHARACTER(len=2),dimension(:),ALLOCATABLE, intent(OUT) :: elemnm
 CHARACTER(len=2),dimension(:),ALLOCATABLE :: elemt
 INTEGER, intent(OUT) :: ntype 
 INTEGER :: ios, i, pos

 line2=line
 ios=0
 i = 0
 ALLOCATE(elemt(20))

 DO WHILE (ios==0)
   i = i + 1
   READ(line2,*,IOSTAT=ios) elemt(i)
   IF (ios==0) then
     elemt(i) = ADJUSTL(elemt(i))
     pos = Index(line2,trim(elemt(i))) + len_trim(elemt(i))
     line2=""
     line2(pos:255)=line(pos:255)
     IF (len_trim(line2)<=0) then
       ios=10
     END IF
   END IF
 END DO
 ntype = i 
 ALLOCATE(elemnm(ntype))
 DO WHILE (i>0)
    elemnm(i) = elemt(i)
    i = i - 1
 END DO
 END SUBROUTINE getelemtnm

 ! calculate number of combinations by choosing k from n
 FUNCTION choose(n, k, err)
    integer :: choose
    integer, intent(in) :: n, k
    integer, optional, intent(out) :: err
    integer :: imax, i, imin, ie
    ie = 0
    if ( (n < 0 ) .or. (k < 0 ) ) then
       write(*, *) "negative in choose"
       choose = 0
       ie = 1
    else
       if ( n < k ) then
          choose = 0
       else if ( n == k ) then
          choose = 1
       else
          imax = max(k, n-k)
          imin = min(k, n-k)
          choose = 1
          do i = imax+1, n
             choose = choose * i
          end do
          do i = 2, imin
             choose = choose / i
          end do
       end if
    end if
    if ( present(err) ) err = ie
  END FUNCTION choose

END PROGRAM XDAT
