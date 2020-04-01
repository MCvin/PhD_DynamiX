program FORMATION

implicit none

!Systemes de reference des satellites
!tous les systemes sont exprimes dans le repere J2000 avec comme origine 
!le centre des miroirs
        type SystemeReference
           real (kind=8), dimension(3) :: position
           real (kind=8), dimension(3) :: attitude
           real (kind=8), dimension(3) :: positionprecedente
           real (kind=8), dimension(3) :: attitudeprecedente
        endtype SystemeReference


!------------------------------------------------------------------------
!-                      Programme principal                             -
!------------------------------------------------------------------------

! IN/OUT variables
        real (kind=8), dimension(2) :: source_J2000
        real (kind=8) :: focale
        integer :: duree

! local variables
        integer :: i, tempsmax
        type (SystemeReference) :: RefPM        !referentiel de pointage du satellite miroir
        type (SystemeReference) :: RefPD        !referentiel de pointage du satellite detecteur
        type (SystemeReference) :: RefSM        !referentiel du satellite miroir (centre de masse)
        type (SystemeReference) :: RefSD        !referentiel du satellite detecteur (centre de masse)

!variables pour la lecture
        character(len=80) :: chaine, mot
        integer :: err, debut, fin

!initialisation des parametres par lecture du fichier parametres.par
    open(1,file='../parametres.par')

    err = 0
    do while (err.eq.0)

       read(1,'(A)',iostat=err) chaine

       if ((chaine(1:1).ne.'#') .and. (err.eq.0)) then
          mot=chaine(1:(index(chaine,':'))-1)
          debut=index(chaine,'[')+1
          fin=index(chaine,']')-1

          select case (mot)
            !lecture des parametres de l'observation
             case ('focale')
                read(chaine(debut:fin),*) focale
          end select

       end if

    end do
    close(1)

    !source=(/ 0.0 , 0.0 /) pour exprimer les mouvements d'attitude en relatif
    source_J2000=(/ 0.0 , 0.0 /)
    !pour ponderer les mouvements
    duree=1
    !nombre d'etats de la formation
    tempsmax=86165


open(12,file='../input/mouvement_SM.data')
open(13,file='../input/mouvement_SD.data')

! Initialisation de la formation si necessaire (a t=0 la duree=0)
        call pointage(source_J2000,focale, RefPM,RefSM,RefPD,RefSD)

    do i=1,tempsmax
! Sauvegarde des positions des deux satellites
        write(12,*) RefSM%position, RefSM%attitude
        write(13,*) RefSD%position, RefSD%attitude

! Derive de la formation
        call derive(duree, RefSM,RefSD)

! Correction de la formation si necessaire
        call correction(RefPM,RefPD, RefSM,RefSD)
    enddo

close(12)
close(13)

contains

!------------------------------------------------------------------------
!-               Initialise la Formation pour une observation           -
!------------------------------------------------------------------------
subroutine pointage(source_J2000,focale, RefPM,RefSM,RefPD,RefSD)

implicit none

! IN/OUT variables
        real (kind=8), dimension(2) :: source_J2000
        real (kind=8) :: focale

        type (SystemeReference) :: RefPM        !referentiel de pointage du satellite miroir
        type (SystemeReference) :: RefSM        !referentiel du satellite miroir (centre de masse)
        type (SystemeReference) :: RefPD        !referentiel de pointage du satellite detecteur
        type (SystemeReference) :: RefSD        !referentiel du satellite detecteur (centre de masse)

! local variables
        real (kind=8), dimension(3) :: erratt_M, erratt_D
        real (kind=8), dimension(3) :: errpos_D


! generation d'erreurs sur les senseurs
        call random_number(erratt_M(1))
        call random_number(erratt_M(2))
        call random_number(erratt_M(3))
        call random_number(erratt_D(1))
        call random_number(erratt_D(2))
        call random_number(erratt_D(3))
        !erreur de l'ordre du millieme de degre sur le senseur stellaire (~arcsec)
        erratt_M=(erratt_M-0.5)/1000.0          !pour centrer sur 0.0
        erratt_D=(erratt_D-0.5)/1000.0          !pour centrer sur 0.0
        call random_number(errpos_D(1))
        call random_number(errpos_D(2))
        call random_number(errpos_D(3))
        !erreur de l'ordre du dixieme de millimetre sur les senseurs RF et optique
        errpos_D=(errpos_D-0.5)/10000.0         !pour centrer sur 0.0

! Initialise le systeme de reference lie au centre de masse du satellite miroir
        RefPM%position=(/ 0.0 , 0.0 , 0.0 /)
        RefPM%attitude(1)=-source_J2000(1)      !-alpha car alpha est dans le sens antitrigo
        RefPM%attitude(2)=source_J2000(2)       !delta
        RefPM%attitude(3)=0.0                   !0 car on suppose une source isotrope

        RefSM%positionprecedente=RefPM%position
        RefSM%attitudeprecedente=RefPM%attitude
        RefSM%position=RefPM%position
        !auquel on ajoute quelques erreurs
        RefSM%attitude=RefPM%attitude+erratt_M

! Initialise le systeme de reference lie au centre de masse du satellite detecteur
        RefPD%position=(/ 0.0 , 0.0 , 0.0 /)
        RefPD%position(3)=-focale-0.396         !-0.5 pour compenser les detecteurs qui sont a +0.5 en amont du plan focal
        RefPD%attitude(1)=-source_J2000(1)      !-alpha car alpha est dans le sens antitrigo
        RefPD%attitude(2)=source_J2000(2)       !delta
        RefPD%attitude(3)=0.0                   !0 car on suppose une source isotrope

        RefSD%positionprecedente=RefPD%position
        RefSD%attitudeprecedente=RefPD%attitude
        !auquel on ajoute quelques erreurs
        RefSD%position=RefPD%position+errpos_D
        RefSD%attitude=RefPD%attitude+erratt_D
        
print*, "position et attitude initiale du satellite Miroir (X,Y,Z):"
print*, RefSM%position
print*, RefSM%attitude
print*, ""
print*, "position et attitude initiale du satellite Detecteur (X,Y,Z):"
print*, RefSD%position
print*, RefSD%attitude
print*, ""

end subroutine pointage


!------------------------------------------------------------------------
!-               Simule la derive de la formation                       -
!------------------------------------------------------------------------
subroutine derive(duree, RefSM,RefSD)

implicit none

! IN/OUT variables
        integer :: duree  !entre les deux appels de procedure

        type (SystemeReference) :: RefSM
        type (SystemeReference) :: RefSD

! local variables
        real (kind=8), dimension(3) :: tampon
        real (kind=8), dimension(3) :: accelattSM, accelattSD, accelposSD
        real (kind=8) :: rotation


!if (duree == 0) then
   call random_number(accelattSM(1))
   call random_number(accelattSM(2))
   call random_number(accelattSM(3))
   accelattSM=(accelattSM-0.5)/10000.0

   call random_number(accelattSD(1))
   call random_number(accelattSD(2))
   call random_number(accelattSD(3))
   accelattSD=(accelattSD-0.5)/10000.0

   call random_number(accelposSD(1))
   call random_number(accelposSD(2))
   call random_number(accelposSD(3))
   accelposSD=(accelposSD-0.5)/100000.0
!else
   !rotation=360.0/84600.0
   !call rotationY(accelattSM,rotation) !pour simuler le vent solaire variant sur l'orbite
   !call rotationY(accelattSD,rotation)
   !call rotationY(accelposSD,rotation)
!endif

if (duree /= 0) then

   ! Derive de l'attitude
        !du satellite miroir par rapport a l'attitude de reference: RefPM%attitude
                !par la pression de radiation solaire (acceleration arbitraire: accelattSM)
        tampon=RefSM%attitude
        RefSM%attitude(1)=RefSM%attitude(1)+((RefSM%attitude(1)-RefSM%attitudeprecedente(1))+accelattSM(1))*duree
        RefSM%attitude(2)=RefSM%attitude(2)+((RefSM%attitude(2)-RefSM%attitudeprecedente(2))+accelattSM(2))*duree
        RefSM%attitude(3)=RefSM%attitude(3)+((RefSM%attitude(3)-RefSM%attitudeprecedente(3))+accelattSM(3))*duree
        RefSM%attitudeprecedente=tampon

        !du satellite detecteur par rapport a l'attitude de reference: RefPD%attitude
                !par la pression de radiation solaire (acceleration arbitraire: accelattSD)
        tampon=RefSD%attitude
        RefSD%attitude(1)=RefSD%attitude(1)+((RefSD%attitude(1)-RefSD%attitudeprecedente(1))+accelattSD(1))*duree
        RefSD%attitude(2)=RefSD%attitude(2)+((RefSD%attitude(2)-RefSD%attitudeprecedente(2))+accelattSD(2))*duree
        RefSD%attitude(3)=RefSD%attitude(3)+((RefSD%attitude(3)-RefSD%attitudeprecedente(3))+accelattSD(3))*duree
        RefSD%attitudeprecedente=tampon

   ! Derive de la position
        !du satellite detecteur par rapport a la position de reference: RefPD%position
                !par la pression de radiation solaire (acceleration arbitraire: accelposSD)
        tampon=RefSD%position
        RefSD%position(1)=RefSD%position(1)+((RefSD%position(1)-RefSD%positionprecedente(1))+accelposSD(1))*duree
        RefSD%position(2)=RefSD%position(2)+((RefSD%position(2)-RefSD%positionprecedente(2))+accelposSD(2))*duree
        RefSD%position(3)=RefSD%position(3)+((RefSD%position(3)-RefSD%positionprecedente(3))+accelposSD(3))*duree
        RefSD%positionprecedente=tampon

!print*, RefSM%attitude
!print*, RefSD%position
!print*, RefSD%attitude

endif

end subroutine derive


!------------------------------------------------------------------------
!-             Correction de la derive de la formation                  -
!------------------------------------------------------------------------
subroutine correction(RefPM,RefPD, RefSM,RefSD)

implicit none

! IN/OUT variables
        type (SystemeReference) :: RefPM
        type (SystemeReference) :: RefPD

        type (SystemeReference) :: RefSM
        type (SystemeReference) :: RefSD

! local variables
        real (kind=8), dimension(3) :: ErrAttSM
        real (kind=8), dimension(3) :: ErrAttSD
        real (kind=8), dimension(3) :: ErrPosSD


! Correction si necessaire de l'attitude
        ErrAttSM=RefSM%attitude-RefPM%attitude
        ErrAttSD=RefSD%attitude-RefPD%attitude
        ErrPosSD=RefSD%position-RefPD%position

        ! du satellite miroir
        ! a completer avec les senseurs
        if (ABS(ErrAttSM(1)) > (20.0/3600.0)) then      !(20.0/3600.0) pour 20arcsec
           RefSM%attitudeprecedente(1)=RefSM%attitude(1)
           RefSM%attitudeprecedente(2)=RefSM%attitude(2)
           RefSM%attitudeprecedente(3)=RefSM%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction d'attitude
           RefSM%attitude(1)=RefSM%attitudeprecedente(1)-ErrAttSM(1)/100.0
           RefSM%attitude(2)=RefSM%attitudeprecedente(2)-ErrAttSM(2)/100.0
           RefSM%attitude(3)=RefSM%attitudeprecedente(3)-ErrAttSM(3)/100.0
        endif
        if (ABS(ErrAttSM(2)) > (20.0/3600.0)) then      !(20.0/3600.0) pour 20arcsec
           RefSM%attitudeprecedente(1)=RefSM%attitude(1)
           RefSM%attitudeprecedente(2)=RefSM%attitude(2)
           RefSM%attitudeprecedente(3)=RefSM%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction d'attitude
           RefSM%attitude(1)=RefSM%attitudeprecedente(1)-ErrAttSM(1)/100.0
           RefSM%attitude(2)=RefSM%attitudeprecedente(2)-ErrAttSM(2)/100.0
           RefSM%attitude(3)=RefSM%attitudeprecedente(3)-ErrAttSM(3)/100.0
        endif
        if (ABS(ErrAttSM(3)) > (1.0/60.0)) then         !(1.0/60.0) pour 1arcmin
           RefSM%attitudeprecedente(3)=RefSM%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction d'attitude
           RefSM%attitude(3)=RefSM%attitudeprecedente(3)-ErrAttSM(3)/100.0
        endif

        ! du satellite detecteur
        ! a completer avec les senseurs
        if (ABS(ErrAttSD(1)) > (2.0/60.0)) then         !(2.0/60.0) pour 2arcmin  
           RefSD%attitudeprecedente(1)=RefSD%attitude(1)
           RefSD%attitudeprecedente(2)=RefSD%attitude(2)
           RefSD%attitudeprecedente(3)=RefSD%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction d'attitude
           RefSD%attitude(1)=RefSD%attitudeprecedente(1)-ErrAttSD(1)/100.0
           RefSD%attitude(2)=RefSD%attitudeprecedente(2)-ErrAttSD(2)/100.0
           RefSD%attitude(3)=RefSD%attitudeprecedente(3)-ErrAttSD(3)/100.0
        endif
        if (ABS(ErrAttSD(2)) > (2.0/60.0)) then         !(2.0/60.0) pour 2arcmin
           RefSD%attitudeprecedente(1)=RefSD%attitude(1)
           RefSD%attitudeprecedente(2)=RefSD%attitude(2)
           RefSD%attitudeprecedente(3)=RefSD%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction d'attitude
           RefSD%attitude(1)=RefSD%attitudeprecedente(1)-ErrAttSD(1)/100.0
           RefSD%attitude(2)=RefSD%attitudeprecedente(2)-ErrAttSD(2)/100.0
           RefSD%attitude(3)=RefSD%attitudeprecedente(3)-ErrAttSD(3)/100.0
        endif
        if (ABS(ErrAttSD(3)) > (2.0/60.0)) then         !(2.0/60.0) pour 2arcmin
           RefSD%attitudeprecedente(3)=RefSD%attitude(3)
           ! /100.0 est un facteur arbitraire d'acceleration du a la correction d'attitude
           RefSD%attitude(3)=RefSD%attitudeprecedente(3)-ErrAttSD(3)/100.0
        endif

! Correction si necessaire de la position

        ! du satellite detecteur
        ! a completer avec les senseurs
        if (ABS(ErrPosSD(1)) > 0.005) then         !0.005 pour 0.5cm sur les X
           RefSD%positionprecedente(1)=RefSD%position(1)
           RefSD%positionprecedente(2)=RefSD%position(2)
           RefSD%positionprecedente(3)=RefSD%position(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction de position
           RefSD%position(1)=RefSD%positionprecedente(1)-ErrPosSD(1)/100.0
           RefSD%position(2)=RefSD%positionprecedente(2)-ErrPosSD(2)/100.0
           RefSD%position(3)=RefSD%positionprecedente(3)-ErrPosSD(3)/100.0
        endif
        if (ABS(ErrPosSD(2)) > 0.005) then         !0.005 pour 0.5cm sur les Y
           RefSD%positionprecedente(1)=RefSD%position(1)
           RefSD%positionprecedente(2)=RefSD%position(2)
           RefSD%positionprecedente(3)=RefSD%position(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction de position
           RefSD%position(1)=RefSD%positionprecedente(1)-ErrPosSD(1)/100.0
           RefSD%position(2)=RefSD%positionprecedente(2)-ErrPosSD(2)/100.0
           RefSD%position(3)=RefSD%positionprecedente(3)-ErrPosSD(3)/100.0
        endif
        if (ABS(ErrPosSD(3)) > 0.005) then         !0.005 pour 0.5cm sur les Z
           RefSD%positionprecedente(1)=RefSD%position(1)
           RefSD%positionprecedente(2)=RefSD%position(2)
           RefSD%positionprecedente(3)=RefSD%position(3)
           ! /100.0 est un facteur arbitraire d'acceleration suite a la manoeuvre de correction de position
           RefSD%position(1)=RefSD%positionprecedente(1)-ErrPosSD(1)/100.0
           RefSD%position(2)=RefSD%positionprecedente(2)-ErrPosSD(2)/100.0
           RefSD%position(3)=RefSD%positionprecedente(3)-ErrPosSD(3)/100.0
        endif

end subroutine correction


end program FORMATION
