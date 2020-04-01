program SIMU

use SENSEURS
use TELESCOPE
use PHOTON
use FPA

implicit none

!parametres observation
        type (Observation) :: Obs
        type (photons) :: photon_D

!parametres satellite miroirs
        type (ParametresSatMiroir) :: Para_SM
        type (Repere) :: RefM

!parametre satellite detecteurs
        type (ParametresSatDetecteur) :: Para_SD
        type (Repere) :: RefD

!parametres internes
        integer :: i, j, etat, n
        real (kind=wp) :: k, l
        real (kind=wp), dimension(:,:), allocatable :: PSF_data
        real (kind=wp), dimension(2) :: centroide
        real (kind=wp), dimension(:), allocatable :: r
        integer :: IFAIL
        EXTERNAL M01CAF
        real (kind=wp) :: HEW, HEW_corr, Pi, focale
        character(8) :: date
        character(10) :: time
        character(26) :: directory
        character(50) :: command


!------------------------------------------------------------------------
!-                      Debut du programme                              -
!------------------------------------------------------------------------
        !effacement des fichiers output precedents
        command='rm output/*.data'
        call system(command)

j=0
k=0
l=0
HEW=0
HEW_corr=0

        !lecture et initialisation des parametres
        call lectureparametres(Obs,Para_SM,Para_SD)

print*
print*, 'Avancement:'

        !on parcours l'ensemble des etats de la formation
        do etat=1,Para_SM%nb_etats
           !affiche l'avancement du calcul en pourcentage
           if (mod(etat,Para_SM%nb_etats/10) .eq. 0) then
              print*, etat*100/Para_SM%nb_etats, '%'
           endif

           !calcul de l'attitude et de la position des reperes instruments    
           call etat_telescope(Para_SM,Para_SD,etat, RefM,RefD)

           !simule les donnees senseurs
           call donnees_senseurs_interpretees(Para_SM,Para_SD,etat)
        
           do n=1,Obs%nb_sources
              if (Obs%luminosites_sources) then
                 do i=1,Obs%nb_photons*Obs%sources(n)%luminosite/(Para_SM%nb_etats)

                    !effectue le transport d'un photon a travers le telescope        
                    call parcoursphoton(Obs,Para_SM,Para_SD,RefM,RefD,n, photon_D)
                    call plan_de_detection(Para_SM,Para_SD,photon_D, j,k,l)

                 enddo
              else
                 do i=1,Obs%nb_photons/(Para_SM%nb_etats*Obs%nb_sources)

                    !effectue le transport d'un photon a travers le telescope        
                    call parcoursphoton(Obs,Para_SM,Para_SD,RefM,RefD,n, photon_D)
                    call plan_de_detection(Para_SM,Para_SD,photon_D, j,k,l)

                 enddo
              endif
           enddo

        enddo

        focale=-Para_SD%pointage%position(3) + Para_SM%pos_miroirs(3) - Para_SD%pos_LED(3) - Para_SD%defocus

        if (j .ne. 0) then
        !caracteristiques de la PSF: HEW
                !lecture des coordonnees des photons composant la tache focale
                REWIND(11)
                allocate(PSF_data(3,j))
                read(11) PSF_data
                !calcul des coordonnees du centre de la tache
                centroide(1)=SUM(PSF_data(1,:))/(MAX(1,j))
                centroide(2)=SUM(PSF_data(2,:))/(MAX(1,j))
                !calcul de la distance de chaque point par rapport au centre de la tache
                allocate(r(j))
                r=SQRT( (PSF_data(1,:)-centroide(1))**2 + (PSF_data(2,:)-centroide(2))**2 )
                !rangement dans l'ordre croissant des distances
                IFAIL=0
                call M01CAF(r,1,j,'A',IFAIL)
                !calcul de la HEW
                Pi=ACOS(-1._wp)
                HEW=ATAN( (r(j/2)*2) / (focale*100) ) *180/Pi * 3600    !focale*100 car les coordonnees des photons sont en cm

        !caracteristiques de la PSF corrigee: HEW
                !lecture des coordonnees des photons composant la tache focale
                REWIND(12)
                read(12) PSF_data
                !calcul des coordonnees du centre de la tache
                centroide(1)=SUM(PSF_data(1,:))/(MAX(1,j))
                centroide(2)=SUM(PSF_data(2,:))/(MAX(1,j))
                !calcul de la distance de chaque point par rapport au centre de la tache
                r=SQRT( (PSF_data(1,:)-centroide(1))**2 + (PSF_data(2,:)-centroide(2))**2 )
                !rangement dans l'ordre croissant des distances
                IFAIL=0
                call M01CAF(r,1,j,'A',IFAIL)
                !calcul de la HEW
                Pi=ACOS(-1._wp)
                HEW_corr=ATAN( (r(j/2)*2) / (focale*100) ) *180/Pi * 3600    !focale*100 car les coordonnees des photons sont en cm
        endif

print*
print*, 'Nombre de photons focalises sur le plan de detection:',j
if (j .ne. 0) then
print*, 'Nombre de photons detectes par le LED:',k,'(',k/j*100,'%)'
print*, 'Nombre de photons detectes par le HED:',l,'(',l/j*100,'%)'
print*, 'Half Energy Width:',HEW,'arcsec'
print*, 'Half Energy Width corrected:',HEW_corr,'arcsec'
endif

write(10,*) j
write(10,*) focale
write(10,*) maxval(Para_SM%rayon)
write(10,*) HEW
write(10,*) HEW_corr
write(10,*) k
write(10,*) l
close(10)
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)

        !sauvegarde des donnees (la methode est un peu sombre...)
        if (Obs%saves) then
           call date_and_time(DATE=date,TIME=time)
           directory='saves/save_'//date//'-'//time
           call system('mkdir '//directory)
           command='cp switch.par '//directory
           call system(command)
           command='cp parametres.par '//directory
           call system(command)
           command='cp -r input '//directory
           call system(command)
           command='cp -r output '//directory
           call system(command)
        endif

        !fichiers temporaires
!!$        close(21)
!!$        close(22)
!!$        close(23)
!!$        close(24)
!!$        close(25)
        close(31)
        close(32)

print*
print*, 'Fin du programme'

end program SIMU
