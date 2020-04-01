!------------------------------------------------------------------------
!-         Monte Carlo                                                  -
!-         Plan de detection de Simbol-X                                -
!-         Detecteurs Si(LED) et Cd(Zn)Te(HED)                          -
!-         Interactions Compton + Photoelectrique                       -
!------------------------------------------------------------------------
module FPA

use OUTILS
use LECTURE

implicit none

    type Eabsorbee
       real (kind=wp), dimension(3) :: position
       integer, dimension(2) :: pixel
       real (kind=wp) :: energie
       logical :: Silicon, Cadmium, Tellurium, Photo, Compton
       logical :: photon_perdu
    endtype Eabsorbee
    type Epixel
       integer, dimension(2) :: pixel
       real (kind=wp) :: energie
    endtype Epixel

    type Detecteur
       logical :: actif                         !detecteur present ou non
       real (kind=wp), dimension(2) :: dim      !dimension x,y du detecteur (en cm)
       real (kind=wp) :: h_top                  !hauteur z du dessus du detecteur (en cm)
       real (kind=wp) :: h_bot                  !hauteur z du dessous du detecteur (en cm)
       real (kind=wp) :: dim_pixel              !dimension des pixels du detecteur (en cm)
       integer, dimension(135,135) :: map       !geometrie du detecteur (HED: zones actives et mortes)
       integer :: nb                            !nombre d'interactions par photon
       type (Eabsorbee), dimension(10) :: Eabs  !liste des energies absorbees, et position des interactions
       type (Epixel), dimension(10) :: Epix     !energies absorbees par pixel
       real (kind=wp), dimension(3) :: err_Pos  !erreur de positionnement du detecteur du a la formation

    endtype Detecteur


contains
!------------------------------------------------------------------------
!-                      Programme principal                             -
!-           ATTENTION: energie en MeV, distance en cm                  -
!------------------------------------------------------------------------
subroutine plan_de_detection(Para_SM,Para_SD,photon, j,k,l)
implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM           !parametres Satellite Miroir
        type (ParametresSatDetecteur) :: Para_SD        !parametres Satellite Detecteur
        type (photons) :: photon                        !position, direction et energie du photon

        integer :: j
        real (kind=wp) :: k, l

! local variables
        type (Detecteur) :: LED         !LED: nb d'interactions, liste des energies absorbees, ...
        type (Detecteur) :: HED         !HED: nb d'interactions, liste des energies absorbees, ...
        integer :: i
        logical :: entreecollimateur, surface_LED, surface_HED
        real (kind=wp) :: x, y
        integer :: f

entreecollimateur=.FALSE.
surface_LED=.FALSE.
surface_HED=.FALSE.

        if (Para_SM%reflexionhyperbole) then
           !===========================================INITIALISATION DES VARIABLES===========================================
           call init(Para_SD, LED,HED,photon)

           !coordonnees du photon sur la face d'entree du collimateur
           x=photon%position(1)+(photon%position(3)-Para_SD%h_col*100)*photon%direction(1)/abs(photon%direction(3))
           y=photon%position(2)+(photon%position(3)-Para_SD%h_col*100)*photon%direction(2)/abs(photon%direction(3))
           !le photon passe t-il dans le collimateur
           if ((x**2+y**2) .lt. (Para_SD%r_col*100)**2) then
              entreecollimateur=.TRUE.
           endif

        endif

        if (entreecollimateur) then
           if (LED%actif) then          !LED ou LED+HED
              !coordonnees du photon a la surface du LED
              photon%position=photon%position+(LED%h_top-photon%position(3))*photon%direction/photon%direction(3)
              if ((ABS(photon%position(1)) .le. LED%dim(1)/2) .and. (ABS(photon%position(2)) .le. LED%dim(2)/2)) then
                 !=>LE PHOTON EST A LA SURFACE D'UN PIXEL DU LED
                 surface_LED=.TRUE.
              endif
           else if (HED%actif) then     !HED
              !coordonnees du photon a la surface du HED
              photon%position=photon%position+(HED%h_top-photon%position(3))*photon%direction/photon%direction(3)
              if ((ABS(photon%position(1)) .le. HED%dim(1)/2) .and. (ABS(photon%position(2)) .le. HED%dim(2)/2)) then
                 !=>LE PHOTON EST A LA SURFACE D'UN PIXEL DU HED
                 surface_HED=.TRUE.
              endif
           endif
        endif

        if (surface_LED .or. surface_HED) then

                !=>ecriture des fichiers position reelles(cm) + position reelles corrigees (a la surface du LED ou du HED)
                j=j+1
                write(11) photon%position(1), photon%position(2), photon%energie*1000
                if (surface_LED) then
                   write(12) photon%position(1)+LED%err_Pos(1), photon%position(2)+LED%err_Pos(2), photon%energie*1000
                else
                   write(12) photon%position(1)+HED%err_Pos(1), photon%position(2)+HED%err_Pos(2), photon%energie*1000
                endif

           !====================================INTERACTION DU PHOTON AVEC LES DETECTEURS=====================================
           call detecteurs(LED,HED,photon)

           !============================================RESOLUTION DES DETECTEURS=============================================
           call resolution('LED',LED)
           call resolution('HED',HED)

                !=>ecriture des fichiers spectres du LED et du HED (sans pixels)
                if (SUM(LED%Epix%energie) .ne. 0) then
                   k=k+1
                   write(13,*) SUM(LED%Epix%energie)
                endif
                if (SUM(HED%Epix%energie) .ne. 0) then
                   l=l+1
                   write(14,*) SUM(HED%Epix%energie)
                endif
                !=>ecriture des fichiers images et spectres du LED et du HED (avec pixels)
                i=1
                do while (LED%Epix(i)%energie .ne. 0)
                   write(15,*) LED%Epix(i)%pixel, LED%Epix(i)%energie
                   i=i+1
                enddo
                i=1
                do while (HED%Epix(i)%energie .ne. 0)
                   write(16,*) HED%Epix(i)%pixel, HED%Epix(i)%energie
                   i=i+1
                enddo

           !=============================================CORRECTION DE L'IMAGERIE=============================================
           call correction(LED,HED)

                !=>ecriture des fichiers images corrigees et spectres du LED et du HED (avec pixels)
                i=1
                do while (LED%Epix(i)%energie .ne. 0)
                   if ((LED%Epix(i)%pixel(1) .ge. 1) .and. (LED%Epix(i)%pixel(2) .ge. 1) .and. &
                       (LED%Epix(i)%pixel(1) .lt. 128) .and. (LED%Epix(i)%pixel(2) .lt. 128)) then
                      write(17,*) LED%Epix(i)%pixel, LED%Epix(i)%energie
                   endif
                   i=i+1
                enddo
                i=1
                do while (HED%Epix(i)%energie .ne. 0)
                   if ((HED%Epix(i)%pixel(1) .ge. 1) .and. (HED%Epix(i)%pixel(2) .ge. 1) .and. &
                       (HED%Epix(i)%pixel(1) .lt. 135) .and. (HED%Epix(i)%pixel(2) .lt. 135)) then
                      write(18,*) HED%Epix(i)%pixel, HED%Epix(i)%energie
                   endif
                   i=i+1
                enddo

!!$           !=>pour produire des images en cours de focalisation
!!$           do f=197,0,-1
!!$              !coordonnees du photon tout le long de la focalisation
!!$              x=photon%position(1)+(photon%position(3)-f/10._wp)*photon%direction(1)/abs(photon%direction(3))
!!$              y=photon%position(2)+(photon%position(3)-f/10._wp)*photon%direction(2)/abs(photon%direction(3))
!!$              j=j+1
!!$              write(11) x, y
!!$           enddo
        endif

end subroutine plan_de_detection

!------------------------------------------------------------------------
!-         Initialisation des variables LED et HED                      -
!-                                                                      -
!------------------------------------------------------------------------
subroutine init(Para_SD, LED,HED,photon)
implicit none

! IN/OUT variables
        type (ParametresSatDetecteur) :: Para_SD!parametres Satellite Detecteur
        type (Detecteur) :: LED                 !LED: nb d'interactions, liste des energies absorbees, ...
        type (Detecteur) :: HED                 !HED: nb d'interactions, liste des energies absorbees, ...
        type (photons) :: photon                !position, direction et energie du photon


!=>Initialisation des variables de LED:
        LED%actif=Para_SD%LED
        !pour passer des metres aux cm
        LED%dim=Para_SD%dimensions_LED*100
        LED%h_top=0
        LED%h_bot=LED%h_top-0.045_wp     !450 micrometres d'epaisseur
        LED%dim_pixel=Para_SD%taille_pixels_LED*100
        LED%err_Pos=Para_SD%err_Pos_LED*100
        LED%nb=0
        LED%Eabs%position(1)=0
        LED%Eabs%position(2)=0
        LED%Eabs%pixel(1)=0
        LED%Eabs%pixel(2)=0
        LED%Eabs%energie=0
        LED%Eabs%Silicon=.FALSE.
        LED%Eabs%Cadmium=.FALSE.
        LED%Eabs%Tellurium=.FALSE.
        LED%Eabs%Photo=.FALSE.
        LED%Eabs%Compton=.FALSE.
        LED%Eabs%photon_perdu=.FALSE.
        LED%Epix%pixel(1)=0
        LED%Epix%pixel(2)=0
        LED%Epix%energie=0

!=>Initialisation des variables de HED:
        HED%actif=Para_SD%HED
        !pour passer des metres aux cm
        HED%dim=Para_SD%dimensions_HED*100
        HED%h_top=Para_SD%pos_HED(3)*100-Para_SD%pos_LED(3)*100+LED%h_bot
        HED%h_bot=HED%h_top-0.2_wp      !2 mm d'epaisseur
        HED%dim_pixel=Para_SD%taille_pixels_HED*100
        HED%map=Para_SD%HED_map
        HED%err_Pos=Para_SD%err_Pos_HED*100
        HED%nb=0
        HED%Eabs%position(1)=0
        HED%Eabs%position(2)=0
        HED%Eabs%pixel(1)=0
        HED%Eabs%pixel(2)=0
        HED%Eabs%energie=0
        HED%Eabs%Silicon=.FALSE.
        HED%Eabs%Cadmium=.FALSE.
        HED%Eabs%Tellurium=.FALSE.
        HED%Eabs%Photo=.FALSE.
        HED%Eabs%Compton=.FALSE.
        HED%Eabs%photon_perdu=.FALSE.
        HED%Epix%pixel(1)=0
        HED%Epix%pixel(2)=0
        HED%Epix%energie=0

!=>Conversion des unites pour le photon:
        !pour passer des metres aux cm
        photon%position=photon%position*100
        !pour passer des keV a MeV
        photon%energie=photon%energie/1000

end subroutine init

!------------------------------------------------------------------------
!-         Interactions du photon incident                              -
!-         avec les detecteurs (LED et HED)                             -
!------------------------------------------------------------------------
subroutine detecteurs(LED,HED,photon)
implicit none

! IN/OUT variables
        type (Detecteur) :: LED         !LED: nb d'interactions, liste des energies absorbees, ...
        type (Detecteur) :: HED         !HED: nb d'interactions, liste des energies absorbees, ...
        type (photons) :: photon        !position, direction et energie du photon

! local variables
        integer, dimension(2) :: pixel


        if (LED%actif .and. HED%actif) then

              call Al_filter(photon)
              !=>LE PHOTON N'A PAS ETE ABSORBE, IL EST A LA SURFACE DU LED
              call inside_LED(LED, photon)
              !print*,'LED',LED%Eabs%energie

              !si un photon se situe entre les detecteurs on le suit jusqu'au prochain (HED ou LED)
              do while ((photon%position(3) .le. LED%h_bot) .and. (photon%position(3) .ge. HED%h_top))

                 !=>LE PHOTON SE DIRIGE VERS LE HED
                 if (photon%direction(3) .lt. 0) then
                    !position du photon a la surface du HED:
                    photon%position=photon%position+(HED%h_top-photon%position(3))*photon%direction/photon%direction(3)
                    if ((ABS(photon%position(1)) .le. HED%dim(1)/2) .and. (ABS(photon%position(2)) .le. HED%dim(2)/2)) then
                       !=>LE PHOTON EST A LA SURFACE DU HED
                       call inside_HED(HED, photon)
                    else
                       !=>LE PHOTON N'EST PAS A LA SURFACE DU HED
                       photon%position=0   !pour sortir de la boucle
                    endif
                    !print*,'HED',HED%Eabs%energie

                 !=>LE PHOTON SE DIRIGE VERS LE LED
                 else if (photon%direction(3) .gt. 0) then
                    !position du photon a la surface du LED:
                    photon%position=photon%position+(LED%h_bot-photon%position(3))*photon%direction/photon%direction(3)
                    if ((ABS(photon%position(1)) .le. LED%dim(1)/2) .and. (ABS(photon%position(2)) .le. LED%dim(2)/2)) then
                       !=>LE PHOTON EST A LA SURFACE DU LED
                       call inside_LED(LED, photon)
                    else
                       !=>LE PHOTON N'EST PAS A LA SURFACE DU HED
                       photon%position=0   !pour sortir de la boucle
                    endif
                    !print*,'LED',LED%Eabs%energie

                 !=>LE PHOTON SE DIRIGE HORIZONTALEMENT
                 else
                    !le photon sort du plan de detection
                    photon%position=0   !pour sortir de la boucle
                 endif

              enddo

        else if (LED%actif) then

           call Al_filter(photon)
           !=>LE PHOTON N'A PAS ETE ABSORBE, IL EST A LA SURFACE DU LED
           call inside_LED(LED, photon)

        else if (HED%actif) then

           call inside_HED(HED, photon)

        endif

end subroutine detecteurs

!------------------------------------------------------------------------
!-         Resolution spatiale et energetique                           -
!-         du detecteur LED et HED                                      -
!------------------------------------------------------------------------
subroutine resolution(detector,ED)
implicit none

! IN/OUT variables
        character (len=3) :: detector                   !detecteur LED ou HED
        type (Detecteur) :: ED                          !nb d'interactions, liste des energies absorbees, ...

! Local variables
        integer :: i, n, indice
        integer, dimension(2) :: min, max, pixel        !coordonnees du pixel ou a eu lieu l'interaction
        integer, dimension(:,:), allocatable :: tab     !tab representant une partie du LED ou du HED et contenant l'indice correspondant
        real (kind=wp) :: G05DDF
        real (kind=wp) :: alea
        real (kind=wp) :: E


        if (SUM(ED%Eabs%energie) .ne. 0) then
           min(1)=MINVAL(ED%Eabs%pixel(1),mask=ED%Eabs%pixel(1) .ne. 0)
           min(2)=MINVAL(ED%Eabs%pixel(2),mask=ED%Eabs%pixel(2) .ne. 0)
           max(1)=MAXVAL(ED%Eabs%pixel(1))
           max(2)=MAXVAL(ED%Eabs%pixel(2))
           !on alloue un tableau representant la partie du detecteur touchee
           allocate(tab(max(1)-min(1)+1,max(2)-min(2)+1))
           tab=0

           n=0
           if (detector .eq. 'LED') then
              do i=1,ED%nb
                 !=>RESOLUTION EN ENERGIE: on pondere l'energie mesuree par une incertitude
                 alea=G05DDF(0._wp,0.000064_wp)                    !resolution du LED: 150 eV FWHM = 0.15keV/2.35 sigma
                 E=ED%Eabs(i)%energie+alea
                 !=>SUEIL DE DETECTION
                 if(E .ge. 0.0001_wp) then                         !seuil de detection du LED: 0.1 keV
                    !=>RESOLUTION SPATIALE: on somme les energies par pixels
                    !dans quel pixel de tab l'energie s'est elle deposee?
                    pixel=ED%Eabs(i)%pixel-min+1
                    if (tab(pixel(1),pixel(2)) .eq. 0) then        !on est en presence d'un nouveau pixel
                       n=n+1
                       tab(pixel(1),pixel(2))=n                    !on lui attribue un nouvel indice que l'on stock dans tab_pixels
                       ED%Epix(n)%pixel=ED%Eabs(i)%pixel           !on remplit la case correspondante avec sa position pixel
                       ED%Epix(n)%energie=E                        !on remplit la case correspondante avec son energie
                    else
                       indice=tab(pixel(1),pixel(2))               !on va chercher l'indice correspondant au pixel
                       ED%Epix(indice)%energie=ED%Epix(indice)%energie + E !on remplit la case correspondante avec son energie
                    endif
                 endif
              enddo

           else if (detector .eq. 'HED') then
              do i=1,ED%nb
                 !=>RESOLUTION EN ENERGIE: on pondere l'energie mesuree par une incertitude
                 alea=G05DDF(0._wp,0.00042_wp)                     !resolution du HED: 1 keV FWHM = 1keV/2.35 sigma
                 E=ED%Eabs(i)%energie+alea
                 !=>SUEIL DE DETECTION
                 if(E .ge. 0.0015_wp) then                         !seuil de detection du HED: 1.5 keV
                    !=>RESOLUTION SPATIALE: on somme les energies par pixels
                    !dans quel pixel de tab l'energie s'est elle deposee?
                    pixel=ED%Eabs(i)%pixel-min+1
                    if (tab(pixel(1),pixel(2)) .eq. 0) then        !on est en presence d'un nouveau pixel
                       n=n+1
                       tab(pixel(1),pixel(2))=n                    !on lui attribue un nouvel indice que l'on stock dans tab_pixels
                       ED%Epix(n)%pixel=ED%Eabs(i)%pixel           !on remplit la case correspondante avec sa position pixel
                       ED%Epix(n)%energie=E                        !on remplit la case correspondante avec son energie
                    else
                       indice=tab(pixel(1),pixel(2))               !on va chercher l'indice correspondant au pixel
                       ED%Epix(indice)%energie=ED%Epix(indice)%energie + E !on remplit la case correspondante avec son energie
                    endif
                 endif
              enddo
           endif
        endif

end subroutine resolution

!------------------------------------------------------------------------
!-         Algorithme de correction pour                                -
!-         l'imagerie du LED et du HED                                  -
!------------------------------------------------------------------------
subroutine correction(LED,HED)
implicit none

! IN/OUT variables
        type (Detecteur) :: LED !LED: nb d'interactions, liste des energies absorbees, ...
        type (Detecteur) :: HED !HED: nb d'interactions, liste des energies absorbees, ...

! Local variables
        integer, dimension(2) :: correction_LED, correction_HED !correction en nb de pixels
        integer :: i

!nb: on corrige uniquement de la position du LED ou du HED dans le referentiel lie aux miroirs,
!on neglige les effets de projection des photons dus a l'inclinaison du plan de detection (max=position/cos(2'))

        !calcul de la correction a appliquer (supposant une grande probabilite au centre du pixel)
        correction_LED(1)=(LED%err_Pos(1)-LED%dim_pixel/2) /LED%dim_pixel +1
        correction_LED(2)=(LED%err_Pos(2)-LED%dim_pixel/2) /LED%dim_pixel +1
        correction_HED(1)=(HED%err_Pos(1)-HED%dim_pixel/2) /HED%dim_pixel +1
        correction_HED(2)=(HED%err_Pos(2)-HED%dim_pixel/2) /HED%dim_pixel +1
        i=1
        do while (LED%Epix(i)%energie .ne. 0)
           LED%Epix(i)%pixel=LED%Epix(i)%pixel+correction_LED
           i=i+1
        enddo
        i=1
        do while (HED%Epix(i)%energie .ne. 0)
           HED%Epix(i)%pixel=HED%Epix(i)%pixel+correction_HED
           i=i+1
        enddo

end subroutine correction

!------------------------------------------------------------------------
!-         Filtre en Aluminium                                          -
!-                                                                      -
!------------------------------------------------------------------------
subroutine Al_filter(photon)
implicit none

! IN/OUT variables
        type (photons) :: photon        !position, direction et energie du photon

! Local variables
        real (kind=wp) :: Al_Photo, Al_Total
        real (kind=wp) :: G05CAF
        real (kind=wp) :: alea
        real (kind=wp) :: l


        !le photon traverse t'il le filtre d' Aluminium?
        call Al_CrossSection(photon%energie, Al_Photo,Al_Total)
        alea=G05CAF()
        l=-LOG(alea)/(Al_Total*2.6941_wp)       !distance avant prochaine interaction (masse volumique de Al=2.6941 g.cm-3)
        if (l .le. 0.00001_wp) then     !=>LE PHOTON NE TRAVERSE PAS LE FILTRE (epaisseur d'aluminium: 100 nm)
           photon%energie=0
        endif

end subroutine Al_filter

!------------------------------------------------------------------------
!-         Interactions rayonnement/matiere                             -
!-         dans le LED                                                  -
!------------------------------------------------------------------------
subroutine inside_LED(LED, photon)
implicit none

! IN/OUT variables
        type (Detecteur) :: LED         !nb d'interactions, liste des energies absorbees, ...
        type (photons) :: photon        !position, direction et energie du photon

! Local variables
        logical :: dedans
        real (kind=wp) :: Si_Photo, Si_Total
        real (kind=wp) :: G05CAF
        real (kind=wp) :: alea
        real (kind=wp) :: l
        integer, dimension(2) :: pixel
        real (kind=wp) :: E                             !energie absorbee


dedans=.TRUE.

        do while ((photon%energie .ne. 0) .and. dedans)

           !selon E0 aller chercher la section efficace du Si
           call Si_CrossSection(photon%energie, Si_Photo,Si_Total)
           !propagation du photon jusqu'a la prochaine interaction
           alea=G05CAF()
           l=-LOG(alea)/(Si_Total*2.329_wp)  !distance avant prochaine interaction (masse volumique du Si=2.329 g.cm-3)
           photon%position=photon%position+l*photon%direction

           !l'interaction du photon se situe t'elle a l'interieur du detecteur?
           if ((ABS(photon%position(1)) .le. LED%dim(1)/2) .and. (ABS(photon%position(2)) .le. LED%dim(2)/2) &
                .and. (photon%position(3) .le. LED%h_top) .and. (photon%position(3) .ge. LED%h_bot)) then
              !dans quel pixel du LED se situe alors le photon:
              pixel=(/ (photon%position(1)+LED%dim(1)/2)/LED%dim_pixel +1 , (photon%position(2)+LED%dim(2)/2)/LED%dim_pixel +1 /)

              !=====================================INTERACTION AVEC LE SILICON=====================================
              LED%Eabs(LED%nb+1)%Silicon=.TRUE.
              alea=G05CAF()
              if (alea .le. Si_Photo/Si_Total) then     !=>EFFET PHOTOELECTRIQUE
                 call photoelectrique('Si',photon,E)
                 LED%Eabs(LED%nb+1)%energie=E
                 LED%Eabs(LED%nb+1)%Photo=.TRUE.
                 LED%Eabs(LED%nb+1)%position=photon%position
                 LED%Eabs(LED%nb+1)%pixel=pixel

              else                                      !=>EFFET COMPTON
                 call Compton(photon,E)
                 LED%Eabs(LED%nb+1)%energie=E
                 LED%Eabs(LED%nb+1)%Compton=.TRUE.
                 LED%Eabs(LED%nb+1)%position=photon%position
                 LED%Eabs(LED%nb+1)%pixel=pixel
              endif
              LED%nb=LED%nb+1   !nombre d'interactions pour un photon dans le LED

           else
              !=>LE PHOTON N'EST PAS DANS UN PIXEL
              LED%Eabs(LED%nb+1)%photon_perdu=.TRUE.
              dedans=.FALSE.
           endif

           !le photon a t'il traverse le detecteur (en direction du HED)?
           if (photon%position(3) .lt. LED%h_bot) then
              !si oui 'repositionner' ce photon juste a la sortie du detecteur
              photon%position=photon%position+(LED%h_bot-photon%position(3))*photon%direction/photon%direction(3)
              photon%position(3)=LED%h_bot      !pour une valeur non approchee pour la condition .le. LED%h_bot de la boucle dans detecteurs
           endif

        enddo

end subroutine inside_LED

!------------------------------------------------------------------------
!-         Interactions rayonnement/matiere                             -
!-         dans le HED                                                  -
!------------------------------------------------------------------------
subroutine inside_HED(HED, photon)
implicit none

! IN/OUT variables
        type (Detecteur) :: HED         !nb d'interactions, liste des energies absorbees, ...
        type (photons) :: photon        !position, direction et energie du photon

! Local variables
        logical :: dedans
        real (kind=wp) :: Cd_Photo, Cd_Total
        real (kind=wp) :: Te_Photo, Te_Total
        real (kind=wp) :: Cd_CdTe, CdTe_Total
        real (kind=wp) :: G05CAF
        real (kind=wp) :: alea
        real (kind=wp) :: l
        integer, dimension(2) :: pixel
        real (kind=wp) :: E                             !energie absorbee


dedans=.TRUE.

        do while ((photon%energie .ne. 0) .and. dedans)

           !selon E0 aller chercher la section efficace de chaque element
           call Cd_CrossSection(photon%energie, Cd_Photo,Cd_Total)
           call Te_CrossSection(photon%energie, Te_Photo,Te_Total)
           !fraction en masse du Cd:0.468358, Te:0.531642
           CdTe_Total=0.468358_wp*Cd_Total+0.531642_wp*Te_Total
           !propagation du photon jusqu'a la prochaine interaction
           alea=G05CAF()
           l=-LOG(alea)/(CdTe_Total*5.85_wp)    !distance avant prochaine interaction (masse volumique du CdTe=5.85 g.cm-3)
           photon%position=photon%position+l*photon%direction

           !l'interaction du photon se situe t'elle a l'interieur du detecteur?
           if ((ABS(photon%position(1)) .le. HED%dim(1)/2) .and. (ABS(photon%position(2)) .le. HED%dim(2)/2) &
                .and. (photon%position(3) .le. HED%h_top) .and. (photon%position(3) .ge. HED%h_bot)) then
              !dans quel secteur du HED se situe alors le photon:
              pixel=(/ (photon%position(1)+HED%dim(1)/2)/HED%dim_pixel +1 , (photon%position(2)+HED%dim(2)/2)/HED%dim_pixel +1 /)
              if (HED%map(pixel(1),pixel(2)) .eq. 1) then       !=>LE PHOTON SE SITUE DANS UN PIXEL

                 !interaction avec le Cd ou le Te: rapport de la section efficace de Cd par rapport au CdTe
                 Cd_CdTe=(0.468358_wp*Cd_Total)/CdTe_Total

                 alea=G05CAF()
                 !=====================================INTERACTION AVEC LE CADMIUM=====================================
                 if (alea .le. Cd_CdTe) then
                    HED%Eabs(HED%nb+1)%Cadmium=.TRUE.
                    alea=G05CAF()
                    if (alea .le. Cd_Photo/Cd_Total) then     !=>EFFET PHOTOELECTRIQUE
                       call photoelectrique('Cd',photon,E)
                       HED%Eabs(HED%nb+1)%energie=E
                       HED%Eabs(HED%nb+1)%Photo=.TRUE.
                       HED%Eabs(HED%nb+1)%position=photon%position
                       HED%Eabs(HED%nb+1)%pixel=pixel

                    else                                      !=>EFFET COMPTON
                       call Compton(photon,E)
                       HED%Eabs(HED%nb+1)%energie=E
                       HED%Eabs(HED%nb+1)%Compton=.TRUE.
                       HED%Eabs(HED%nb+1)%position=photon%position
                       HED%Eabs(HED%nb+1)%pixel=pixel
                    endif

                    !====================================INTERACTION AVEC LE TELLURIUM====================================
                 else
                    HED%Eabs(HED%nb+1)%Tellurium=.TRUE.
                    alea=G05CAF()
                    if (alea .le. Te_Photo/Te_Total)then      !=>EFFET PHOTOELECTRIQUE
                       call photoelectrique('Te',photon,E)
                       HED%Eabs(HED%nb+1)%energie=E
                       HED%Eabs(HED%nb+1)%Photo=.TRUE.
                       HED%Eabs(HED%nb+1)%position=photon%position
                       HED%Eabs(HED%nb+1)%pixel=pixel

                    else                                      !=>EFFET COMPTON
                       call Compton(photon,E)
                       HED%Eabs(HED%nb+1)%energie=E
                       HED%Eabs(HED%nb+1)%Compton=.TRUE.
                       HED%Eabs(HED%nb+1)%position=photon%position
                       HED%Eabs(HED%nb+1)%pixel=pixel
                    endif

                 endif
                 HED%nb=HED%nb+1        !nombre d'interactions pour un photon dans le HED

              else
                 !=>LE PHOTON N'EST PAS DANS UN PIXEL
                 HED%Eabs(HED%nb+1)%photon_perdu=.TRUE.
                 dedans=.FALSE.
              endif

           else
              !=>LE PHOTON N'EST PAS DANS UN PIXEL
              HED%Eabs(HED%nb+1)%photon_perdu=.TRUE.
              dedans=.FALSE.
           endif

           !le photon a t'il traverse le detecteur (en direction du LED)?
           if (photon%position(3) .gt. HED%h_top) then
              !si oui 'repositionner' ce photon juste a la sortie du detecteur
              photon%position=photon%position+(HED%h_top-photon%position(3))*photon%direction/photon%direction(3)
              photon%position(3)=HED%h_top      !pour une valeur non approchee pour la condition .ge. HED%h_top de la boucle dans detecteurs
           endif

        enddo

end subroutine inside_HED

!------------------------------------------------------------------------
!-         Effet Photoelectrique                                        -
!-         (photon de fluorescence + electron Auger)                    -
!------------------------------------------------------------------------
subroutine photoelectrique(element,photon,Ea)
implicit none

! IN/OUT variables
        character (len=2) :: element    !element avec lequel le photon interagit
        type (photons) :: photon        !position, direction et energie du photon
        real (kind=wp) :: Ea            !energie absorbee


! Local variables
        real (kind=wp) :: G05CAF, G05DAF
        real (kind=wp) :: alea
        real (kind=wp) :: Pi
        real (kind=wp) :: alpha, delta


Pi=ACOS(-1._wp)

        !=======================================INTERACTION AVEC LE SILICON=======================================
        if (element .eq. 'Si') then

           alea=G05CAF()
           !=>COUCHE K si energie>Kedge et alea>rapport L/K
           if ((photon%energie .ge. 0.001839_wp) .and. (alea .gt. 0.096_wp)) then
              Ea=photon%energie-0.001839_wp             !Energie absorbee = energie du photon - energie de liaison 1s
              !une fois l'electron arrache de la couche K la desexcitation de l'atome entraine:
              alea=G05CAF()
              !emission d'un photon de fluorescence
              if (alea .lt. 0.91_wp) then
                 photon%energie=0.00174_wp              !energie de transition moyenne vers la couche K
                 Ea=Ea+(0.001839_wp-0.00174_wp)         !complement absorbe
                 !nouvelle direction du photon
                 alpha=G05DAF(-Pi,Pi)
                 delta=G05DAF(0,2*Pi)
                 call anglesolide(photon%direction,alpha,delta)
              !ejection d'un electron Auger
              else
                 Ea=Ea+0.001839_wp                      !toute l'energie est absorbee
                 photon%energie=0                       !il n'y a plus de photon
              endif

           !=>COUCHE L sinon
           else
              Ea=photon%energie                         !on fait l'approximation que toute l'energie est absorbee
              photon%energie=0                          !il n'y a plus de photon
           endif

        !=======================================INTERACTION AVEC LE CADMIUM=======================================
        else if (element .eq. 'Cd') then

           alea=G05CAF()
           !=>COUCHE K si energie>Kedge et alea>rapport L/K
           if ((photon%energie .ge. 0.02671_wp) .and. (alea .gt. 0.159_wp)) then
              Ea=photon%energie-0.02671_wp              !Energie absorbee = energie du photon - energie de liaison 1s
              !une fois l'electron arrache de la couche K la desexcitation de l'atome entraine:
              alea=G05CAF()
              !emission d'un photon de fluorescence
              if (alea .lt. 0.91_wp) then
                 photon%energie=0.02317_wp              !energie de transition moyenne vers la couche K
                 Ea=Ea+(0.02671_wp-0.02317_wp)          !complement absorbe
                 !nouvelle direction du photon
                 alpha=G05DAF(-Pi,Pi)
                 delta=G05DAF(0,2*Pi)
                 call anglesolide(photon%direction,alpha,delta)
              !ejection d'un electron Auger
              else
                 Ea=Ea+0.02671_wp                       !toute l'energie est absorbee
                 photon%energie=0                       !il n'y a plus de photon
              endif

           !=>COUCHE L sinon
           else
              Ea=photon%energie                         !on fait l'approximation que toute l'energie est absorbee
              photon%energie=0                          !il n'y a plus de photon
           endif

        !======================================INTERACTION AVEC LE TELLURIUM======================================
        else if (element .eq. 'Te') then

           alea=G05CAF()
           !=>COUCHE K si energie>Kedge et alea>rapport L/K
           if ((photon%energie .ge. 0.03181_wp) .and. (alea .gt. 0.165_wp)) then
              Ea=photon%energie-0.03181_wp              !Energie absorbee = energie du photon - energie de liaison 1s
              !une fois l'electron arrache de la couche K la desexcitation de l'atome entraine:
              alea=G05CAF()
              !emission d'un photon de fluorescence
              if (alea .lt. 0.91_wp) then
                 photon%energie=0.02747_wp              !energie de transition moyenne vers la couche K
                 Ea=Ea+(0.03181_wp-0.02747_wp)          !complement absorbe
                 !nouvelle direction du photon
                 alpha=G05DAF(-Pi,Pi)
                 delta=G05DAF(0,2*Pi)
                 call anglesolide(photon%direction,alpha,delta)
              !ejection d'un electron Auger
              else
                 Ea=Ea+0.03181_wp                       !toute l'energie est absorbee
                 photon%energie=0                       !il n'y a plus de photon
              endif

           !=>COUCHE L sinon
           else
              Ea=photon%energie                         !on fait l'approximation que toute l'energie est absorbee
              photon%energie=0                          !il n'y a plus de photon
           endif

        endif


end subroutine photoelectrique

!------------------------------------------------------------------------
!-         Effet Compton                                                -
!-         (diffusion du photon incident)                               -
!------------------------------------------------------------------------
subroutine Compton(photon,Ea)
implicit none

! IN/OUT variables
        type (photons) :: photon        !position, direction et energie du photon
        real (kind=wp) :: Ea            !energie absorbee

! Local variables
        real (kind=wp) :: a                     !rapport entre l'energie du photon et l'energie incidente
        real (kind=wp) :: a1, a2, a3, a4        !coefficient du polynome F(X)
        real (kind=wp) :: Xmin, Xmax, X, dX     !ensemble de definition
        real (kind=wp) :: Fmin, Fmax, F         !bornes de F(X) sur l'ensemble de definition

        real (kind=wp) :: G05CAF, G05DAF
        real (kind=wp) :: alea
        integer :: i
        real (kind=wp) :: Pi
        real (kind=wp) :: theta, phi

Pi=ACOS(-1._wp)

        !soit la fonction de repartition de la distribution Compton en energie:
        !F(X)=(mc^2)/E1( X^2(1/2) + X(2mc^2/E1+(mc^2/E1)^2) + lnX(1-2mc^2/E1-2(mc^2/E1)^2) - 1/X(mc^2/E1)^2)
        !avec X=E2/E1
        a=0.511_wp/photon%energie
        a1=0.5_wp
        a2=2*a+a**2
        a3=1-2*a-2*a**2
        a4=-a**2

        !bornes de la fonction sur l'ensemble de definition [Xmin,Xmax]
        Xmin=1/(1+2/a)
        Xmax=1
        !F(X)=a*( a1*X**2 + a2*X + a3*LOG(X) + a4/X )
        Fmin=a*( a1*Xmin**2 + a2*Xmin + a3*LOG(Xmin) + a4/Xmin )
        Fmax=a*( a1*Xmax**2 + a2*Xmax + a3*LOG(Xmax) + a4/Xmax )

        !tirage uniforme sur l'intervalle [Fmin,Fmax]
        alea=G05DAF(Fmin,Fmax)
        !dicotomie pour trouver l'antecedent X (F(X)=alea) de cette valeur:
        dX=(Xmax-Xmin)/2
        X=Xmax-dX
        do i=1,10
           F=a*( a1*X**2 + a2*X + a3*LOG(X) + a4/X )
           dX=dX/2
           if (F .gt. alea) then
              X=X-dX
           else
              X=X+dX
           endif
        enddo

        !energie deposee
        Ea=photon%energie*(1-X)
        !energie du photon diffuse
        photon%energie=photon%energie*X
        !nouvelle direction du photon diffuse
        theta=ACOS(1-a*(1/X-1))
        phi=G05DAF(0,2*Pi)
        call anglesolide(photon%direction,theta,phi)


end subroutine Compton

!------------------------------------------------------------------------
!-         Surface efficace                                             -
!-         de l' Aluminum Z=13 [0.001,1.0]MeV                           -
!------------------------------------------------------------------------
subroutine Al_CrossSection(E0, aPhoto,aTotal)
implicit none

! IN/OUT variables
       real (kind=wp) :: E0
       real (kind=wp) :: aPhoto, aTotal

! Local variables
       real (kind=wp) :: aCompt
       real (kind=wp), dimension(27) :: E, Compt, Photo
Data E     /0.001000_wp, 0.001500_wp, 0.001560_wp, 0.001560_wp, 0.0020000_wp, 0.0030000_wp, 0.0040_wp, &
            0.005000_wp, 0.006000_wp, 0.008000_wp, 0.010000_wp, 0.0150000_wp, 0.0200000_wp, 0.0300_wp, &
            0.040000_wp, 0.050000_wp, 0.060000_wp, 0.080000_wp, 0.1000000_wp, 0.1500000_wp, 0.2000_wp, &
            0.300000_wp, 0.400000_wp, 0.500000_wp, 0.600000_wp, 0.8000000_wp, 1.0000000_wp/
Data Compt /0.014300_wp, 0.024800_wp, 0.025900_wp, 0.025900_wp, 0.0337000_wp, 0.0473000_wp, 0.0581_wp, &
            0.067900_wp, 0.077000_wp, 0.092900_wp, 0.106000_wp, 0.1270000_wp, 0.1370000_wp, 0.1460_wp, &
            0.149000_wp, 0.150000_wp, 0.148000_wp, 0.144000_wp, 0.1390000_wp, 0.1270000_wp, 0.1170_wp, &
            0.102000_wp, 0.091600_wp, 0.083700_wp, 0.077500_wp, 0.0681000_wp, 0.0613000_wp/
Data Photo /1180.000_wp, 400.0000_wp, 360.0000_wp, 3960.000_wp, 2260.0000_wp, 787.00000_wp, 359.00_wp, &
            192.0000_wp, 114.0000_wp, 49.50000_wp, 25.60000_wp, 7.5100000_wp, 3.1000000_wp, 0.8720_wp, &
            0.350000_wp, 0.172000_wp, 0.095600_wp, 0.037800_wp, 0.0184000_wp, 0.0049900_wp, 0.0020_wp, &
            0.000574_wp, 0.000248_wp, 0.000134_wp, 0.000084_wp, 0.0000425_wp, 0.0000264_wp/
       integer :: i
       real (kind=wp) :: coeff


!Quel est le mu/rho de l' Aluminum pour un photon d'energie E0?
       if (E0 .eq. 0.00156_wp) then    !point double de la courbe: Kedge
          aCompt=0.0259_wp
          aPhoto=3960.0_wp
          aTotal=aCompt+aPhoto
       else                             !interpolation sur la droite log
          i=2
          do while (E(i) .lt. E0)       !on cherche l'interval [i-1,i]
             i=i+1
          enddo
          coeff=LOG10(Photo(i)/Photo(i-1))/LOG10(E(i)/E(i-1))
          aPhoto=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Photo(i-1)))
          coeff=LOG10(Compt(i)/Compt(i-1))/LOG10(E(i)/E(i-1))
          aCompt=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Compt(i-1)))
          aTotal=aCompt+aPhoto
       endif

end subroutine Al_CrossSection

!------------------------------------------------------------------------
!-         Surface efficace                                             -
!-         du Silicon Z=14 [0.001,1.0]MeV                               -
!------------------------------------------------------------------------
subroutine Si_CrossSection(E0, aPhoto,aTotal)
implicit none

! IN/OUT variables
       real (kind=wp) :: E0
       real (kind=wp) :: aPhoto, aTotal

! Local variables
       real (kind=wp) :: aCompt
       real (kind=wp), dimension(27) :: E, Compt, Photo
Data E     /0.001000_wp, 0.001500_wp, 0.001839_wp, 0.001839_wp, 0.0020000_wp, 0.0030000_wp, 0.00400_wp, &
            0.005000_wp, 0.006000_wp, 0.008000_wp, 0.010000_wp, 0.0150000_wp, 0.0200000_wp, 0.03000_wp, &
            0.040000_wp, 0.050000_wp, 0.060000_wp, 0.080000_wp, 0.1000000_wp, 0.1500000_wp, 0.20000_wp, &
            0.300000_wp, 0.400000_wp, 0.500000_wp, 0.600000_wp, 0.8000000_wp, 1.0000000_wp/
Data Compt /0.013200_wp, 0.023900_wp, 0.030800_wp, 0.030800_wp, 0.0339000_wp, 0.0496000_wp, 0.06130_wp, &
            0.071100_wp, 0.079800_wp, 0.095100_wp, 0.108000_wp, 0.1290000_wp, 0.1400000_wp, 0.15000_wp, &
            0.153000_wp, 0.154000_wp, 0.153000_wp, 0.148000_wp, 0.1430000_wp, 0.1310000_wp, 0.12100_wp, &
            0.106000_wp, 0.094800_wp, 0.086600_wp, 0.080200_wp, 0.0705000_wp, 0.0634000_wp/
Data Photo /1570.000_wp, 533.0000_wp, 307.0000_wp, 3190.000_wp, 2770.0000_wp, 977.00000_wp, 451.000_wp, &
            244.0000_wp, 146.0000_wp, 63.80000_wp, 33.10000_wp, 9.8500000_wp, 4.0900000_wp, 1.16000_wp, &
            0.469000_wp, 0.231000_wp, 0.129000_wp, 0.051200_wp, 0.0250000_wp, 0.0068100_wp, 0.00274_wp, &
            0.000788_wp, 0.000341_wp, 0.000185_wp, 0.000116_wp, 0.0000585_wp, 0.0000364_wp/
       integer :: i
       real (kind=wp) :: coeff


!Quel est le mu/rho du Silicon pour un photon d'energie E0?
       if (E0 .eq. 0.001839_wp) then    !point double de la courbe: Kedge
          aCompt=0.0308_wp
          aPhoto=3190.0_wp
          aTotal=aCompt+aPhoto
       else                             !interpolation sur la droite log
          i=2
          do while (E(i) .lt. E0)       !on cherche l'interval [i-1,i]
             i=i+1
          enddo
          coeff=LOG10(Photo(i)/Photo(i-1))/LOG10(E(i)/E(i-1))
          aPhoto=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Photo(i-1)))
          coeff=LOG10(Compt(i)/Compt(i-1))/LOG10(E(i)/E(i-1))
          aCompt=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Compt(i-1)))
          aTotal=aCompt+aPhoto
       endif

end subroutine Si_CrossSection

!------------------------------------------------------------------------
!-         Surface efficace                                             -
!-         du Cadmium Z=48 [0.005,1.0]MeV                               -
!------------------------------------------------------------------------
subroutine Cd_CrossSection(E0, aPhoto,aTotal)
implicit none

! IN/OUT variables
       real (kind=wp) :: E0
       real (kind=wp) :: aPhoto, aTotal

! Local variables
       real (kind=wp) :: aCompt
       real (kind=wp), dimension(22) :: E, Compt, Photo
Data E     /0.00500_wp, 0.0060_wp, 0.0080_wp, 0.0100_wp, 0.0150_wp, 0.0200_wp, 0.02671_wp, &
            0.02671_wp, 0.0300_wp, 0.0400_wp, 0.0500_wp, 0.0600_wp, 0.0800_wp, 0.10000_wp, &
            0.15000_wp, 0.2000_wp, 0.3000_wp, 0.4000_wp, 0.5000_wp, 0.6000_wp, 0.80000_wp, 1.00000_wp/
Data Compt /0.03550_wp, 0.0417_wp, 0.0528_wp, 0.0621_wp, 0.0790_wp, 0.0899_wp, 0.09970_wp, &
            0.09970_wp, 0.1030_wp, 0.1100_wp, 0.1140_wp, 0.1150_wp, 0.1150_wp, 0.11300_wp, &
            0.10700_wp, 0.0998_wp, 0.0885_wp, 0.0800_wp, 0.0734_wp, 0.0682_wp, 0.06010_wp, 0.05410_wp/
Data Photo /764.000_wp, 475.00_wp, 222.00_wp, 122.00_wp, 40.000_wp, 17.900_wp, 7.93000_wp, &
            49.8000_wp, 36.900_wp, 17.200_wp, 9.3600_wp, 5.6400_wp, 2.5000_wp, 1.32000_wp, &
            0.40700_wp, 0.1770_wp, 0.0561_wp, 0.0257_wp, 0.0144_wp, 0.0092_wp, 0.00475_wp, 0.00297_wp/
       integer :: i
       real (kind=wp) :: coeff


!Quel est le mu/rho du Cadmium pour un photon d'energie E0?
       if (E0 .eq. 0.02671_wp) then     !point double de la courbe: Kedge
          aCompt=0.0997_wp
          aPhoto=49.800_wp
          aTotal=aCompt+aPhoto
       else                             !interpolation sur la droite log
          i=2
          do while (E(i) .lt. E0)       !on cherche l'interval [i-1,i]
             i=i+1
          enddo
          coeff=LOG10(Photo(i)/Photo(i-1))/LOG10(E(i)/E(i-1))
          aPhoto=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Photo(i-1)))
          coeff=LOG10(Compt(i)/Compt(i-1))/LOG10(E(i)/E(i-1))
          aCompt=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Compt(i-1)))
          aTotal=aCompt+aPhoto
       endif

end subroutine Cd_CrossSection

!------------------------------------------------------------------------
!-         Surface efficace                                             -
!-         du Tellurium Z=52 [0.005,1.0]MeV                             -
!------------------------------------------------------------------------
subroutine Te_CrossSection(E0, aPhoto,aTotal)
implicit none

! IN/OUT variables
       real (kind=wp) :: E0
       real (kind=wp) :: aPhoto, aTotal

! Local variables
       real (kind=wp) :: aCompt
       real (kind=wp), dimension(22) :: E, Compt, Photo
Data E     /0.00500_wp, 0.00600_wp, 0.0080_wp, 0.0100_wp, 0.0150_wp, 0.0200_wp, 0.03000_wp, &
            0.03181_wp, 0.03181_wp, 0.0400_wp, 0.0500_wp, 0.0600_wp, 0.0800_wp, 0.10000_wp, &
            0.15000_wp, 0.20000_wp, 0.3000_wp, 0.4000_wp, 0.5000_wp, 0.6000_wp, 0.80000_wp, 1.00000_wp/
Data Compt /0.03460_wp, 0.04010_wp, 0.0498_wp, 0.0582_wp, 0.0742_wp, 0.0846_wp, 0.09700_wp, &
            0.09840_wp, 0.09840_wp, 0.1030_wp, 0.1070_wp, 0.1090_wp, 0.1090_wp, 0.10800_wp, &
            0.10100_wp, 0.09490_wp, 0.0842_wp, 0.0763_wp, 0.0700_wp, 0.0650_wp, 0.05730_wp, 0.05160_wp/
Data Photo /897.000_wp, 568.000_wp, 267.00_wp, 147.00_wp, 48.900_wp, 22.100_wp, 7.07000_wp, &
            5.99000_wp, 36.4000_wp, 20.100_wp, 11.000_wp, 6.6900_wp, 3.0000_wp, 1.59000_wp, &
            0.49900_wp, 0.21900_wp, 0.0700_wp, 0.0322_wp, 0.0182_wp, 0.0116_wp, 0.00602_wp, 0.00376_wp/
       integer :: i
       real (kind=wp) :: coeff


!Quel est le mu/rho du Cadmium pour un photon d'energie E0?
       if (E0 .eq. 0.03181_wp) then     !point double de la courbe: Kedge
          aCompt=0.0984_wp
          aPhoto=36.400_wp
          aTotal=aCompt+aPhoto
       else                             !interpolation sur la droite log
          i=2
          do while (E(i) .lt. E0)       !on cherche l'interval [i-1,i]
             i=i+1
          enddo
          coeff=LOG10(Photo(i)/Photo(i-1))/LOG10(E(i)/E(i-1))
          aPhoto=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Photo(i-1)))
          coeff=LOG10(Compt(i)/Compt(i-1))/LOG10(E(i)/E(i-1))
          aCompt=10._wp**(coeff*LOG10(E0/E(i-1))+LOG10(Compt(i-1)))
          aTotal=aCompt+aPhoto
       endif

end subroutine Te_CrossSection

!------------------------------------------------------------------------
!-      Distribution du vecteur photon dans l'angle solide              -
!------------------------------------------------------------------------
subroutine anglesolide(vecteur,alpha,delta)

implicit none

! IN/OUT variables
        real (kind=wp), dimension(3) :: vecteur
        real (kind=wp) :: alpha, delta

! local variables
        real (kind=wp), dimension(3) :: axe, temp

!produit vectoriel de vecteur par n'importe quoi pour trouver axe de rotation orthogonal:
        !si vecteur est colineaire a (0,0,1) alors
        if ((vecteur(1) .eq. 0) .and. (vecteur(2) .eq. 0)) then
           axe(1)=0             !vecteur^(1,0,0)
           axe(2)=vecteur(3)
           axe(3)=-vecteur(2)
        !si vecteur n'est pas colineaire a (0,0,1) alors
        else
           axe(1)=vecteur(2)    !vecteur^(0,0,1)
           axe(2)=-vecteur(1)
           axe(3)=0
        endif

        temp=vecteur
!rotation de alpha autour de l'axe orthogonal au vecteur
        call rotation_rad(axe,vecteur,alpha)
!rotation de delta autour du vecteur d'origine
        call rotation_rad(temp,vecteur,delta)

end subroutine anglesolide


end module FPA
