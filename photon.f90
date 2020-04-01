module PHOTON

use OUTILS
use LECTURE

implicit none


contains


!------------------------------------------------------------------------
!-	Genere un photon sur la face d'entree des miroirs dont la       -
!-   direction provient d'une source dans le repere RefM des miroirs    -
!------------------------------------------------------------------------
subroutine tirage_photon(Para_SM,Obs,RefM,n, photon)

implicit none

!IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM
        type (Observation) :: Obs
        type (Repere) :: RefM
        integer :: n

        type (photons) :: photon

!local variables
        real (kind=wp) :: x, y
        real (kind=wp) :: alpha, delta
        integer :: i
        real (kind=wp) :: G05DAF
        real (kind=wp) :: Emin, Emax, E
        real (kind=wp) :: indice
        real (kind=wp) :: Fmin, Fmax

!genere une position pour le photon sur la face d'entree des miroirs
        !(dans un carre de max Para_SM%rayon*2 de cote pour couvrir tous les miroirs)
        x=G05DAF(-maxval(Para_SM%rayon),maxval(Para_SM%rayon))
        !x=0.0 !pour faire un profil 2D gaussien
        y=G05DAF(-maxval(Para_SM%rayon),maxval(Para_SM%rayon))

        photon%position(1)=x
        photon%position(2)=y
        photon%position(3)=0.3_wp  !(+0.3m du milieu du module miroir)

        !ce photon va t-il rencontrer un miroir, si oui lequel?
        Para_SM%numero_miroir=0
        Para_SM%entreemiroirs=.FALSE.
        do i=1,Para_SM%nb_miroirs
           if (SQRT(x**2._wp+y**2._wp) .lt. Para_SM%rayon(i)) then
              if (Para_SM%numero_miroir==0) then
                 Para_SM%numero_miroir=i
                 Para_SM%entreemiroirs=.TRUE.
              endif
           endif
        enddo

if (Para_SM%entreemiroirs) then
!genere la direction du photon
        !passage des coordonnees equatoriales de la source dans le repere cartesien du pointage du telescope
        call conv_equat_to_instr(Para_SM%pointage%attitude,Obs%sources(n)%coordonnees, photon%direction)

        !modification de la direction du photon dans le repere local instrument etant donne les derives de celui-ci
        call rotationZ(photon%direction,-RefM%attitude(3))
        call rotationY(photon%direction,-RefM%attitude(2))
        call rotationX(photon%direction,-RefM%attitude(1))

!genere une energie pour le photon
	Emin=Obs%minmax_energie(1)
	Emax=Obs%minmax_energie(2)

        if (Obs%energie_unique) then
            !mono-energie des photons
            photon%energie=Obs%energie_value

        elseif (Obs%energie_uniforme) then
            !distribution uniforme
                E=G05DAF(Emin,Emax)
            photon%energie=E

        elseif (Obs%energie_powerlaw) then
            !distribution selon une loi de probabilite f(x)=x^(-indice)
	    	indice=Obs%sources(n)%indice
            	!primitive F(x)=1/(-indice+1)*x^(-indice+1)
            	!calculer les bornes de F(x) sur le domaine de definition (bande en energie):
	    	Fmax=1._wp/(-indice+1._wp)*Emax**(-indice+1._wp)
	    	Fmin=1._wp/(-indice+1._wp)*Emin**(-indice+1._wp)
            	!effectuer un tirage aleatoire uniforme sur [Fmin,Fmax]:
                E=G05DAF(Fmin,Fmax)
		!calculer la fonction inverse F-1(x)=((-indice+1)*x)^(1/(-indice+1))
		!calculer les images de x par la fonction F-1(x):
		E=((-indice+1._wp)*E)**(1._wp/(-indice+1._wp))
            photon%energie=E

        endif
endif

end subroutine tirage_photon


!------------------------------------------------------------------------
!-	Introduction d'erreurs simulant le mauvais assemblage        	-
!-      des coquilles en alignement (tilt)                              -
!------------------------------------------------------------------------
subroutine tilt(Para_SM, photon)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM

        type (photons) :: photon

!introduction d'erreurs sur la direction du photon incident pour simuler:
    !les defauts d'alignement des shells (tilt)
    if (Para_SM%assembly) then
    if (Para_SM%entreemiroirs) then
        call rotationX(photon%direction,-Para_SM%tilt(1,Para_SM%numero_miroir))
        call rotationY(photon%direction,-Para_SM%tilt(2,Para_SM%numero_miroir))
    endif
    endif
!normalement il faudrait aussi effectuer un changement de position du photon du au tilt des shells,
!mais cela impliquerait qu'il interagirait avec une autre shell (inconnue) alors que statistiquement
!une shell a autant de probabilite de recevoir un photon, qu'elle soit tiltee ou non. Comme on exprime
!la position du photon dans son repere, on ne modifie pas sa position

end subroutine tilt


!------------------------------------------------------------------------
!-	Calcul de la deviation du photon par les miroirs        	-
!-      dans le repere M du satellite miroir                            -
!------------------------------------------------------------------------
subroutine miroir(Para_SM,Obs,photon, photonR2)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM
        type (Observation) :: Obs
        type (photons) :: photon

        type (photons) :: photonR2

! local variables
        integer :: m
        type (photons) :: photonR1
	real (kind=wp), dimension(3) :: abc_photon, abc_parabole, abc_hyperbole
	real (kind=wp)  :: a, b, c, delta
	real (kind=wp)  :: z1, z2, r1, r2
        real (kind=wp), dimension(3) :: normale, tangente, projnormale
        real (kind=wp) :: G05DAF, G05DDF
        real (kind=wp), dimension(3) :: F
        real (kind=wp) :: A
        real (kind=wp) :: incidence
        real (kind=wp) :: R
        logical :: energie
        real (kind=wp) :: proba

! Initialisation des variables locales
    m=Para_SM%numero_miroir
    Para_SM%reflexionparabole=.FALSE.
    Para_SM%reflexionhyperbole=.FALSE.
    energie=Obs%energie_unique .or. Obs%energie_uniforme .or. Obs%energie_powerlaw

if (Para_SM%entreemiroirs) then
!reflexion sur le miroir parabolique
    !point d'intersection du photon avec le miroir
        abc_photon(1)=(photon%direction(2)/abs(photon%direction(3)))**2._wp &
                      + (photon%direction(1)/abs(photon%direction(3)))**2._wp
        abc_photon(2)=2._wp*photon%position(2)*(photon%direction(2)/abs(photon%direction(3))) &
                      + 2._wp*photon%position(1)*(photon%direction(1)/abs(photon%direction(3)))
        abc_photon(3)=photon%position(2)**2._wp + photon%position(1)**2._wp
        !print *, abc_photon
        abc_parabole(1)=0._wp
        abc_parabole(2)=-2._wp*Para_SM%paraboles(m)%parametre
        abc_parabole(3)=2._wp*Para_SM%paraboles(m)%parametre*(Para_SM%paraboles(m)%origine+photon%position(3)) &
                        -Para_SM%paraboles(m)%parametre**2._wp
        !print *, abc_parabole

        a=abc_photon(1)-abc_parabole(1)
        b=abc_photon(2)-abc_parabole(2)
        c=abc_photon(3)-abc_parabole(3)
        delta=b**2._wp - 4._wp*a*c
        if (a == 0._wp) then
           z1=-c/b
           z2=0._wp
        elseif (delta >= 0._wp) then
           z1=(-b-SQRT(delta))/(2._wp*a)
           z2=(-b+SQRT(delta))/(2._wp*a)
        else
           z1=0._wp
           z2=0._wp
        endif

        !on exprime les solutions dans le repere miroirs
        r1=photon%position(3)-z1
        r2=photon%position(3)-z2
        !print*
        !print*, "solutions pour l'intersection 1, z="
        !print*,photon%position(3),z1,z2,a
        !print*, r1,r2,delta

        !la solution qui nous interesse est sur le miroir
        if ((r1 < 0.3_wp) .and. (r1 > 0._wp)) then
           Para_SM%reflexionparabole=.TRUE.
           photonR1%position(1)=photon%position(1)+z1*(photon%direction(1)/abs(photon%direction(3)))
           photonR1%position(2)=photon%position(2)+z1*(photon%direction(2)/abs(photon%direction(3)))
           photonR1%position(3)=r1
        elseif 	((r2 < 0.3_wp) .and. (r2 > 0._wp)) then
           Para_SM%reflexionparabole=.TRUE.
           photonR1%position(1)=photon%position(1)+z2*(photon%direction(1)/abs(photon%direction(3)))
           photonR1%position(2)=photon%position(2)+z2*(photon%direction(2)/abs(photon%direction(3)))
           photonR1%position(3)=r2
        else	
           !sortir du calcul pour ce photon
           Para_SM%reflexionparabole=.FALSE.
        endif

    !reflexion du photon en ce point
        !on norme le vecteur direction reflechi du photon
        !photon%direction=photon%direction/NORME(photon%direction)

        if (Para_SM%reflexionparabole) then
           normale(1)=-photonR1%position(1)
           normale(2)=-photonR1%position(2)
           normale(3)=Para_SM%paraboles(m)%parametre
           tangente(1)=photonR1%position(2)
           tangente(2)=-photonR1%position(1)
           tangente(3)=0._wp

           if (Para_SM%figure) then
           !introduction d'erreur sur la normale pour simuler les defauts de surface
              !erreurs sur X et Y
              call rotationZ(normale,G05DDF(0._wp,Para_SM%forme))!*310.d0 pour un profil artificiel Gaussien
              !erreurs sur Z
              call rotation_deg(tangente,normale,G05DDF(0._wp,Para_SM%forme))
           endif

           if ((Para_SM%scattering) .and. (energie)) then
           !introduction d'erreur sur la normale pour simuler le XRS
              !erreurs sur Z (le XRS se situe essentiellement dans le plan d'incidence)
              call rotation_deg(tangente,normale,G05DDF(0._wp,Para_SM%spread*photon%energie))
           endif

           !calcul la nouvelle direction du photon suite a la 1ere reflexion
           projnormale=-normale*dot_product(photon%direction,normale) &
                       /(normale(1)**2._wp+normale(2)**2._wp+normale(3)**2._wp)
           photonR1%direction=photon%direction + 2._wp*projnormale

           !print*
           !print*, "position 1 du photon reflechi"
           !print*, photonR1%position
           !print*
           !print*, "vecteur normal au miroir au point d'intersection 1"
           !print*, normale
           !print*
           !print*, "direction 1 du photon reflechi"
           !print*, photonR1%direction

           !probabilite pour le photon d'etre reflechi en fonction de son energie et de son angle d'incidence
           if ( (energie) .and. (.not. Para_SM%reflexiontotale) ) then
               incidence=ACOSD( DOT_PRODUCT(photon%direction,normale) / (NORME(photon%direction)*NORME(normale)) ) - 90._wp
!incidence=0.08
               if ((incidence .ge. 0) .and. (incidence .le. 0.5)) then  !au dessus de 0.5 on sort forcement du plan de detection
write(32) incidence
                  !le tableau Obs%coeff_reflec(1000,1000) va de 0.1-100 keV et de 0-30 arcmin
                  R=Obs%coeff_reflec(INT(photon%energie/100._wp*1000._wp),INT(incidence*60._wp/30._wp*999._wp)+1)
                  proba=G05DAF(0._wp,1._wp)
                  if (proba .le. R) then
                     Para_SM%reflexionparabole=.TRUE.
                  else
                     !sortir du calcul pour ce photon
                     Para_SM%reflexionparabole=.FALSE.
                  endif

                  !print*
                  !print*,'miroir numero:',m
                  !print*,'energie:',photon%energie,'indice tab:',INT(photon%energie/100._wp*1000._wp)
                  !print*,'incidence:',incidence*60._wp,'indice tab:',INT(incidence*60._wp/30._wp*999._wp)+1
                  !print*,'coefficient de reflexion:',R
                  !print*,'photon reflechi ?',Para_SM%reflexionparabole
               endif
           endif
           photonR1%energie=photon%energie
        endif


!reflexion sur le miroir hyperbolique	
    !point d'intersection du photon avec le miroir
        if (Para_SM%reflexionparabole) then
           abc_photon(1)=(photonR1%direction(2)/abs(photonR1%direction(3)))**2._wp &
                         + (photonR1%direction(1)/abs(photonR1%direction(3)))**2._wp
           abc_photon(2)=2._wp*photonR1%position(2)*(photonR1%direction(2)/abs(photonR1%direction(3))) &
                         + 2._wp*photonR1%position(1)*(photonR1%direction(1)/abs(photonR1%direction(3)))
           abc_photon(3)=photonR1%position(2)**2._wp + photonR1%position(1)**2._wp
           !print *, abc_photon
           abc_hyperbole(1)=Para_SM%hyperboles(m)%e2-1._wp
           abc_hyperbole(2)=-2._wp*(Para_SM%hyperboles(m)%origine+photonR1%position(3))*(Para_SM%hyperboles(m)%e2-1._wp) &
                            -2._wp*Para_SM%hyperboles(m)%parametre
           abc_hyperbole(3)=(Para_SM%hyperboles(m)%origine+photonR1%position(3))**2._wp*(Para_SM%hyperboles(m)%e2-1._wp) &
                            +2._wp*Para_SM%hyperboles(m)%parametre*(Para_SM%hyperboles(m)%origine+photonR1%position(3)) &
                            -Para_SM%hyperboles(m)%parametre**2._wp
           !print *, abc_parabole

           a=abc_photon(1)-abc_hyperbole(1)
           b=abc_photon(2)-abc_hyperbole(2)
           c=abc_photon(3)-abc_hyperbole(3)
           delta=b**2._wp - 4._wp*a*c
           if (a == 0._wp) then
              z1=-c/b
              z2=0._wp
           elseif (delta >= 0._wp) then
              z1=(-b-SQRT(delta))/(2._wp*a)
              z2=(-b+SQRT(delta))/(2._wp*a)
           else
              z1=0._wp
              z2=0._wp
           endif

           !on exprime les solutions dans le repere miroirs
           r1=photonR1%position(3)-z1
           r2=photonR1%position(3)-z2
           !print*
           !print*, "solutions pour l'intersection 2, z="
           !print*, r1,r2

           !la solution qui nous interesse est sur le miroir
           if ((r1 < 0._wp) .and. (r1 > -0.3_wp)) then
              Para_SM%reflexionhyperbole=.TRUE.
              photonR2%position(1)=photonR1%position(1)+z1*(photonR1%direction(1)/abs(photonR1%direction(3)))
              photonR2%position(2)=photonR1%position(2)+z1*(photonR1%direction(2)/abs(photonR1%direction(3)))
              photonR2%position(3)=r1
           elseif ((r2 < 0._wp) .and. (r2 > -0.3_wp)) then
              Para_SM%reflexionhyperbole=.TRUE.
              photonR2%position(1)=photonR1%position(1)+z2*(photonR1%direction(1)/abs(photonR1%direction(3)))
              photonR2%position(2)=photonR1%position(2)+z2*(photonR1%direction(2)/abs(photonR1%direction(3)))
              photonR2%position(3)=r2
           else	
              !sortir du calcul pour ce photon
              Para_SM%reflexionhyperbole=.FALSE.
           endif

    !reflexion du photon en ce point
           !on norme le vecteur direction reflechi du photon
           !photonR1%direction=photonR1%direction/NORME(photonR1%direction)

           if (Para_SM%reflexionhyperbole) then
              normale(1)=-photonR2%position(1)
              normale(2)=-photonR2%position(2)
              normale(3)=(Para_SM%hyperboles(m)%origine+photonR2%position(3))*(Para_SM%hyperboles(m)%e2-1._wp) &
                         +Para_SM%hyperboles(m)%parametre
              tangente(1)=photonR2%position(2)
              tangente(2)=-photonR2%position(1)
              tangente(3)=0._wp

              if (Para_SM%figure) then
              !introduction d'erreur sur la normale pour simuler les defauts de surface
                 !erreurs sur X et Y
                 call rotationZ(normale,G05DDF(0._wp,Para_SM%forme))!*310.d0 pour un profil artificiel Gaussien
                 !erreurs sur Z
                 call rotation_deg(tangente,normale,G05DDF(0._wp,Para_SM%forme))
              endif

              if ((Para_SM%scattering) .and. (energie)) then
              !introduction d'erreur sur la normale pour simuler le XRS
                 !erreurs sur Z (le XRS se situe essentiellement dans le plan d'incidence)
                 call rotation_deg(tangente,normale,G05DDF(0._wp,Para_SM%spread*photon%energie))
              endif

              !calcul la nouvelle direction du photon suite a la 2eme reflexion
              projnormale=-normale*dot_product(photonR1%direction,normale) &
                          /(normale(1)**2._wp+normale(2)**2._wp+normale(3)**2._wp)
              photonR2%direction=photonR1%direction + 2._wp*projnormale

              !print*
              !print*, "position 2 du photon reflechi"
              !print*, photonR2%position
              !print*
              !print*, "vecteur normal au miroir au point d'intersection 2"
              !print*, normale
              !print*
              !print*, "direction 2 du photon reflechi"
              !print*, photonR2%direction

              !probabilite pour le photon d'etre reflechi en fonction de son energie et de son angle d'incidence
              if ( (energie) .and. (.not. Para_SM%reflexiontotale) ) then
                  incidence=ACOSD( DOT_PRODUCT(photonR1%direction,normale) / (NORME(photonR1%direction)*NORME(normale)) ) - 90._wp
!incidence=0.08
               if ((incidence .ge. 0) .and. (incidence .le. 0.5)) then  !au dessus de 0.5 on sort forcement du plan de detection
write(32) incidence
                     !le tableau Obs%coeff_reflec(1000,1000) va de 0.1-100 keV et de 0-30 arcmin
                     R=Obs%coeff_reflec(INT(photonR1%energie/100._wp*1000._wp),INT(incidence*60._wp/30._wp*999._wp)+1)
                     proba=G05DAF(0._wp,1._wp)
                     if (proba .le. R) then
                        Para_SM%reflexionhyperbole=.TRUE.
                     else
                        !sortir du calcul pour ce photon
                        Para_SM%reflexionhyperbole=.FALSE.
                     endif

                     !print*
                     !print*,'miroir numero:',m
                     !print*,'energie:',photonR1%energie,'indice tab:',INT(photonR1%energie/100._wp*1000._wp)
                     !print*,'incidence:',incidence*60._wp,'indice tab:',INT(incidence*60._wp/30._wp*999._wp)+1
                     !print*,'coefficient de reflexion:',R
                     !print*,'photon reflechi ?',Para_SM%reflexionhyperbole
                  endif
              endif
              photonR2%energie=photon%energie
           endif
        endif

    !print*,Para_SM%reflexionparabole,Para_SM%reflexionhyperbole

endif


end subroutine miroir


!------------------------------------------------------------------------
!-	Introduction d'erreurs simulant le mauvais assemblage        	-
!-      des coquilles selon les axes X, Y et Z                          -
!------------------------------------------------------------------------
subroutine XYZshellserrors(Para_SM, photon)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM

        type (photons) :: photon


!introduction d'erreurs sur la position du photon reflechi pour simuler:
    if (Para_SM%assembly) then
    if (Para_SM%reflexionhyperbole) then
        !les defauts selon X et Y des shells (decenter)
        photon%position(1)=photon%position(1)+Para_SM%delta(1,Para_SM%numero_miroir)
        photon%position(2)=photon%position(2)+Para_SM%delta(2,Para_SM%numero_miroir)
        !les defauts selon Z des shells (defocus)
        photon%position(3)=photon%position(3)+Para_SM%delta(3,Para_SM%numero_miroir)
    endif
    endif

end subroutine XYZshellserrors


!------------------------------------------------------------------------
!-      Exprime les coordonnees du photon dans le repere RefD           -
!-      du satellite detecteur                                          -
!------------------------------------------------------------------------
subroutine reperedetecteurs(Para_SM,photon,RefM,RefD, photon_D)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM
        type (photons) :: photon
        type (Repere) :: RefM
        type (Repere) :: RefD

        type (photons) :: photon_D

! local variables
        real (kind=wp), dimension(3) :: deriveM
        real (kind=wp), dimension(3) :: RefD_position_temp


!changement de repere pour le photon
    if (Para_SM%reflexionhyperbole) then
        !derive en attitude des miroirs par rapport aux detecteurs
        deriveM=RefM%attitude - RefD%attitude

            !derive d'attitude supplementaire du au mauvais alignement des shells (tilt)
            if (Para_SM%assembly) then
                deriveM(1)=deriveM(1)+Para_SM%tilt(1,Para_SM%numero_miroir)
                deriveM(2)=deriveM(2)+Para_SM%tilt(2,Para_SM%numero_miroir)
            endif

        RefD_position_temp=RefD%position    !pour ne pas ecraser RefD%position (plusieurs photons pour le meme RefD%position)

            !la derive supplementaire (tilt) induit un changement de position de RefD exprime dans RefM
            if (Para_SM%assembly) then
                call rotationY(RefD_position_temp,-Para_SM%tilt(2,Para_SM%numero_miroir))
                call rotationX(RefD_position_temp,-Para_SM%tilt(1,Para_SM%numero_miroir))
            endif

        !position du photon dans le repere RefD
	photon_D%position = photon%position-RefD_position_temp
        call rotationZ(photon_D%position,deriveM(3))
        call rotationY(photon_D%position,deriveM(2))
        call rotationX(photon_D%position,deriveM(1))

        !direction du photon dans le repere RefD
        photon_D%direction=photon%direction
        call rotationZ(photon_D%direction,deriveM(3))
        call rotationY(photon_D%direction,deriveM(2))
        call rotationX(photon_D%direction,deriveM(1))

        !energie du photon
        photon_D%energie=photon%energie

    endif

end subroutine reperedetecteurs


!------------------------------------------------------------------------
!-                      Programme principal                             -
!------------------------------------------------------------------------
subroutine parcoursphoton(Obs,Para_SM,Para_SD,RefM,RefD,n, photon_D)

implicit none

! IN/OUT variables
        type (Observation) :: Obs
        type (ParametresSatMiroir) :: Para_SM
        type (ParametresSatDetecteur) :: Para_SD
        type (Repere) :: RefM                       !referentiel des miroirs
        type (Repere) :: RefD                       !referentiel des detecteurs
        integer :: n                                !numero de la source en cours de simulation

        type (photons) :: photon_D

! local variables
        type (photons) :: photon_M
        type (photons) :: photonreflechi_M


! Generation d'un photon
        call tirage_photon(Para_SM,Obs,RefM,n, photon_M)

! Simulation d'erreur d'alignement des coquilles
        call tilt(Para_SM, photon_M)

! Transport du photon a travers les miroirs
        call miroir(Para_SM,Obs,photon_M, photonreflechi_M)

! Simulation d'erreur de positionnement des coquilles
        call XYZshellserrors(Para_SM, photonreflechi_M)

! Position du photon dans le repere des detecteurs
        call reperedetecteurs(Para_SM,photonreflechi_M,RefM,RefD, photon_D)

! Optionnel: pour visualiser l'origine des photons n'ayant pas ete reflechis
        if (Para_SM%reflexionhyperbole) then
            if (Para_SM%assembly) then
                photon_M%position=photon_M%position+Para_SM%delta(:,Para_SM%numero_miroir)
                call rotationX(photon_M%position,-Para_SM%tilt(1,Para_SM%numero_miroir))
                call rotationY(photon_M%position,-Para_SM%tilt(2,Para_SM%numero_miroir))
                write(31) photon_M%position
            else
                write(31) photon_M%position
            endif
        endif

end subroutine parcoursphoton


end module PHOTON
