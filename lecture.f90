module LECTURE

implicit none

    integer, parameter :: wp = kind(1.d0)                                !precision des nombres reels du code

    type parabole
       real (kind=wp) :: parametre
       real (kind=wp) :: origine
       real (kind=wp) :: focale
    endtype parabole
    type hyperbole
       real (kind=wp) :: parametre
       real (kind=wp) :: origine
       real (kind=wp) :: focale
       real (kind=wp) :: e2
    endtype hyperbole
    type Repere
       real (kind=wp), dimension(3) :: position
       real (kind=wp), dimension(3) :: attitude
    endtype Repere
    type photons
       real (kind=wp), dimension(3) :: position
       real (kind=wp), dimension(3) :: direction
       real (kind=wp) :: energie
    endtype photons

        type Source
            real (kind=wp), dimension(2) :: coordonnees                  !coordonnees equatoriales de la source
            real (kind=wp) :: luminosite                                 !luminosite de la source
            real (kind=wp) :: indice                                     !indice spectral de la source

        endtype Source

        type Observation
            logical :: saves                                             !activation des sauvegardes
            real (kind=wp), dimension(2) :: pointage                     !pointage de reference du telescope
            integer :: nb_sources                                        !nombre de sources simulees
            logical :: luminosites_sources                               !activation de la luminosite propre a chaque source
            type (Source), dimension(1000) :: sources                    !parametres des sources
            integer :: nb_photons                                        !nombre total de photons simules
            logical :: energie_unique                                    !mono-energie des photons
            real (kind=wp) :: energie_value                              !valeur de l'energie des photons
            logical :: energie_uniforme                                  !distribution uniforme de l'energie photonique
            logical :: energie_powerlaw                                  !distribution en powerlaw de l'energie photonique
            real (kind=wp), dimension(2) :: minmax_energie               !bande en energie de l'instrument
            real (kind=wp), dimension(:,:), allocatable :: coeff_reflec  !tableau contenant les coefficients de reflection
            integer :: seed                                              !graine pour la statistique

        endtype Observation

        type ParametresSatMiroir
            type (Repere) :: pointage                                    !pointage de reference du telescope
            real (kind=wp), dimension(3) :: pos_miroirs                  !position du module miroir sur le MSC
            real (kind=wp), dimension(2,100) :: tilt                     !erreur d'alignement de chaque miroir
            real (kind=wp), dimension(3,100) :: delta                    !erreur de positionnement de chaque miroir
            real (kind=wp)  :: forme                                     !defaut de forme des miroirs
            real (kind=wp)  :: spread                                    !defaut de surface des miroirs
            real (kind=wp), dimension(3,3) :: pos_cibles                 !position des cibles sur le MSC
            logical :: STR                                               !activation des defauts STR
            real (kind=wp), dimension(3) :: biais_sens_att               !biais du STR
            real (kind=wp) :: noise_sens_att                             !noise du STR
            logical, dimension(3) :: derive_att                          !activation des derives du MSC
            integer :: nb_etats                                          !nombre d'etats de la formation
            type (Repere), dimension(:), allocatable :: mouvement        !position et attitude du MSC
            integer :: nb_miroirs                                        !nombre de miroirs
            real (kind=wp)  :: rayon_max                                 !rayon du plus grand shell
            logical :: entreemiroirs                                     !le photon est-il rentre dans les miroirs
            integer :: numero_miroir                                     !quel miroir le photon va t-il rencontrer
            real (kind=wp), dimension(:), allocatable :: rayon           !rayon des shells au niveau de la face d'entree
            type (parabole), dimension(:), allocatable :: paraboles      !parametres du profil parabolique des miroirs
            type (hyperbole), dimension(:), allocatable :: hyperboles    !parametres du profil hyperbolique des miroirs
            logical :: reflexionparabole                                 !y a t-il eu reflexion sur la partie parabolique
            logical :: reflexionhyperbole                                !y a t-il eu reflexion sur la partie hyperbolique
            logical :: assembly                                          !activation des defauts d'assemblage des miroirs
            logical :: figure                                            !activation des defauts de forme des miroirs
            logical :: scattering                                        !activation des defauts de surface des miroirs
            logical :: reflexiontotale                                   !reflexion totale ou selon des coefficients de reflexion
            real (kind=wp), dimension(3) :: err_Att                      !erreur d'attitude du MSC (calculee avec les donnees senseurs)

        endtype ParametresSatMiroir


        type ParametresSatDetecteur
            type (Repere) :: pointage                                    !pointage de reference du telescope
            real (kind=wp) :: h_col                                      !hauteur du collimateur
            real (kind=wp) :: r_col                                      !rayon du collimateur
            real (kind=wp), dimension(3) :: pos_LED                      !position du LED sur le DSC
            real (kind=wp), dimension(3) :: pos_HED                      !position du HED sur le DSC
            real (kind=wp), dimension(3) :: pos_sens                     !position du senseur optique sur le DSC
            logical :: STR                                               !activation des defauts STR
            real (kind=wp), dimension(3) :: biais_sens_att               !biais du STR
            real (kind=wp) :: noise_sens_att                             !noise du STR
            logical :: CLM                                               !activation du CLM
            real (kind=wp), dimension(3) :: biais_sens_pos               !biais du CLM
            real (kind=wp) :: noise_sens_pos                             !noise du CLM
            logical, dimension(3) :: derive_pos, derive_att              !activation des derives du DSC
            integer :: nb_etats                                          !nombre d'etats de la formation
            real (kind=wp) :: defocus                                    !defocalisation forcee du DSC
            type (Repere), dimension(:), allocatable :: mouvement        !position et attitude du DSC
            logical :: LED                                               !activation du LED
            logical :: HED                                               !activation du HED
            real (kind=wp), dimension(2) :: dimensions_LED               !dimensions du LED
            real (kind=wp), dimension(2) :: dimensions_HED               !dimensions du HED
            real (kind=wp), dimension(2) :: nb_pixels_LED                !nombre de pixels du LED
            real (kind=wp), dimension(2) :: nb_pixels_HED                !nombre de pixels du HED
            real (kind=wp) :: taille_pixels_LED                          !taille des pixels du LED
            real (kind=wp) :: taille_pixels_HED                          !taille des pixels du HED
            integer, dimension(135,135) :: HED_map                       !geometrie du HED
            real (kind=wp), dimension(3) :: err_Att                      !erreur d'attitude du DSC (calculee avec les donnees senseurs)
            real (kind=wp), dimension(3) :: err_Pos_LED                  !erreur de positionnement du LED du a la formation (calculee avec les donnees senseurs)
            real (kind=wp), dimension(3) :: err_Pos_HED                  !erreur de positionnement du HED du a la formation (calculee avec les donnees senseurs)

        endtype ParametresSatDetecteur


contains

!------------------------------------------------------------------------
!-                      Programme principal                             -
!------------------------------------------------------------------------
subroutine lectureparametres(Obs,Para_SM,Para_SD)

implicit none

!parametres internes
        character(len=100) :: chaine, mot, switch
        character(len=100) :: file_mirrors_para, file_mirrors_hyper, file_mirrors_entree
        character(len=100) :: file_coeffreflection
        character(len=100) :: file_mouvements_MSC, file_mouvements_DSC
        character(len=100) :: file_PSF_para, file_PSF_reelle, file_PSF_corrige
        character(len=100) :: file_spec_LED_total, file_spec_HED_total
        character(len=100) :: file_spec_LED_pixels, file_spec_HED_pixels
        character(len=100) :: file_spec_LED_pixels_corr, file_spec_HED_pixels_corr

        integer :: err, debut, fin, debut1, fin1, debut2, fin2, n
        type (Repere), dimension(:), allocatable :: etats_MSC, etats_DSC
        real (kind=wp) :: tilt_param
        real (kind=wp), dimension(3) :: delta_param
        integer :: i
        real (kind=wp) :: G05DDF

!parametres observation
        type (Observation) :: Obs

!parametres satellite miroirs
        type (ParametresSatMiroir) :: Para_SM

!parametres satellite detecteurs
        type (ParametresSatDetecteur) :: Para_SD

!configuration de la simulation
    open(1,file='switch.par')

    err=0
    do while (err.eq.0)

       read(1,'(A)',iostat=err) chaine

       if ((chaine(1:1).ne.'#') .and. (err.eq.0)) then
          mot=chaine(1:(index(chaine,':'))-1)
          debut=index(chaine,'[')+1
          fin=index(chaine,']')-1
          switch=chaine(debut:fin)

          select case (mot)
            !lecture des parametres de l'observation
             case ('sauvegarde des donnees')                
                if (switch.eq.'YES') then
                    Obs%saves=.TRUE.
                else
                    Obs%saves=.FALSE.
                endif
             case ('luminosite propre a chaque source')                
                if (switch.eq.'YES') then
                    Obs%luminosites_sources=.TRUE.
                else
                    Obs%luminosites_sources=.FALSE.
                endif
             case ('mono-energie des photons')
                if (switch.eq.'YES') then
                    Obs%energie_unique=.TRUE.
                else
                    Obs%energie_unique=.FALSE.
                endif
             case ("distribution uniforme de l'energie photonique")
                if (switch.eq.'YES') then
                    Obs%energie_uniforme=.TRUE.
                else
                    Obs%energie_uniforme=.FALSE.
                endif
             case ("distribution en powerlaw de l'energie photonique")
                if (switch.eq.'YES') then
                    Obs%energie_powerlaw=.TRUE.
                else
                    Obs%energie_powerlaw=.FALSE.
                endif

            !lecture des parametres du satellite miroir
             case ("derive de l'attitude du MSC selon x")                
                if (switch.eq.'YES') then
                    Para_SM%derive_att(1)=.TRUE.
                else
                    Para_SM%derive_att(1)=.FALSE.
                endif
             case ("derive de l'attitude du MSC selon y")                
                if (switch.eq.'YES') then
                    Para_SM%derive_att(2)=.TRUE.
                else
                    Para_SM%derive_att(2)=.FALSE.
                endif
             case ("derive de l'attitude du MSC selon z")                
                if (switch.eq.'YES') then
                    Para_SM%derive_att(3)=.TRUE.
                else
                    Para_SM%derive_att(3)=.FALSE.
                endif
             case ("defauts d'assemblage des miroirs")
                if (switch.eq.'YES') then
                    Para_SM%assembly=.TRUE.
                else
                    Para_SM%assembly=.FALSE.
                endif
             case ('defauts de forme des miroirs')
                if (switch.eq.'YES') then
                    Para_SM%figure=.TRUE.
                else
                    Para_SM%figure=.FALSE.
                endif
             case ('defauts de surface des miroirs')
                if (switch.eq.'YES') then
                    Para_SM%scattering=.TRUE.
                else
                    Para_SM%scattering=.FALSE.
                endif
             case ('defauts STR du MSC')
                if (switch.eq.'YES') then
                    Para_SM%STR=.TRUE.
                else
                    Para_SM%STR=.FALSE.
                endif
             case ('reflection totale')
                if (switch.eq.'YES') then
                    Para_SM%reflexiontotale=.TRUE.
                else
                    Para_SM%reflexiontotale=.FALSE.
                endif

            !lecture des parametres du satellite detecteur
             case ('derive de la position du DSC selon x')                
                if (switch.eq.'YES') then
                    Para_SD%derive_pos(1)=.TRUE.
                else
                    Para_SD%derive_pos(1)=.FALSE.
                endif
             case ('derive de la position du DSC selon y')                
                if (switch.eq.'YES') then
                    Para_SD%derive_pos(2)=.TRUE.
                else
                    Para_SD%derive_pos(2)=.FALSE.
                endif
             case ('derive de la position du DSC selon z')                
                if (switch.eq.'YES') then
                    Para_SD%derive_pos(3)=.TRUE.
                else
                    Para_SD%derive_pos(3)=.FALSE.
                endif
             case ("derive de l'attitude du DSC selon x")                
                if (switch.eq.'YES') then
                    Para_SD%derive_att(1)=.TRUE.
                else
                    Para_SD%derive_att(1)=.FALSE.
                endif
             case ("derive de l'attitude du DSC selon y")                
                if (switch.eq.'YES') then
                    Para_SD%derive_att(2)=.TRUE.
                else
                    Para_SD%derive_att(2)=.FALSE.
                endif
             case ("derive de l'attitude du DSC selon z")                
                if (switch.eq.'YES') then
                    Para_SD%derive_att(3)=.TRUE.
                else
                    Para_SD%derive_att(3)=.FALSE.
                endif
             case ('defauts STR du DSC')
                if (switch.eq.'YES') then
                    Para_SD%STR=.TRUE.
                else
                    Para_SD%STR=.FALSE.
                endif
             case ('defauts CLM du DSC')
                if (switch.eq.'YES') then
                    Para_SD%CLM=.TRUE.
                else
                    Para_SD%CLM=.FALSE.
                endif
             case ('Low Energie Detector')
                if (switch.eq.'YES') then
                    Para_SD%LED=.TRUE.
                else
                    Para_SD%LED=.FALSE.
                endif
             case ('High Energie Detector')
                if (switch.eq.'YES') then
                    Para_SD%HED=.TRUE.
                else
                    Para_SD%HED=.FALSE.
                endif

          end select

       end if

    end do
    close(1)

print*
print*, 'Sauvegarde donnees:', Obs%saves
print*, 'Luminositees des sources:', Obs%luminosites_sources
print*, "Valeur unique de l'energie des photons", Obs%energie_unique
print*, "Distribution de l'energie des photons,"
print*, 'loi uniforme:', Obs%energie_uniforme, 'power law:', Obs%energie_powerlaw
print*
print*, 'Derive attitude MSC:', Para_SM%derive_att
print*, "Defauts d'assemblage des miroirs:", Para_SM%assembly
print*, "Defauts de forme des miroirs:", Para_SM%figure
print*, 'Defauts de surface des miroirs:',Para_SM%scattering
print*, 'Defauts STR du MSC:', Para_SM%STR
print*, 'Reflection totale', Para_SM%reflexiontotale
print*
print*, 'Derive position DSC:', Para_SD%derive_pos
print*, 'Derive attitude DSC:', Para_SD%derive_att
print*, 'Defauts STR du DSC:', Para_SD%STR
print*, 'Defauts CLM du DSC:', Para_SD%CLM
print*
print*, 'Detecteur basse energie:', Para_SD%LED
print*, 'Detecteur haute energie:', Para_SD%HED


!initialisation des parametres par lecture du fichier parametres.par
    open(1,file='parametres.par')

    n=1
    err = 0
    do while (err.eq.0)

       read(1,'(A)',iostat=err) chaine

       if ((chaine(1:1).ne.'#') .and. (err.eq.0)) then
          mot=chaine(1:(index(chaine,':'))-1)
          debut=index(chaine,'[')+1
          fin=index(chaine,']')-1
          debut1=index(chaine,'(')+1
          fin1=index(chaine,')')-1
          debut2=index(chaine,'{')+1
          fin2=index(chaine,'}')-1

          select case (mot)

            !lecture des parametres de l'observation
             case ('pointage')
                read(chaine(debut:fin),*) Obs%pointage
                !conversion des arcmin du fichier parametre.par en degrees pour le code
                Obs%pointage=Obs%pointage/60._wp
             case ('source')
                read(chaine(debut:fin),*) Obs%sources(n)%coordonnees
                !conversion des arcmin du fichier parametre.par en degrees pour le code
                Obs%sources(n)%coordonnees=Obs%sources(n)%coordonnees/60._wp
                read(chaine(debut1:fin1),*) Obs%sources(n)%luminosite
                read(chaine(debut2:fin2),*) Obs%sources(n)%indice
                n=n+1
                Obs%nb_sources=n-1
             case ('nombre de photons')
                read(chaine(debut:fin),*) Obs%nb_photons
             case ('bande en energie')
                read(chaine(debut:fin),*) Obs%minmax_energie
             case ('energie fixe')
                read(chaine(debut:fin),*) Obs%energie_value
             case ('graine')
                read(chaine(debut:fin),*) Obs%seed

            !lecture des parametres du satellite miroir
             case ('nombres de miroirs')
                read(chaine(debut:fin),*) Para_SM%nb_miroirs
             case ('position du module miroir')
                read(chaine(debut:fin),*) Para_SM%pos_miroirs
             case ("erreur d'alignement des miroirs")
                read(chaine(debut:fin),*) tilt_param
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                tilt_param=tilt_param/3600.0
             case ('erreur de positionnement des miroirs')
                read(chaine(debut:fin),*) delta_param
             case ('erreur de forme des miroirs')
                read(chaine(debut:fin),*) Para_SM%forme
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SM%forme=Para_SM%forme/3600._wp
             case ('polissage des miroirs')
                read(chaine(debut:fin),*) Para_SM%spread
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SM%spread=Para_SM%spread/3600._wp
             case ('position de la cible 1')
                read(chaine(debut:fin),*) Para_SM%pos_cibles(1,:)
             case ('position de la cible 2')
                read(chaine(debut:fin),*) Para_SM%pos_cibles(2,:)
             case ('position de la cible 3')
                read(chaine(debut:fin),*) Para_SM%pos_cibles(3,:)
             case ('biais de la mesure du STR du MSC')
                read(chaine(debut:fin),*) Para_SM%biais_sens_att
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SM%biais_sens_att=Para_SM%biais_sens_att/3600._wp
             case ('bruit pixel de la mesure du STR du MSC')
                read(chaine(debut:fin),*) Para_SM%noise_sens_att
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SM%noise_sens_att=Para_SM%noise_sens_att/3600._wp

            !lecture des parametres du satellite detecteur
             case ('hauteur du collimateur')
                read(chaine(debut:fin),*) Para_SD%h_col
             case ('rayon du collimateur')
                read(chaine(debut:fin),*) Para_SD%r_col
             case ('position du LED')
                read(chaine(debut:fin),*) Para_SD%pos_LED
             case ('position du HED')
                read(chaine(debut:fin),*) Para_SD%pos_HED
             case ('position du CLM du DSC')
                read(chaine(debut:fin),*) Para_SD%pos_sens
             case ('biais de la mesure du STR du DSC')
                read(chaine(debut:fin),*) Para_SD%biais_sens_att
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SD%biais_sens_att=Para_SD%biais_sens_att/3600._wp
             case ('bruit pixel de la mesure du STR du DSC')
                read(chaine(debut:fin),*) Para_SD%noise_sens_att
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SD%noise_sens_att=Para_SD%noise_sens_att/3600._wp
             case ('biais de la mesure du CLM du DSC')
                read(chaine(debut:fin),*) Para_SD%biais_sens_pos
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SD%biais_sens_pos=Para_SD%biais_sens_pos/3600._wp
             case ('bruit pixel de la mesure du CLM du DSC')
                read(chaine(debut:fin),*) Para_SD%noise_sens_pos
                !conversion des arcsec du fichier parametre.par en degrees pour le code
                Para_SD%noise_sens_pos=Para_SD%noise_sens_pos/3600._wp
             case ('dimensions du LED')
                read(chaine(debut:fin),*) Para_SD%dimensions_LED
             case ('dimensions du HED')
                read(chaine(debut:fin),*) Para_SD%dimensions_HED
             case ('nombre de pixels du LED')
                read(chaine(debut:fin),*) Para_SD%nb_pixels_LED
             case ('nombre de pixels du HED')
                read(chaine(debut:fin),*) Para_SD%nb_pixels_HED
             case ('taille des pixels du LED')
                read(chaine(debut:fin),*) Para_SD%taille_pixels_LED
             case ('taille des pixels du HED')
                read(chaine(debut:fin),*) Para_SD%taille_pixels_HED

             !lecture des parametres du vol en formation
             case ("nombre d'etats de la formation")
                read(chaine(debut:fin),*) Para_SM%nb_etats
                Para_SD%nb_etats=Para_SM%nb_etats
             case ("defocalisation")
                read(chaine(debut:fin),*) Para_SD%defocus

             !lecture des fichiers d'entree
             case ("fichier parametre face d'entree des miroirs")
                read(chaine(debut:fin),'(A)') file_mirrors_entree
             case ("fichier parametre des miroirs paraboliques")
                read(chaine(debut:fin),'(A)') file_mirrors_para
             case ("fichier parametre des miroirs hyperboliques")
                read(chaine(debut:fin),'(A)') file_mirrors_hyper
             case ("fichier parametre des coefficients de reflexion")
                read(chaine(debut:fin),'(A)') file_coeffreflection
             case ("fichier donnees des mouvements du MSC")
                read(chaine(debut:fin),'(A)') file_mouvements_MSC
             case ("fichier donnees des mouvements du DSC")
                read(chaine(debut:fin),'(A)') file_mouvements_DSC

             !lecture des fichiers de sortie   
             case ("fichier parametre PSF")
                read(chaine(debut:fin),'(A)') file_PSF_para
             case ("fichier donnees PSF reelle")
                read(chaine(debut:fin),'(A)') file_PSF_reelle
             case ("fichier donnees PSF corrigee")
                read(chaine(debut:fin),'(A)') file_PSF_corrige
             case ("fichier donnees spectre LED total")
                read(chaine(debut:fin),'(A)') file_spec_LED_total
             case ("fichier donnees spectre HED total")
                read(chaine(debut:fin),'(A)') file_spec_HED_total
             case ("fichier donnees image et spectre LED")
                read(chaine(debut:fin),'(A)') file_spec_LED_pixels
             case ("fichier donnees image et spectre HED")
                read(chaine(debut:fin),'(A)') file_spec_HED_pixels
             case ("fichier donnees image corrigee et spectre LED")
                read(chaine(debut:fin),'(A)') file_spec_LED_pixels_corr
             case ("fichier donnees image corrigee et spectre HED")
                read(chaine(debut:fin),'(A)') file_spec_HED_pixels_corr

          end select

       end if

    end do
    close(1)

print*
print*, 'Nombres de sources dans le champ de vue:', Obs%nb_sources
print*, 'Energie:'
if (Obs%energie_unique) then
print*, Obs%energie_value, 'keV'
else if (Obs%energie_uniforme .or. Obs%energie_powerlaw) then
print*, Obs%minmax_energie, 'keV'
endif
print*, "Nombre d'etats de la formation:",Para_SM%nb_etats
print*, 'Nombre de photons envoyes:',Obs%nb_photons
print*, 'Graine utilisee pour la statistique:',Obs%seed

!initialisation de la graine pour les fonctions random du code
    call G05CBF(Obs%seed)

!lecture des profils des miroirs
    allocate(Para_SM%rayon(Para_SM%nb_miroirs))
    allocate(Para_SM%paraboles(Para_SM%nb_miroirs))
    allocate(Para_SM%hyperboles(Para_SM%nb_miroirs))
    !on rentre les parametres miroirs dans des tableaux
    open(1,file=file_mirrors_entree)
    open(2,file=file_mirrors_para)
    open(3,file=file_mirrors_hyper)
    read(1,*) Para_SM%rayon
    read(2,*) Para_SM%paraboles
    read(3,*) Para_SM%hyperboles
    close(1)
    close(2)
    close(3)

!initialisation des defauts d'alignement et de positionnement des miroirs
    if (Para_SM%assembly) then
        do i=1,Para_SM%nb_miroirs
            Para_SM%tilt(1,i)=G05DDF(0._wp,tilt_param)
            Para_SM%tilt(2,i)=G05DDF(0._wp,tilt_param)
            Para_SM%delta(1,i)=G05DDF(0._wp,delta_param(1))
            Para_SM%delta(2,i)=G05DDF(0._wp,delta_param(2))
            Para_SM%delta(3,i)=G05DDF(0._wp,delta_param(3))
        enddo
    endif

!lecture des coefficients de reflexion
    if ((Obs%energie_unique) .or. (Obs%energie_uniforme) .or. (Obs%energie_powerlaw)) then
        allocate(Obs%coeff_reflec(1000,1000))
        !on rentre les coefficients de reflexion dans un tableau
        open(1,file=file_coeffreflection)
        read(1,*) Obs%coeff_reflec
        close(1)
        !on transpose les lignes et colonnes car a la lecture le sens de la matrice change: l'energie devient ligne et l'incidence devient colonne
        Obs%coeff_reflec=TRANSPOSE(Obs%coeff_reflec)
    endif

!definition de la geometrie du HED
    Para_SD%HED_map=1
    !pixels morts: 1 pixel mort entre chaque module caliste de 16*16 pixels
        !en vertical
        Para_SD%HED_map(:,17)=0
        Para_SD%HED_map(:,34)=0
        Para_SD%HED_map(:,51)=0
        Para_SD%HED_map(:,68)=0
        Para_SD%HED_map(:,85)=0
        Para_SD%HED_map(:,102)=0
        Para_SD%HED_map(:,119)=0
        !en horizontal
        Para_SD%HED_map(17,:)=0
        Para_SD%HED_map(34,:)=0
        Para_SD%HED_map(51,:)=0
        Para_SD%HED_map(68,:)=0
        Para_SD%HED_map(85,:)=0
        Para_SD%HED_map(102,:)=0
        Para_SD%HED_map(119,:)=0


!initialisation des reperes de reference (lies au pointage de la source)
    ! =>les positions sont exprimees en relatif / au MSC
    ! =>l'attitude est exprimee en absolu en cartesien dans le ciel J2000
    !du satellite miroir
        Para_SM%pointage%position = (/ 0._wp , 0._wp , 0._wp /)

        Para_SM%pointage%attitude(1)=Obs%pointage(1)    !alpha
        Para_SM%pointage%attitude(2)=Obs%pointage(2)    !delta
        Para_SM%pointage%attitude(3)=0._wp              !0 car on suppose une source isotrope

    !du satellite detecteur
        Para_SD%pointage%position = (/ 0._wp , 0._wp , 0._wp /)
        if (Para_SD%LED) then  !on focalise sur le LED
           Para_SD%pointage%position(3)=-Para_SM%hyperboles(1)%focale &
                                        +Para_SM%pos_miroirs(3)       & !pour compenser les miroirs en aval
                                        -Para_SD%pos_LED(3)           & !pour compenser les detecteurs qui sont en amont du plan focal
                                        -Para_SD%defocus                !pour simuler une defocalisation (sans mouvements)
        else if (Para_SD%HED) then  !on focalise sur le HED
           Para_SD%pointage%position(3)=-Para_SM%hyperboles(1)%focale &
                                        +Para_SM%pos_miroirs(3)       & !pour compenser les miroirs en aval
                                        -Para_SD%pos_HED(3)           & !pour compenser les detecteurs qui sont en amont du plan focal
                                        -Para_SD%defocus                !pour simuler une defocalisation (sans mouvements)
        endif
        Para_SD%pointage%attitude(1)=Obs%pointage(1)    !alpha
        Para_SD%pointage%attitude(2)=Obs%pointage(2)    !delta
        Para_SD%pointage%attitude(3)=0._wp              !0 car on suppose une source isotrope 


!lecture des mouvements des satellites: les mouvements sont exprimes en relatif / aux reperes de reference
    !on lit les etats de la formation
        allocate(etats_MSC(Para_SM%nb_etats))
        open(1,file=file_mouvements_MSC)
        read(1,*) etats_MSC
        close(1)
        allocate(etats_DSC(Para_SD%nb_etats))
        open(1,file=file_mouvements_DSC)
        read(1,*) etats_DSC
        close(1)
!print*,maxval(etats_MSC(:)%attitude(1))*3600*3*20
!!$print*,maxval(etats_DSC%position(3))+20.396
!!$print*,minval(etats_DSC%position(3))+20.396
!!$print*,maxval(etats_DSC%attitude(1)),maxval(etats_DSC%attitude(2)),maxval(etats_DSC%attitude(3))
!!$print*,minval(etats_DSC%attitude(1)),minval(etats_DSC%attitude(2)),minval(etats_DSC%attitude(3))
!!$print*,maxval(etats_MSC%attitude(1)),maxval(etats_MSC%attitude(2)),maxval(etats_MSC%attitude(3))
!!$print*,minval(etats_MSC%attitude(1)),minval(etats_MSC%attitude(2)),minval(etats_MSC%attitude(3))
    !par defaut les satellites restent immobiles
        allocate(Para_SM%mouvement(Para_SM%nb_etats))
        allocate(Para_SD%mouvement(Para_SD%nb_etats))
        ! =>la positition est exprimee en relatif / au MSC
        Para_SM%mouvement(:)%position(1) = Para_SM%pointage%position(1)
        Para_SM%mouvement(:)%position(2) = Para_SM%pointage%position(2)
        Para_SM%mouvement(:)%position(3) = Para_SM%pointage%position(3)
        ! =>l'attitude est exprimee en relatif pour ne pas dependre du pointage
        Para_SM%mouvement(:)%attitude(1) = 0._wp !+ Obs%minmax_energie(1)/60._wp
        Para_SM%mouvement(:)%attitude(2) = 0._wp !+ Obs%minmax_energie(1)/60._wp
        Para_SM%mouvement(:)%attitude(3) = 0._wp !+ Obs%minmax_energie(1)/60._wp
        ! =>la positition est exprimee en relatif / au MSC
        Para_SD%mouvement(:)%position(1) = Para_SD%pointage%position(1)
        Para_SD%mouvement(:)%position(2) = Para_SD%pointage%position(2)
        Para_SD%mouvement(:)%position(3) = Para_SD%pointage%position(3)
        ! =>l'attitude est exprimee en relatif pour ne pas dependre du pointage
        Para_SD%mouvement(:)%attitude(1) = 0._wp !+ Obs%energie_value!+ 20/3600._wp
        Para_SD%mouvement(:)%attitude(2) = 0._wp !+ Obs%energie_value!/3600._wp
        Para_SD%mouvement(:)%attitude(3) = 0._wp !+ Obs%energie_value!+ Obs%energie_value!+ 20/3600._wp

    !mouvements du satellite miroir
        ! =>l'attitude est exprimee en relatif pour ne pas dependre du pointage
        if (Para_SM%derive_att(1)) then
            Para_SM%mouvement(:)%attitude(1) = Para_SM%mouvement(:)%attitude(1) + etats_MSC(:)%attitude(1)!*3*Obs%minmax_energie(1)!*3*20
        end if
        if (Para_SM%derive_att(2)) then
            Para_SM%mouvement(:)%attitude(2) = Para_SM%mouvement(:)%attitude(2) + etats_MSC(:)%attitude(2)!*3!*Obs%minmax_energie(1)!*3*20
        end if
        if (Para_SM%derive_att(3)) then
            Para_SM%mouvement(:)%attitude(3) = Para_SM%mouvement(:)%attitude(3) + etats_MSC(:)%attitude(3)!*3!*Obs%minmax_energie(1)!*3*20
        end if

print*,maxval(Para_SM%mouvement%attitude(1)),maxval(Para_SM%mouvement%attitude(2)),maxval(Para_SM%mouvement%attitude(3))
print*,minval(Para_SM%mouvement%attitude(1)),minval(Para_SM%mouvement%attitude(2)),minval(Para_SM%mouvement%attitude(3))

    !mouvements du satellite detecteur
        ! =>la positition est exprimee en relatif / au MSC
        if (Para_SD%derive_pos(1)) then
            Para_SD%mouvement(:)%position(1) = Para_SD%mouvement(:)%position(1) + etats_DSC(:)%position(1)
        end if
        if (Para_SD%derive_pos(2)) then
            Para_SD%mouvement(:)%position(2) = Para_SD%mouvement(:)%position(2) + etats_DSC(:)%position(2)
        end if
        if (Para_SD%derive_pos(3)) then
            Para_SD%mouvement(:)%position(3) = Para_SD%mouvement(:)%position(3) + etats_DSC(:)%position(3)
        end if
        ! =>l'attitude est exprimee en relatif pour ne pas dependre du pointage
        if (Para_SD%derive_att(1)) then
            Para_SD%mouvement(:)%attitude(1) = Para_SD%mouvement(:)%attitude(1) + etats_DSC(:)%attitude(1)
        end if
        if (Para_SD%derive_att(2)) then
            Para_SD%mouvement(:)%attitude(2) = Para_SD%mouvement(:)%attitude(2) + etats_DSC(:)%attitude(2)
        end if
        if (Para_SD%derive_att(3)) then
            Para_SD%mouvement(:)%attitude(3) = Para_SD%mouvement(:)%attitude(3) + etats_DSC(:)%attitude(3)
        end if


!ouverture des fichiers de sortie pour preparer a l'ecriture
    open(10,file=file_PSF_para)
    open(11,file=file_PSF_reelle,form='binary')
    open(12,file=file_PSF_corrige,form='binary')
    open(13,file=file_spec_LED_total)
    open(14,file=file_spec_HED_total)
    open(15,file=file_spec_LED_pixels)
    open(16,file=file_spec_HED_pixels)
    open(17,file=file_spec_LED_pixels_corr)
    open(18,file=file_spec_HED_pixels_corr)

    !fichiers temporaires
!!$    open(21,file='output/CLM1.data')
!!$    open(22,file='output/CLM2.data')
!!$    open(23,file='output/CLM3.data')
!!$    open(24,file='output/CLM4.data')
!!$    open(25,file='output/CLM5.data')
    open(31,file='output/pos.data',form='binary')
    open(32,file='output/incidences.data',form='binary')


end subroutine lectureparametres


end module LECTURE
