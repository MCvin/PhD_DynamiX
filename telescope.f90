module TELESCOPE

use OUTILS
use LECTURE

implicit none


contains


!------------------------------------------------------------------------
!-          Exprime la position et l'attitude des instruments           -
!-                  (Miroirs, plan de Detection)                        -
!-           dans leur repere de reference (RefPM et RefPD)             -
!-                  (position et attitude relatives)                    -
!------------------------------------------------------------------------
subroutine instruments(Para_SM,Para_SD,etat, RefM,RefD)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM
        type (ParametresSatDetecteur) :: Para_SD
        integer :: etat

        type (Repere) :: RefM         !referentiel des miroirs
        type (Repere) :: RefD         !referentiel des detecteurs

! Local variables
        real (kind=wp), dimension(3) :: deriveSM, deriveSD


! geometrie du satellite miroir
            !derive en attitude dans le referentiel de pointage
            deriveSM=Para_SM%mouvement(etat)%attitude

    !module miroir
        !attitude et position du module miroir (RefM) dans le referentiel du satellite (RefSM)
            !l'attitude relative du module miroir est egale aux derives (relatives)
            RefM%attitude=deriveSM

            !position du module miroir (RefM) dans RefSM
            RefM%position = Para_SM%pos_miroirs    !ici le module miroir est centre sur le centre de masse

            !position de RefM dans RefPM
            !si le module miroir n'est pas centre sur le centre de masse, sa position est affectee
            !par la rotation du satellite; mais ceci est inutile dans la mesure ou l'on va exprimer
            !la position des detecteurs dans RefSM et pas dans RefPM


! geometrie du satellite detecteur
            !derive en attitude par rapport au referentiel de pointage
            deriveSD=Para_SD%mouvement(etat)%attitude

    !detecteurs
        !attitude et position du plan de detection (RefD) dans le referentiel de pointage (RefPD)
            !l'attitude relative du FPA est egale aux derives (relatives)
            RefD%attitude=deriveSD

            !position des detecteurs (RefD) dans RefSD (on prend comme position du FPA la position du LED)
            RefD%position = Para_SD%pos_LED

            !position de RefD dans RefPD
            !si les detecteurs ne sont pas centre sur le centre de masse, leur position est affectee
            !par la rotation du satellite
            call rotationZ(RefD%position,deriveSD(3))
            call rotationY(RefD%position,deriveSD(2))
            call rotationX(RefD%position,deriveSD(1))

        !attitude et position du plan de detection (RefD) dans referentiel des miroirs (RefM)
            !position de RefD dans RefPM
            RefD%position = RefD%position + Para_SD%mouvement(etat)%position

            !position de RefD dans RefSM
            call rotationZ(RefD%position,-deriveSM(3))
            call rotationY(RefD%position,-deriveSM(2))
            call rotationX(RefD%position,-deriveSM(1))

            !position de RefD dans RefM
            RefD%position = RefD%position - RefM%position


end subroutine instruments


!------------------------------------------------------------------------
!-                      Programme principal                             -
!------------------------------------------------------------------------
subroutine etat_telescope(Para_SM,Para_SD,etat, RefM,RefD)

implicit none

! IN/OUT variables
        type (ParametresSatMiroir) :: Para_SM
        type (ParametresSatDetecteur) :: Para_SD
        integer :: etat

        type (Repere) :: RefM         !referentiel des miroirs
        type (Repere) :: RefD         !referentiel des detecteurs


! Position et Attitude des Instruments a bord (miroirs, detecteurs, senseurs)
        call instruments(Para_SM,Para_SD,etat, RefM,RefD)


end subroutine etat_telescope


end module TELESCOPE
