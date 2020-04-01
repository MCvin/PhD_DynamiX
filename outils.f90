module OUTILS

implicit none

    integer, parameter :: vp = kind(1.d0)


contains


!------------------------------------------------------------------------
!-      Outil de calcul: Matrices de rotations (sens trigo)             -
!------------------------------------------------------------------------
subroutine rotationX(vecteur,ax)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: vecteur
        real (kind=vp) :: ax

! local variables
        real (kind=vp), dimension(3,3) :: Rx

!initialisation de la matrice de rotation autour de l'axe X

        Rx(1,:)=(/  1._vp ,       0._vp ,      0._vp  /)
        Rx(2,:)=(/  0._vp ,    COSD(ax) ,  -SIND(ax)  /)
        Rx(3,:)=(/  0._vp ,    SIND(ax) ,   COSD(ax)  /)

        vecteur=MATMUL(Rx,vecteur)

end subroutine rotationX


subroutine rotationY(vecteur,ay)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: vecteur
        real (kind=vp) :: ay

! local variables
        real (kind=vp), dimension(3,3) :: Ry

!initialisation de la matrice de rotation autour de l'axe Y

        Ry(1,:)=(/  COSD(ay) ,   0._vp ,   SIND(ay)  /)
        Ry(2,:)=(/     0._vp ,   1._vp ,      0._vp  /)
        Ry(3,:)=(/ -SIND(ay) ,   0._vp ,   COSD(ay)  /)

        vecteur=MATMUL(Ry,vecteur)

end subroutine rotationY


subroutine rotationZ(vecteur,az)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: vecteur
        real (kind=vp) :: az

! local variables
        real (kind=vp), dimension(3,3) :: Rz

!initialisation de la matrice de rotation autour de l'axe Z

        Rz(1,:)=(/  COSD(az) ,  -SIND(az) ,   0._vp  /)
        Rz(2,:)=(/  SIND(az) ,   COSD(az) ,   0._vp  /)
        Rz(3,:)=(/     0._vp ,      0._vp ,   1._vp  /)

        vecteur=MATMUL(Rz,vecteur)

end subroutine rotationZ


subroutine rotation_deg(axe,vecteur,a)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: axe
        real (kind=vp), dimension(3) :: vecteur
        real (kind=vp) :: a

! local variables
        real (kind=vp), dimension(3,3) :: R

!l'axe doit etre normalise
        axe=axe/norme(axe)

!initialisation de la matrice de rotation autour de l'axe Z

        R(1,1)=axe(1)**2._vp+COSD(a)*(1._vp-axe(1)**2._vp)
        R(1,2)=axe(1)*axe(2)*(1._vp-COSD(a))-axe(3)*SIND(a)
        R(1,3)=axe(1)*axe(3)*(1._vp-COSD(a))+axe(2)*SIND(a)
        R(2,1)=axe(1)*axe(2)*(1._vp-COSD(a))+axe(3)*SIND(a)
        R(2,2)=axe(2)**2._vp+COSD(a)*(1._vp-axe(2)**2._vp)
        R(2,3)=axe(2)*axe(3)*(1._vp-COSD(a))-axe(1)*SIND(a)
        R(3,1)=axe(1)*axe(3)*(1._vp-COSD(a))-axe(2)*SIND(a)
        R(3,2)=axe(2)*axe(3)*(1._vp-COSD(a))+axe(1)*SIND(a)
        R(3,3)=axe(3)**2._vp+COSD(a)*(1._vp-axe(3)**2._vp)

        vecteur=MATMUL(R,vecteur)

end subroutine rotation_deg

subroutine rotation_rad(axe,vecteur,a)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: axe
        real (kind=vp), dimension(3) :: vecteur
        real (kind=vp) :: a

! local variables
        real (kind=vp), dimension(3,3) :: R

!l'axe doit etre normalise
        axe=axe/norme(axe)

!initialisation de la matrice de rotation autour de l'axe Z

        R(1,1)=axe(1)**2._vp+COS(a)*(1._vp-axe(1)**2._vp)
        R(1,2)=axe(1)*axe(2)*(1._vp-COS(a))-axe(3)*SIN(a)
        R(1,3)=axe(1)*axe(3)*(1._vp-COS(a))+axe(2)*SIN(a)
        R(2,1)=axe(1)*axe(2)*(1._vp-COS(a))+axe(3)*SIN(a)
        R(2,2)=axe(2)**2._vp+COS(a)*(1._vp-axe(2)**2._vp)
        R(2,3)=axe(2)*axe(3)*(1._vp-COS(a))-axe(1)*SIN(a)
        R(3,1)=axe(1)*axe(3)*(1._vp-COS(a))-axe(2)*SIN(a)
        R(3,2)=axe(2)*axe(3)*(1._vp-COS(a))+axe(1)*SIN(a)
        R(3,3)=axe(3)**2._vp+COS(a)*(1._vp-axe(3)**2._vp)

        vecteur=MATMUL(R,vecteur)

end subroutine rotation_rad

!------------------------------------------------------------------------
!-      Outil de calcul: Norme d'un vecteur                             -
!------------------------------------------------------------------------
function norme(vecteur)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) ,intent(in) :: vecteur

        real (kind=vp) :: norme


!Calcul de la norme d'un vecteur

        norme=SQRT(vecteur(1)**2._vp+vecteur(2)**2._vp+vecteur(3)**2._vp)

end function norme


!------------------------------------------------------------------------
!-      Conversion des coordonnees equatoriales d'une source en         -
!-      coordonnees instrument                                          -
!------------------------------------------------------------------------
subroutine conv_equat_to_instr(pointage,source, photon)

implicit none

! IN/OUT variables
        real (kind=vp), dimension(3) :: pointage !en radec
        real (kind=vp), dimension(2) :: source   !en radec

        real (kind=vp), dimension(3) :: photon   !en cartesien

! local variables
        real (kind=vp), dimension(3,3) :: axes_pntg
        real(kind=vp),dimension(3,3) :: matrice
        real(kind=vp),dimension(3) :: tmp
        real(kind=vp),dimension(3) :: vect_src


!coordonnees cartesiennes de chaque axe du repere instrument dans repere equatorial
    !coordonnees des axes au pointage ra=0 et dec=0
    axes_pntg(1,:)=(/ 0._vp ,  0._vp , 1._vp /)
    axes_pntg(2,:)=(/ 0._vp , -1._vp , 0._vp /)
    axes_pntg(3,:)=(/ 1._vp ,  0._vp , 0._vp /)
    !coordonnees des axes au pointage ra et dec
    call rotationZ(axes_pntg(1,:),pointage(1))
    call rotationZ(axes_pntg(2,:),pointage(1))
    call rotationZ(axes_pntg(3,:),pointage(1))
    call rotation_deg(axes_pntg(2,:),axes_pntg(1,:),pointage(2))
    call rotation_deg(axes_pntg(2,:),axes_pntg(2,:),pointage(2))
    call rotation_deg(axes_pntg(2,:),axes_pntg(3,:),pointage(2))


!calcul de la matrice de passage
    matrice(1,:)=axes_pntg(1,:)
    matrice(2,:)=axes_pntg(2,:)
    matrice(3,:)=axes_pntg(3,:)

!coordonnees de la source passee en coordonnees cart
    vect_src=cartesien(source(1),source(2))

!calcul du vecteur direction de la source dans le repere de pointage de l'instrument
    photon=MATMUL(matrice,vect_src)
    photon=-photon

end subroutine conv_equat_to_instr

function cartesien(ra,dec)
    !passage des coordonnees equatoriales (en degres) en cartesien
    real(kind=vp),intent(in) :: ra,dec
    real(kind=vp),dimension(3) :: cartesien

    cartesien(1)=COSD(dec)*COSD(ra)
    cartesien(2)=COSD(dec)*SIND(ra)
    cartesien(3)=SIND(dec)

end function cartesien

end module OUTILS
