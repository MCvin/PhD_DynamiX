program MIROIRS

    integer, parameter :: wp = kind(1.q0)       !precision des nombres reels du code

!local variables
        character(len=80) :: chaine, mot
        integer :: err, debut, fin

        integer :: nb_miroirs
        real (kind=wp) :: focale
        real (kind=wp), dimension(:), allocatable :: r
        real (kind=wp) :: inc
        integer :: i, n

!parametres de sortie miroirs paraboliques
        real (kind=wp) :: rayon_entree
        real (kind=wp) :: p_para
        real (kind=wp) :: df_p

!parametres de sortie miroirs hyperboliques
        real (kind=wp) :: p_hyper
        real (kind=wp) :: df_h
        real (kind=wp) :: e2

!initialisation des parametres par lecture du fichier parametres.par
    open(1,file='parametres.par')
    err = 0
    do while (err.eq.0)

       read(1,'(A)',iostat=err) chaine

       if ((chaine(1:1).ne.'#') .and. (err.eq.0)) then
          mot=chaine(1:(index(chaine,':'))-1)
          debut=index(chaine,'[')+1
          fin=index(chaine,']')-1

          select case (mot)
            !lecture des parametres de l'observation
             case ('nombres de miroirs')
                read(chaine(debut:fin),*) nb_miroirs
             case ('focale')
                read(chaine(debut:fin),*) focale

          end select

       end if

    end do
    close(1)

!initialisation des rayons des coquilles miroirs
    open(1,file='modeles/input/mirrors.par')
    allocate(r(1:nb_miroirs))
    read(1,*) r !r en mm dans fichier input
    r=r/1000    !r en metres (simu fonctionne avec des metres)
    close(1)


!initialisation des parametres miroirs

        df_h=focale
        !on calcul les inclinaisons des miroirs a partir de la focale du miroir exterieur
        !on se place a la jointure des 2 types de miroirs qui possede une solution
        !sur la parabole et sur l'hyperbole



open(1,file='input/mirrors_para.par')
open(2,file='input/mirrors_hyper.par')
open(3,file='input/mirrors_entree.par')


do i=nb_miroirs,1,-1

   !si on veut etudier une coquille en particulier
   !n=100
   n=i
   !on veut df_h constante donc on calcul inc en fonction de df_h et du rayon de la coquille
   inc=ATAND(r(n)/df_h)/4
   !puis on calcul la df_p correspondante
   df_p=r(n)/TAND(2*inc)

   p_para=r(n)*TAND(inc)

   p_hyper=(r(n)**2-r(n)*df_h*TAND(3*inc))/(df_h+r(n)*TAND(3*inc))

   e2=(TAND(3*inc)*r(n) - p_hyper)/(df_h+p_hyper) +1

   rayon_entree=SQRT(2*(df_p+p_para+0.3)*p_para-p_para**2)


   write(1,*) p_para, df_p+p_para, df_p
   write(2,*) p_hyper, df_h+p_hyper, df_h, e2
   write(3,*) rayon_entree

enddo

close(1)
close(2)
close(3)

end program MIROIRS
