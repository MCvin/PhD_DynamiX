program MIROIRS

!local variables
        character(len=80) :: chaine, mot
        integer :: err, debut, fin
        integer :: nb_miroirs
        real (kind=16) :: focale
        real (kind=16) :: r, rayon_ext
        real (kind=16) :: inc_ext, inc_int, inc, inc_increment
        integer :: i

!parametres de sortie miroirs paraboliques
        real (kind=16) :: rayon_entree
        real (kind=16) :: p_para
        real (kind=16) :: df_p

!parametres de sortie miroirs hyperboliques
        real (kind=16) :: p_hyper
        real (kind=16) :: df_h
        real (kind=16) :: e2

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
             case ('nombres de miroirs')
                read(chaine(debut:fin),*) nb_miroirs
             case ('focale')
                read(chaine(debut:fin),*) focale
             case ('rayon du miroir exterieur')
                read(chaine(debut:fin),*) rayon_ext
          end select

       end if

    end do
    close(1)

!initialisation des parametres miroirs

        df_h=focale
        !on calcul les inclinaisons des miroirs a partir de la focale du miroir exterieur
        !on se place a la jointure des 2 types de miroirs qui possede une solution
        !sur la parabole et sur l'hyperbole
        inc_ext=(ATAND(rayon_ext/focale))/4
        !a verifier: verifie avec les italiens ok
        inc_int=inc_ext/2


open(1,file='../input/mirrors_para_old.par')
open(2,file='../input/mirrors_hyper_old.par')
open(3,file='../input/mirrors_entree_old.par')

inc_increment=(inc_ext-inc_int)/(nb_miroirs-1)

do i=1,nb_miroirs

   inc=inc_int+(i-1)*inc_increment
   
   !on veut df_h constante donc on calcul r en fonction de df_h et de l'angle d'attaque
   r=df_h*TAND(4*inc)
   !puis on calcul la df_p correspondante
   df_p=r/TAND(2*inc)

   p_para=r*TAND(inc)

   p_hyper=(r**2-r*df_h*TAND(3*inc))/(df_h+r*TAND(3*inc))

   e2=(TAND(3*inc)*r - p_hyper)/(df_h+p_hyper) +1

   rayon_entree=SQRT(2*(df_p+p_para+0.3)*p_para-p_para**2)

   write(1,*) p_para, df_p+p_para, df_p
   write(2,*) p_hyper, df_h+p_hyper, df_h, e2
   write(3,*) rayon_entree

enddo

close(1)
close(2)
close(3)

end program MIROIRS
