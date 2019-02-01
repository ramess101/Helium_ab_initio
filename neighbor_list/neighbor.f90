!********************************************************************************
MODULE Neighbor
!********************************************************************************

USE Global_Variables
USE Type_Definitions
USE IO_Utilities
USE File_Names
USE Energy_Routines
USE Pair_Nrg_Routines

IMPLICIT NONE


CONTAINS

!******************************************************************************

RECURSIVE SUBROUTINE Rectangle(A,B,is,C,D,js,this_box)
INTEGER, INTENT(IN) :: A,B,C,D, is,js, this_box
INTEGER :: BA, DC, BAm, DCm
INTEGER :: bottom, top, i

BA = B - A
DC = D - C

IF (BA > neighbor_list_granularity(this_box) .AND. &
    DC > neighbor_list_granularity(this_box)) THEN
        BAm = A + BA/2
        DCm = C + DC/2

        !$OMP TASK
       CALL Rectangle(A, BAm, is, C, DCm, js, this_box)
        !$OMP END TASK
        !$OMP TASK
       CALL Rectangle(BAm, B, is, DCm, D,js, this_box)
        !$OMP END TASK
        !$OMP TASKWAIT
        !$OMP TASK
       CALL Rectangle(A, BAm, is, DCm, D, js, this_box)
        !$OMP END TASK
        !$OMP TASK
       CALL Rectangle(BAm, B, is, C, DCm, js, this_box)
        !$OMP END TASK
        !$OMP TASKWAIT

ELSE
        CALL Neighbor_Build_Square(A,B,is,C,D, js,this_box)

END IF

END SUBROUTINE Rectangle

!******************************************************************************

RECURSIVE SUBROUTINE Triangle(A, B, is, this_box)
INTEGER, INTENT(IN) :: is, A,B
INTEGER :: midpoint, length, this_box
length = B-A
IF (length > neighbor_list_granularity(this_box)) THEN

       midpoint = A + length/2
        !$OMP TASK
       CALL Triangle(A, midpoint, is, this_box)
        !$OMP END TASK
        !$OMP TASK
       CALL Triangle(midpoint, B, is, this_box)
        !$OMP END TASK
        !$OMP TASKWAIT
       CALL Rectangle(A, midpoint, is, midpoint, B, is, this_box)
ELSE

       CALL Neighbor_Build_Triangle(A,B, is, this_box)
END IF



END SUBROUTINE Triangle

!******************************************************************************

SUBROUTINE Neighbor_Initialize(this_box)
INTEGER :: A,B,C,D
INTEGER :: this_box
INTEGER :: is, js
 
!$OMP PARALLEL 
!$OMP SINGLE
!$OMP TASK UNTIED
CALL Neighbor_Reset_Count(this_box)

DO is=1, nspecies
    DO js = is, nspecies

       IF (is==js) THEN
          !$OMP TASK 
          CALL Triangle(1, nmols(is,this_box)+1, is, this_box)
          !$OMP END TASK
!       ELSE
!
!          !$OMP TASK 
!          CALL Rectangle(1, nmols(is,this_box), is, &
!                         1, nmols(js,this_box), js, this_box)
!          !$OMP END TASK
!
       END IF
    
    END DO
END DO
!$OMP END TASK
!$OMP END SINGLE
!$OMP END PARALLEL

END SUBROUTINE Neighbor_Initialize

!******************************************************************************

SUBROUTINE Neighbor_Build_Square(A,B,is,C,D,js,this_box)
INTEGER :: this_box, A, B, C, D
INTEGER :: im, is, jm, js, locate_1, locate_2
REAL(DP) :: rxij, ryij, rzij, rxijp, ryijp, rzijp, rijsq

DO im=A,B - 1
   locate_1 = locate(im, is, this_box)
   DO jm = C, D - 1
      locate_2 = locate(jm, js, this_box)
          rxijp = molecule_list(locate_1,is)%xcom - &
                  molecule_list(locate_2,js)%xcom 
          ryijp = molecule_list(locate_1,is)%ycom - &
                  molecule_list(locate_2,js)%ycom 
          rzijp = molecule_list(locate_1,is)%zcom - &
                  molecule_list(locate_2,js)%zcom 

          CALL Minimum_Image_Separation &
               (this_box, rxijp, ryijp, rzijp, rxij, ryij, rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij

          IF (rijsq <= neighbor_rcuthigh_sq(is,this_box) .AND. &
              rijsq >= neighbor_rcutlow_sq(is,this_box)) THEN
                  CALL Neighbor_Append(locate_1,is,locate_2,js,this_box)
          END IF
          IF (rijsq <= neighbor_rcuthigh_sq(js,this_box) .AND. &
              rijsq >= neighbor_rcutlow_sq(js,this_box)) THEN
                  CALL Neighbor_Append(locate_2,js,locate_1,is,this_box)
          END IF

   END DO
END DO


END SUBROUTINE Neighbor_Build_Square

!******************************************************************************

SUBROUTINE Neighbor_Build_Triangle(A,B, is, this_box)
INTEGER :: this_box, A, B, locate_1, locate_2
INTEGER :: im, is, ia, jm, ja
REAL(DP) :: rxij, ryij, rzij, rxijp, ryijp, rzijp, rijsq


DO im=A,B-1

    locate_1 = locate(im, is, this_box)
    
       DO jm = im, B-1

          locate_2 = locate(jm, is, this_box)

             IF (locate_1==locate_2) CYCLE

             rxijp = molecule_list(locate_1,is)%xcom - &
                     molecule_list(locate_2,is)%xcom 
             ryijp = molecule_list(locate_1,is)%ycom - &
                     molecule_list(locate_2,is)%ycom 
             rzijp = molecule_list(locate_1,is)%zcom - &
                     molecule_list(locate_2,is)%zcom 

             CALL Minimum_Image_Separation &
                  (this_box, rxijp, ryijp, rzijp, rxij, ryij, rzij)

             rijsq = rxij*rxij + ryij*ryij + rzij*rzij
             IF (rijsq <= neighbor_rcuthigh_sq(is,this_box) .AND. &
                 rijsq >= neighbor_rcutlow_sq(is,this_box)) THEN
                     CALL Neighbor_Append(locate_1,is, &
                                          locate_2,is,this_box)
                     !IF (locate_1/=locate_2) &
                     CALL Neighbor_Append(locate_2,is, &
                                          locate_1,is,this_box)
             END IF

       END DO


END DO


END SUBROUTINE Neighbor_Build_Triangle

!******************************************************************************

SUBROUTINE Neighbor_Build_Molecule(im, is, this_box)
INTEGER :: this_box, A, B, locate_1, locate_2
INTEGER :: im, is, ia, jm, ja
REAL(DP) :: rxij, ryij, rzij, rxijp, ryijp, rzijp, rijsq

locate_1 = im

molecule_list(locate_1,is)%nbr_neighbors = 0

DO jm = 1,nmols(is,this_box)

   locate_2 = locate(jm, is, this_box)

      IF (locate_1==locate_2) CYCLE

      rxijp = molecule_list(locate_1,is)%xcom - &
              molecule_list(locate_2,is)%xcom 
      ryijp = molecule_list(locate_1,is)%ycom - &
              molecule_list(locate_2,is)%ycom 
      rzijp = molecule_list(locate_1,is)%zcom - &
              molecule_list(locate_2,is)%zcom 

      CALL Minimum_Image_Separation &
           (this_box, rxijp, ryijp, rzijp, rxij, ryij, rzij)

      rijsq = rxij*rxij + ryij*ryij + rzij*rzij
      IF (rijsq <= neighbor_rcuthigh_sq(is,this_box) .AND. &
          rijsq >= neighbor_rcutlow_sq(is,this_box)) THEN
              CALL Neighbor_Append(locate_1,is, &
                                   locate_2,is,this_box)
              !IF (locate_1/=locate_2) &
              CALL Neighbor_Append(locate_2,is, &
                                   locate_1,is,this_box)
      END IF

END DO


END SUBROUTINE Neighbor_Build_Molecule


SUBROUTINE Neighbor_Remove_From_Others(im,is, this_box)
INTEGER, INTENT(IN) ::im, is, this_box
INTEGER :: i,j,k, neighbor_im, neighbor_is

!!$OMP PARALLEL DO DEFAULT(SHARED) &
!!$OMP PRIVATE(i,j,k,neighbor_im,neighbor_is) &
!!$OMP SCHEDULE(STATIC)
DO i=1, molecule_list(im,is)%nbr_neighbors 

   neighbor_im = molecule_list(im,is)%neighbor_list(i)%im
   neighbor_is = molecule_list(im,is)%neighbor_list(i)%is

   DO j=1, molecule_list(neighbor_im, neighbor_is)%nbr_neighbors
      IF (molecule_list(neighbor_im, neighbor_is)%neighbor_list(j)%im==&
          im) THEN
         DO k=j+1,molecule_list(neighbor_im, neighbor_is)%nbr_neighbors
            molecule_list(neighbor_im, neighbor_is)%neighbor_list(k-1)%im = &
                   molecule_list(neighbor_im, neighbor_is)%neighbor_list(k)%im
         END DO
         molecule_list(neighbor_im,neighbor_is)%nbr_neighbors = &
            molecule_list(neighbor_im,neighbor_is)%nbr_neighbors - 1
      END IF
   END DO


END DO
!!$OMP END PARALLEL DO
   
END SUBROUTINE Neighbor_Remove_From_Others

!SUBROUTINE Neighbor_Reset(im,is)
!INTEGER, INTENT(IN) ::im, is
!INTEGER :: i, ia
!
!DO ia = 1, natoms(is)
!   atom_list(ia,im,is)%nbr_neighbors = &
!           atom_list(ia,im,is)%nbr_neighbors_old
!   
!   DO i=1, atom_list(ia,im,is)%nbr_neighbors_old
!   
!           atom_list(ia,im,is)%neighbor_list(i)%ia = &
!                   atom_list(ia,im,is)%neighbor_list_old(i)%ia
!           atom_list(ia,im,is)%neighbor_list(i)%im = &
!                   atom_list(ia,im,is)%neighbor_list_old(i)%im
!           atom_list(ia,im,is)%neighbor_list(i)%is = &
!                   atom_list(ia,im,is)%neighbor_list_old(i)%is
!
!   END DO
!   
!
!END DO
!
!END SUBROUTINE Neighbor_Reset
!
!SUBROUTINE Neighbor_Reset_Box(this_box)
!INTEGER, INTENT(IN) ::this_box
!INTEGER :: i, ia, im, is, locate_1
!
!DO is=1, nspecies
!   DO im=1, nmols(is,this_box)
!      locate_1 = locate(im, is, this_box)
!         DO ia = 1, natoms(is)
!            atom_list(ia,locate_1,is)%nbr_neighbors = &
!                    atom_list(ia,locate_1,is)%nbr_neighbors_old
!            
!            DO i=1, atom_list(ia,im,is)%nbr_neighbors_old
!            
!                    atom_list(ia,locate_1,is)%neighbor_list(i)%ia = &
!                            atom_list(ia,locate_1,is)%neighbor_list_old(i)%ia
!                    atom_list(ia,locate_1,is)%neighbor_list(i)%im = &
!                            atom_list(ia,locate_1,is)%neighbor_list_old(i)%im
!                    atom_list(ia,locate_1,is)%neighbor_list(i)%is = &
!                            atom_list(ia,locate_1,is)%neighbor_list_old(i)%is
!            
!            END DO
!            
!         END DO
!   END DO
!END DO
!
!END SUBROUTINE Neighbor_Reset_Box
!
!SUBROUTINE Neighbor_Update_Others(im,is, this_box)
!INTEGER, INTENT(IN) ::im, is, this_box
!INTEGER :: ia, i, neighbor_ia, neighbor_im, neighbor_is
!
!DO ia = 1, natoms(is)
!
!   DO i=1, atom_list(ia,im,is)%nbr_neighbors 
!   
!           neighbor_ia = atom_list(ia,im,is)%neighbor_list(i)%ia
!           neighbor_im = atom_list(ia,im,is)%neighbor_list(i)%im
!           neighbor_is = atom_list(ia,im,is)%neighbor_list(i)%is
!
!           CALL Neighbor_Append(neighbor_ia,neighbor_im,neighbor_is, &
!                               ia,im,is, this_box) 
!   
!   END DO
!   
!
!END DO
!
!END SUBROUTINE Neighbor_Update_Others
!
!
!
!SUBROUTINE Neighbor_Save(im,is)
!INTEGER, INTENT(IN) ::im, is
!INTEGER :: i, ia
!
!DO ia = 1, natoms(is)
!   atom_list(ia,im,is)%nbr_neighbors_old = &
!           atom_list(ia,im,is)%nbr_neighbors
!   
!   DO i=1, atom_list(ia,im,is)%nbr_neighbors 
!   
!           atom_list(ia,im,is)%neighbor_list_old(i)%ia = &
!                   atom_list(ia,im,is)%neighbor_list(i)%ia
!           atom_list(ia,im,is)%neighbor_list_old(i)%im = &
!                   atom_list(ia,im,is)%neighbor_list(i)%im
!           atom_list(ia,im,is)%neighbor_list_old(i)%is = &
!                   atom_list(ia,im,is)%neighbor_list(i)%is
! 
!   END DO
!   
!   atom_list(ia,im,is)%nbr_neighbors = 0         
!
!END DO
!END SUBROUTINE Neighbor_Save
!
!SUBROUTINE Neighbor_Save_Box(this_box)
!INTEGER, INTENT(IN) ::this_box
!INTEGER :: i, ia, im, is, locate_1
!
!DO is = 1, nspecies
!   DO im = 1 , nmols(is,this_box)
!      locate_1 = locate(im,is,this_box)
!      DO ia = 1, natoms(is)
!      
!          atom_list(ia,locate_1,is)%nbr_neighbors_old = &
!                 atom_list(ia,locate_1,is)%nbr_neighbors
!       
!          DO i=1, atom_list(ia,locate_1, is)%nbr_neighbors
!         
!                 atom_list(ia,locate_1,is)%neighbor_list_old(i)%ia = &
!                         atom_list(ia,locate_1,is)%neighbor_list(i)%ia
!                 atom_list(ia,locate_1,is)%neighbor_list_old(i)%im = &
!                         atom_list(ia,locate_1,is)%neighbor_list(i)%im
!                 atom_list(ia,locate_1,is)%neighbor_list_old(i)%is = &
!                         atom_list(ia,locate_1,is)%neighbor_list(i)%is
!          END DO
!          atom_list(ia,locate_1,is)%nbr_neighbors = 0         
!      END DO
!   END DO
!END DO
!
!END SUBROUTINE Neighbor_Save_Box
!
SUBROUTINE Neighbor_Reset_Count(this_box)
INTEGER, INTENT(IN) ::this_box
INTEGER :: i, ia, im, is, locate_1

DO is = 1, nspecies
   DO im = 1 , nmols(is,this_box)
      locate_1 = locate(im,is,this_box)
      molecule_list(locate_1,is)%nbr_neighbors = 0
   END DO
END DO

END SUBROUTINE Neighbor_Reset_Count


SUBROUTINE Neighbor_Append(im , is, jm, js, this_box)
INTEGER :: im, is, js, jm, i1
INTEGER :: j, k, this_box

molecule_list(im,is)%nbr_neighbors = molecule_list(im,is)%nbr_neighbors + 1

IF (molecule_list(im,is)%nbr_neighbors > neighbor_list_size(is,this_box)) THEN
      err_msg = ""
      err_msg(1) = "Too many neighbors"
      err_msg(2) = "Increase size of neighbor list"
      CALL Clean_Abort(err_msg,'Neighbor_Append')
END IF

i1 = molecule_list(im,is)%nbr_neighbors
molecule_list(im,is)%neighbor_list(i1)%im = jm
molecule_list(im,is)%neighbor_list(i1)%is = js

END SUBROUTINE Neighbor_Append


SUBROUTINE PrintList_(im,is)
INTEGER :: im,is,i
DO i=1, molecule_list(im,is)%nbr_neighbors
        WRITE(*,*) molecule_list(im,is)%neighbor_list(i)%im
END DO
END SUBROUTINE PrintList_


SUBROUTINE Check_Neighbor(im,is,jm,js,is_neighbor)
INTEGER :: im,is,jm,js, i
LOGICAL :: is_neighbor

is_neighbor = .FALSE.
DO i=1, molecule_list(im,is)%nbr_neighbors
        IF (molecule_list(im,is)%neighbor_list(i)%im == jm .AND. &
            molecule_list(im,is)%neighbor_list(i)%is == js) THEN
           is_neighbor = .TRUE.
        END IF
END DO

END SUBROUTINE


END MODULE Neighbor
