#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

      subroutine filter_margin( var_adj, my_value, myThid )

C     this subroutine is created by Andrey.Vlasenko and used to remove spike sensitivity
C     in the sensitivity variables in the adjoint model
C     it is used in the study of "A Pilot Climate Sensitivity Study Using the CEN Coupled Adjoint Model (CESAM)"
C     "Seasonal Prediction of Northern European Winter Air Temperatures From SST Anomalies Based on Sensitivity Estimates",
C     and "Arctic Ocean-Sea ice reanalysis for the period 2007-2016 using the adjoint method"

C     var_adj - is a variable used as input(for filtering) and output(after filtering)
C     my_value - is the threshhold value for filtering.

      IMPLICIT NONE
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

    

      real*8 var_adj (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
      real*8 my_value

      real*8 a1
      real*8 a2
      real*8 a3
      real*8 a4
      real*8 a5
      real*8 a6
      real*8 a7
      real*8 a8
      real*8 aux (1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)  
      real*8 b1
      real*8 b2
      real*8 b3
      real*8 b4
      real*8 b5
      real*8 b6
      real*8 b7
      real*8 b8
      real*8 b9
      real*8 dev1,dev2,dev3
      INTEGER i, j, bi, bj
      INTEGER myThid
      real*8 ma
      real*8 mi
      real*8 mm
      real*8 s
      real*8 s1,s2,s3,s01,s02


! averaged over the subtile
        dev1=0. 
        ma=0.
        do bj = myByLo(myThid),myByHi(myThid)
            do bi = myBxLo(myThid),myBxHi(myThid)
              do j = 1-oly, sNy+oly
                do i = 1-olx, sNx+olx
                  if(abs(var_adj(i,j,bi,bj)) .gt. 10.e-20)then
                    dev1=dev1+1.;
                    ma=ma+abs(var_adj(i,j,bi,bj))
                  endif                   
                enddo
               enddo
             enddo
         enddo

         CALL GLOBAL_sum_R8 (ma,  myThid  )
         CALL GLOBAL_sum_R8 (dev1,  myThid  )
C do the average
C        IF(dev1)
         ma=ma/(dev1+1.);      
        do bj = myByLo(myThid),myByHi(myThid)
            do bi = myBxLo(myThid),myBxHi(myThid)
              do j = 1, sNy
                do i = 1, sNx
                  if(abs(var_adj(i,j,bi,bj)) .gt. 30.*ma)then
            var_adj(i,j,bi,bj)= 0. !sign(30.*ma,var_adj(i,j,bi,bj))
                  endif                   
                enddo
               enddo
             enddo
         enddo
             aux = var_adj
      do bj = myByLo(myThid),myByHi(myThid)
            do bi = myBxLo(myThid),myBxHi(myThid)
              do j = 1-oly, sNy+oly
                do i = 1-olx, sNx+olx
                   s = 0
                   s1 = 0
                   dev1 = 0
                   dev2 = 0
                   s2 = 0.
                   dev3=0
                   s3=0
                   a1 = var_adj(i+1,j,bi,bj)
                   a2 = var_adj(i-1,j,bi,bj)
                   a3 = var_adj(i,j+1,bi,bj)
                   a4 = var_adj(i,j-1,bi,bj)
                   a5 = var_adj(i+1,j+1,bi,bj)
                   a6 = var_adj(i-1,j-1,bi,bj)
                   a7 = var_adj(i+1,j-1,bi,bj)
                   a8 = var_adj(i-1,j+1,bi,bj)
                   call nonzero_count1( a1,b1,dev1 )
                   call nonzero_count1( a2,b2,dev1 )
                   call nonzero_count1( a3,b3,dev1 )
                   call nonzero_count1( a4,b4,dev1 )
                   call nonzero_count1( a5,b5,dev2 )
                   call nonzero_count1( a6,b6,dev2 )
                   call nonzero_count1( a7,b7,dev2 )
                   call nonzero_count1( a8,b8,dev2 )

                   b9 = abs(aux(i,j,bi,bj))
                   if(dev1.eq.0) dev1=1.
                   if(dev2.eq.0) dev2=1.
                   if (dev1 .ge. 1 ) then
                      s01 = (a1+a2+a3+a4)/dev1
                      s02 = (a5+a6+a7+a8)/dev2
                      s=(a1+a2+a3+a4+a5+a6+a7+a8)/(dev1+dev2)

                 s1 = my_value*(b1+b2+b3+b4)/dev1
                 s2 = my_value*(b5+b6+b7+b8)/dev2
                 s3 = my_value*(b1+b2+b3+b4+b5+b6+b7+b8)/(dev2+dev1)
                      if (b9 .gt. s3 ) then
                          aux(i,j,bi,bj) = 0.
                      else  ! b9.gt. s3
                   if(b9 .gt. s1 
     $.and. b9.lt. s2 )then
                         aux(i,j,bi,bj) = 0.
                   elseif(b9 .lt. s1 
     $.and. b9.gt. s2 )then
                      aux(i,j,bi,bj) = 0.                      
                      endif
                   endif
                 endif
                end do
             end do

          end do
       end do
!       var_adj(1:sNx,1:sNy,:,:) = aux(1:sNx,1:sNy,:,:)
      end subroutine



      subroutine nonzero_count1( a, b, dev )

      implicit none


      double precision, intent(inout) :: a
      double precision, intent(out) :: b
      double precision, intent(out) :: dev

      if(a.ne.a) a = 0.

      if (a .lt. 0) then
        b = -(1.*a)
      else
        b = a
      endif
      if (b .gt. 10e-30) then
        dev = dev+1.
      endif
      end subroutine nonzero_count1
