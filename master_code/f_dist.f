C FILE: F_DIST.F
      SUBROUTINE f_dist(n1,n2,x1,y1,z1,x2,y2,z2,result)
C
C     CALCULATRE DISTANCES USING FORTRAN
C
      DIMENSION x1(*), y1(*), z1(*), x2(*), y2(*), z2(*), result(*)
      INTEGER n1,n2,k
      REAL x1,y1,z1,x2,y2,z2,result,dx,dy,dz
Cf2py intent(in, out, overwrite) result
C
      k = 0
      DO j = 1,n2
         DO i = 1,n1
            k = k + 1
            dx = x2(j) - x1(i)
            dy = y2(j) - y1(i)
            dz = z2(j) - z1(i)
            result(k) = dx*dx + dy*dy + dz*dz
         ENDDO
      ENDDO
      END
C END FILE F_DIST.F

