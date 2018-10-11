!!
!!
!!

SUBROUTINE tdma(a, b, c, rhs, n)

  INTEGER :: n
  DOUBLE PRECISION, DIMENSION(n) :: a, b, c, rhs
  DOUBLE PRECISION :: m

  INTEGER :: i

  !! Forward elimination
  DO i = 2, n
     m = a(i) / b(i - 1)
     b(i) = b(i) - m * c(i - 1)
     rhs(i) = rhs(i) - m * rhs(i - 1)
  ENDDO

  !! Backward substitution
  rhs(n) = rhs(n) / b(n)
  DO i = n - 1, 1, -1
     rhs(i) = (rhs(i) - c(i) * rhs(i + 1)) / b(i)
  ENDDO
  
ENDSUBROUTINE tdma
