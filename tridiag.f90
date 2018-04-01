MODULE TRIDIAGONAL_ROUTINES

	IMPLICIT NONE

	PUBLIC :: Tridiagonal_Solver, Periodic_Tridiagonal_Solver, &
		  Block_Tridiagonal_Solver

	CONTAINS

	SUBROUTINE Tridiagonal_Solver(imin, imax, a, b, c, r, flag)

!		Incoming Variables
		INTEGER,INTENT(IN)			::	imin, imax
		REAL,INTENT(INOUT),DIMENSION(imin:)	::	a, b, c
		REAL,INTENT(INOUT),DIMENSION(imin:)	::	r
		LOGICAL,INTENT(INOUT)				::	flag

!		Local Variables
		INTEGER	::	i
		REAL		::	factor

!		First Executable Statement
		flag = .true.
		DO i = imin,imax-3
			IF (b(i) .EQ. 0.0_8) THEN
!				print*, "ERROR: Zero encountered on main diagonal at ",i
!				STOP
				flag = .false.
				RETURN
			ENDIF
			factor = a(i+1) / b(i)
			b(i+1) = b(i+1) - factor * c(i)
			r(i+1) = r(i+1) - factor * r(i)
		ENDDO

		IF (b(imax-2) .EQ. 0.0_8) THEN
!			print*, "ERROR: Zero encountered on main diagonal at ",i
!			STOP
			flag = .false.
			RETURN
		ENDIF
		factor = a(imax-1) / b(imax-2)
		b(imax-1) = b(imax-1) - factor * c(imax-2)
		r(imax-1) = r(imax-1) - factor * r(imax-2)

		factor = a(imax) / b(imax-1)
		b(imax) = b(imax) - factor * c(imax-1)
		r(imax) = r(imax) - factor * r(imax-1)

		r(imax) = r(imax) / b(imax)
		DO i = imax-1, imin, -1
			r(i) = (r(i) - c(i)*r(i+1)) / b(i)
		ENDDO

		flag = .true.
	END SUBROUTINE Tridiagonal_Solver

	SUBROUTINE Periodic_Tridiagonal_Solver(imin, imax, a, b, c, r, flag)

!		Incoming Variables
		INTEGER,INTENT(IN)			::	imin, imax
		REAL,INTENT(INOUT),DIMENSION(imin:)	::	a, b, c
		REAL,INTENT(INOUT),DIMENSION(imin:)	::	r
		LOGICAL,INTENT(INOUT)				::	flag

!		Local Variables
		INTEGER	::	i
		REAL		::	factor

!		First Executable Statement
		flag = .true.
		DO i = imin,imax-3
			IF (b(i) .EQ. 0.0_8) THEN
!				print*, "ERROR: Zero encountered on main diagonal at ",i
!				STOP
				flag = .false.
				RETURN
			ENDIF
			factor = a(i+1) / b(i)
			b(i+1) = b(i+1) - factor * c(i)
			a(i+1) = - factor * a(i)
			r(i+1) = r(i+1) - factor * r(i)
		ENDDO

		IF (b(imax-2) .EQ. 0.0_8) THEN
!			print*, "ERROR: Zero encountered on main diagonal at ",i
!			STOP
			flag = .false.
			RETURN
		ENDIF
		factor = a(imax-1) / b(imax-2)
		b(imax-1) = b(imax-1) - factor * c(imax-2)
		c(imax-1) = c(imax-1) - factor * a(imax-2)
		r(imax-1) = r(imax-1) - factor * r(imax-2)

		DO i = imin,imax-3
			factor = c(imax) / b(i)
			c(imax) = - factor * c(i)
			b(imax) = b(imax) - factor * a(i)
			r(imax) = r(imax) - factor * r(i)
		ENDDO

		factor = c(imax) / b(imax-2)
		a(imax) = a(imax) - factor * c(imax-2)
		b(imax) = b(imax) - factor * a(imax-2)
		r(imax) = r(imax) - factor * r(imax-2)

		factor = a(imax) / b(imax-1)
		b(imax) = b(imax) - factor * c(imax-1)
		r(imax) = r(imax) - factor * r(imax-1)

		r(imax) = r(imax) / b(imax)
		r(imax-1) = (r(imax-1) - c(imax-1) * r(imax)) / b(imax-1)
		DO i = imax-2, imin, -1
			r(i) = (r(i) - a(i)*r(imax) - c(i)*r(i+1)) / b(i)
		ENDDO
		flag = .true.

	END SUBROUTINE Periodic_Tridiagonal_Solver

	SUBROUTINE Block_Tridiagonal_Solver(imin, imax, nVars, a, b, c, r, flag)

!		Incoming Variables
		INTEGER,INTENT(IN)				:: imin, imax, nVars
		LOGICAL,INTENT(INOUT)					:: flag
		REAL,INTENT(INOUT),DIMENSION(imin:,1:)	:: R
		REAL,INTENT(INOUT),DIMENSION(imin:,1:,1:)	:: A, B, C

!		Local Variables
		INTEGER				:: i, j, k, l
		REAL					:: term
		REAL,DIMENSION(1:nVars,1:nVars)	:: Binv, F

!		First Executable Statement
		DO i = imin, imax-1
			CALL Matrix_Inverse(nVars, B(i,:,:), Binv, flag)
			IF (flag .EQ. .FALSE.)	THEN
				RETURN
			ENDIF
			DO j = 1,nVars
				DO k = 1,nVars
					F(j,k) = 0.0_8
					DO l = 1,nVars
						F(j,k) = F(j,k) + A(i+1,j,l) * Binv(l,k)
					ENDDO
				ENDDO
			ENDDO
			DO j = 1,nVars
				DO k = 1,nVars
					term = 0.0_8
					DO l = 1, nVars
						term = term + F(j,l) * C(i,l,k)
					ENDDO
					B(i+1,j,k) = B(i+1,j,k) - term
				ENDDO
			ENDDO
			DO j = 1,nVars
				term = 0.0_8
				DO l = 1,nVars
					term = term + F(j,l) * R(i,l)
				ENDDO
				R(i+1,j) = R(i+1,j) - term
			ENDDO
		ENDDO
		CALL LU_Decomposition(nVars, B(imax,:,:), R(imax,:), flag)
		IF (flag .EQ. .FALSE.)	THEN
			RETURN
		ENDIF
		DO i = imax-1,imin,-1
			DO j = 1,nVars
				term = 0.0_8
				DO l = 1,nVars
					term = term + C(i,j,l) * R(i+1,l)
				ENDDO
				R(i,j) = R(i,j) - term
			ENDDO
			CALL LU_Decomposition(nVars, B(i,:,:), R(i,:), flag)
			IF (flag .EQ. .FALSE.)	THEN
				RETURN
			ENDIF
		ENDDO
		flag = .true.

	END SUBROUTINE Block_Tridiagonal_Solver

	SUBROUTINE Matrix_Inverse(n, LHS, RHS, flag)

!		Incoming Variables
		INTEGER,INTENT(IN)			:: n
		REAL,INTENT(IN),DIMENSION(1:,1:)	:: LHS
		REAL,INTENT(INOUT),DIMENSION(1:,1:)	:: RHS
		LOGICAL,INTENT(INOUT)				:: flag

!		Local variables
		INTEGER				:: row,col,dindex
		REAL,DIMENSION(1:n,1:2*n)		:: arr2
		REAL					:: tempval, wval

!		First Executable statement
		DO row = 1,n
			DO col = 1,n
				arr2(row,col) = LHS(row,col)
			ENDDO
			DO col = n+1,2*n
				IF (col-n .EQ. row) THEN
					arr2(row,col) = 1.0_8
				ELSE
					arr2(row,col) = 0.0_8
				ENDIF
			ENDDO
		ENDDO

		DO dindex = 1,n
			IF ((dindex .EQ. n) .AND. (arr2(dindex,dindex) .EQ. 0.0_8)) THEN
				flag = .FALSE.
				RETURN
			ELSEIF (arr2(dindex,dindex) .EQ. 0.0_8) THEN
				CALL swaprows(2*n, arr2, dindex, dindex+1)
			ENDIF
			IF (arr2(dindex,dindex) .EQ. 0.0_8) THEN
				flag = .FALSE.
				RETURN
			ELSE
				tempval = arr2(dindex,dindex)
				DO col = 1,2*n
					arr2(dindex,col) = arr2(dindex,col) / tempval
				ENDDO
			ENDIF
			DO row = dindex+1,n
				wval = arr2(row,dindex)
				DO col = 1,2*n
					arr2(row,col) = arr2(row,col) - wval*arr2(dindex,col)
				ENDDO
			ENDDO
		ENDDO

		DO dindex = n,1,-1
			DO row = dindex-1,1,-1
				wval = arr2(row,dindex)
				DO col = 1,2*n
					arr2(row,col) = arr2(row,col) - wval*arr2(dindex,col)
				ENDDO
			ENDDO
		ENDDO

		DO row = 1,n
			DO col = 1,n
				RHS(row,col) = arr2(row,col+n)
			ENDDO
		ENDDO
		flag = .TRUE.

	END SUBROUTINE Matrix_Inverse

	SUBROUTINE SWAPROWS(n, A, i, j)
!		Incoming Variables
		REAL,INTENT(INOUT),DIMENSION(1:,1:)	:: A
		INTEGER,INTENT(IN)			:: n, i, j

!		Local Variables
		INTEGER		:: m
		REAL		:: temp

!		First Executable statement
		DO m = 1,n
			temp = A(i,m)
			A(i,m) = A(j,m)
			A(j,m) = temp
		ENDDO

	END SUBROUTINE SWAPROWS

	SUBROUTINE LU_Decomposition(n, A, R, flag)

!		Incoming Variables
		INTEGER,INTENT(IN)			:: n
		REAL,INTENT(INOUT),DIMENSION(1:,1:)	:: A
		REAL,INTENT(INOUT),DIMENSION(1:)	:: R
		LOGICAL,INTENT(INOUT)			:: flag

!		Local Variables
		INTEGER		:: i, j, k
		REAL		:: factor, term

!		First Executable Statement
		flag = .true.
		DO i = 1, n
			IF (A(i,i) .EQ. 0.0_8) THEN
				flag = .FALSE.
				RETURN
			ENDIF
			DO j = i+1, n
				factor = A(j,i) / A(i,i)
				A(j,i) = 0.0_8
				DO k = i+1, n
					A(j,k) = A(j,k) - factor * A(i,k)
				ENDDO
				R(j) = R(j) - factor * R(i)
			ENDDO
		ENDDO
		R(n) = R(n) / A(n,n)
		DO i = n-1,1,-1
			term = 0.0_8
			DO j = i+1,n
				term = term + A(i,j) * R(j)
			ENDDO
			R(i) = (R(i) - term) / A(i,i)
		ENDDO
		flag = .true.

	END SUBROUTINE LU_Decomposition

END MODULE TRIDIAGONAL_ROUTINES
