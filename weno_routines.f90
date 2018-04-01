MODULE WENO_ROUTINES

	USE PARAMS_GLOBAL
	USE TRIDIAGONAL_ROUTINES
	IMPLICIT NONE

	PUBLIC :: WENO5_SCALAR, CRWENO5_SCALAR

	CONTAINS

	SUBROUTINE WENO5_SCALAR(Q_Cell, Q_Face_Left, Q_Face_Right, imin, imax, int_junk1, real_junk1, real_junk2, real_junk3, fmin, fmax, ibmin, ibmax, mdim)

!		Incoming Variables
		INTEGER,INTENT(IN) :: imin,imax, mdim
		REAL,INTENT(INOUT),DIMENSION(imin:,1:) :: Q_Cell
		REAL,INTENT(INOUT),DIMENSION(imin:,1:) :: Q_Face_Left,Q_Face_Right
		INTEGER,INTENT(INOUT) :: int_junk1, ibmin, ibmax
		REAL,INTENT(INOUT),DIMENSION(1:) :: fmin, fmax
		REAL,INTENT(INOUT) ::  real_junk1, real_junk2, real_junk3

!		Local Variables
		INTEGER,PARAMETER ::	nmv=5
		INTEGER		::	i, n
		REAL		::	blank_fac, eps, tau
		REAL		::	w1_opt, w2_opt, w3_opt
		REAL		::	f1, f2, f3, b1, b2, b3
		REAL		::	a1, a2, a3, w1, w2, w3
		REAL 		:: 	rho, u, v, E, p
		REAL		::	Q_m2, Q_m1, Q_c, Q_p1, Q_p2
		LOGICAL		:: 	Borges, Carpenter, No_Limiting, Mapped

!		First executable statement
		eps = 1e-6
		Mapped = .FALSE.
		Borges = .FALSE.
		Carpenter = .FALSE.
		No_Limiting = .FALSE.

		DO n = 1,nmv
!			Boundary treatment
!			No need to reconstruct Q_Face_Right(imin) and
!			Q_Face_Left(imax)

!			Left Boundary
			Q_Face_Left(imin,n) =  (1.0_8/3.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) - (1.0_8/6.0_8)*Q_Cell(imin+2,n)
			Q_Face_Left(imin+1,n) = -(1.0_8/6.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) + (1.0_8/3.0_8)*Q_Cell(imin+2,n)
			Q_Face_Right(imin+1,n) = (1.0_8/3.0_8)*Q_Cell(imin+3,n) - (7.0_8/6.0_8)*Q_Cell(imin+2,n) + (11.0_8/6.0_8)*Q_Cell(imin+1,n)
!			Right Boundary
			Q_Face_Right(imax,n) = (1.0_8/3.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) - (1.0_8/6.0_8)*Q_Cell(imax-2,n)
			Q_Face_Right(imax-1,n) = -(1.0_8/6.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) + (1.0_8/3.0_8)*Q_Cell(imax-2,n)
			Q_Face_Left(imax-1,n) = (1.0_8/3.0_8)*Q_Cell(imax-3,n) - (7.0_8/6.0_8)*Q_Cell(imax-2,n) + (11.0_8/6.0_8)*Q_Cell(imax-1,n)

!			Interior
			DO i = imin+2, imax-2
!				Defining stencil points
				Q_m2 = Q_Cell(i-2, n)
				Q_m1 = Q_Cell(i-1, n)
				Q_c  = Q_Cell(i, n)
				Q_p1 = Q_Cell(i+1, n)
				Q_p2 = Q_Cell(i+2, n)
!				Candidate Stencils
				f1 =  (1.0_8/3.0_8)*Q_c + (5.0_8/6.0_8)*Q_p1 - (1.0_8/6.0_8)*Q_p2
				f2 = -(1.0_8/6.0_8)*Q_m1 + (5.0_8/6.0_8)*Q_c + (1.0_8/3.0_8)*Q_p1
				f3 =  (1.0_8/3.0_8)*Q_m2 - (7.0_8/6.0_8)*Q_m1 + (11.0_8/6.0_8)*Q_c
!				Smoothness Indicators
				b1 = (13.0_8/12.0_8) * (Q_c-2.0_8*Q_p1+Q_p2)*(Q_c-2.0_8*Q_p1+Q_p2) &
					+ (1.0_8/4.0_8) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2)*(3.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				b2 = (13.0_8/12.0_8) * (Q_m1-2.0_8*Q_c+Q_p1)*(Q_m1-2.0_8*Q_c+Q_p1) &
					+ (1.0_8/4.0_8) * (Q_m1-Q_p1)*(Q_m1-Q_p1)
				b3 = (13.0_8/12.0_8) * (Q_m2-2.0_8*Q_m1+Q_c)*(Q_m2-2.0_8*Q_m1+Q_c) &
					+ (1.0_8/4.0_8) * (Q_m2-4.0_8*Q_m1+3.0_8*Q_c)*(Q_m2-4.0_8*Q_m1+3.0_8*Q_c)
				IF (Borges) THEN
!					Borges definition of tau
					tau = abs (b1 - b3)
				ELSEIF (Carpenter) THEN
!					Yamaleev-Carpenter definition of tau
					tau = (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2) &
					    * (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				ELSE
					tau = 0.0_8
				ENDIF
!				Optimal Weights
				w1_opt = 0.3
				w2_opt = 0.6
				w3_opt = 0.1
!				WENO limiting
				IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
					a1 = w1_opt * (1.0_8 + tau / (eps+b1))
					a2 = w2_opt * (1.0_8 + tau / (eps+b2))
					a3 = w3_opt * (1.0_8 + tau / (eps+b3))
				ELSE
					a1 = w1_opt / ((eps+b1)*(eps+b1))
					a2 = w2_opt / ((eps+b2)*(eps+b2))
					a3 = w3_opt / ((eps+b3)*(eps+b3))
				ENDIF
!				Making the weights convex
				w1 = a1 / (a1 + a2 + a3)
				w2 = a2 / (a1 + a2 + a3)
				w3 = a3 / (a1 + a2 + a3)
				IF (Mapped) THEN
!					Mapping the WENO weights
					a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
					a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
					a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
				ENDIF

				Q_Face_Left(i,n) = w1*f1 + w2*f2 + w3*f3

!				Defining stencil points
				Q_m2 = Q_Cell(i+2, n)
				Q_m1 = Q_Cell(i+1, n)
				Q_c  = Q_Cell(i, n)
				Q_p1 = Q_Cell(i-1, n)
				Q_p2 = Q_Cell(i-2, n)
!				Candidate Stencils
				f1 =  (1.0_8/3.0_8)*Q_c + (5.0_8/6.0_8)*Q_p1 - (1.0_8/6.0_8)*Q_p2
				f2 = -(1.0_8/6.0_8)*Q_m1 + (5.0_8/6.0_8)*Q_c + (1.0_8/3.0_8)*Q_p1
				f3 =  (1.0_8/3.0_8)*Q_m2 - (7.0_8/6.0_8)*Q_m1 + (11.0_8/6.0_8)*Q_c
!				Smoothness indicators
				b1 = (13.0_8/12.0_8) * (Q_c-2.0_8*Q_p1+Q_p2)*(Q_c-2.0_8*Q_p1+Q_p2) &
					+ (1.0_8/4.0_8) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2)*(3.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				b2 = (13.0_8/12.0_8) * (Q_m1-2.0_8*Q_c+Q_p1)*(Q_m1-2.0_8*Q_c+Q_p1) &
					+ (1.0_8/4.0_8) * (Q_m1-Q_p1)*(Q_m1-Q_p1)
				b3 = (13.0_8/12.0_8) * (Q_m2-2.0_8*Q_m1+Q_c)*(Q_m2-2.0_8*Q_m1+Q_c) &
					+ (1.0_8/4.0_8) * (Q_m2-4.0_8*Q_m1+3.0_8*Q_c)*(Q_m2-4.0_8*Q_m1+3.0_8*Q_c)
				IF (Borges) THEN
!					Borges definition of tau
					tau = abs (b1 - b3)
				ELSEIF (Carpenter) THEN
!					Yamaleev-Carpenter definition of tau
					tau = (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2) &
					    * (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				ELSE
					tau = 0.0_8
				ENDIF
!				Optimal Weights
				w1_opt = 0.3
				w2_opt = 0.6
				w3_opt = 0.1
!				WENO limiting
				IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
					a1 = w1_opt * (1.0_8 + tau / (eps+b1))
					a2 = w2_opt * (1.0_8 + tau / (eps+b2))
					a3 = w3_opt * (1.0_8 + tau / (eps+b3))
				ELSE
					a1 = w1_opt / ((eps+b1)*(eps+b1))
					a2 = w2_opt / ((eps+b2)*(eps+b2))
					a3 = w3_opt / ((eps+b3)*(eps+b3))
				ENDIF
!				Making the weights convex
				w1 = a1 / (a1 + a2 + a3)
				w2 = a2 / (a1 + a2 + a3)
				w3 = a3 / (a1 + a2 + a3)
				IF (Mapped) THEN
!					Mapping the WENO weights
					a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
					a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
					a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
				ENDIF

				Q_Face_Right(i,n) = w1*f1 + w2*f2 + w3*f3
			ENDDO !i loop

!			Overset correction
!			DO i = imin+2, imax-2
!				blank_fac  = iba(i)*iba(i-1)*iba(i+1)*iba(i-2)*iba(i+2)
!				Q_Face_Left(i,n) = Q_Face_Left(i,n)*blank_fac+(1.-blank_fac)*Q_Cell(i,n)
!				Q_Face_Right(i,n) = Q_Face_Right(i,n)*blank_fac+(1.-blank_fac)*Q_Cell(i,n)
!			ENDDO !i loop
		ENDDO	!n loop

	END SUBROUTINE WENO5_SCALAR

	SUBROUTINE CRWENO5_SCALAR(Q_Cell, Q_Face_Left, Q_Face_Right, imin, imax, int_junk1, real_junk1, real_junk2, real_junk3, fmin, fmax, ibmin, ibmax, mdim)

!		Incoming Variables
		INTEGER,INTENT(IN) :: imin,imax, mdim
		REAL,INTENT(INOUT),DIMENSION(imin:,1:) :: Q_Cell
		REAL,INTENT(INOUT),DIMENSION(imin:,1:) :: Q_Face_Left,Q_Face_Right
		INTEGER,INTENT(INOUT) :: int_junk1, ibmin, ibmax
		REAL,INTENT(INOUT),DIMENSION(1:) :: fmin, fmax
		REAL,INTENT(INOUT) ::  real_junk1, real_junk2, real_junk3

!		Local Variables
		INTEGER,PARAMETER ::	nmv=5
		INTEGER		::	i, n
		REAL		::	blank_fac, eps, tau
		REAL		::	w1_opt, w2_opt, w3_opt
		REAL		::	f1, f2, f3, b1, b2, b3
		REAL		::	a1, a2, a3, w1, w2, w3
		REAL 		:: 	rho, u, v, E, p
		REAL		::	Q_m2, Q_m1, Q_c, Q_p1, Q_p2, Q_p3
		LOGICAL		:: 	Borges, Carpenter, No_Limiting, Mapped, Variable_Eps
		REAL		::	vareps_min, vareps_max, term, IS_p
		REAL,ALLOCATABLE,DIMENSION(:) :: Q_Copy
		REAL,ALLOCATABLE	:: 	a(:), b(:), c(:), r(:)
		LOGICAL				::	tridiag_flag

!		First executable statement
		eps = 1e-6
		Mapped = .FALSE.
		Borges = .FALSE.
		Carpenter = .FALSE.
		No_Limiting = .FALSE.
		IS_p = 2.0

		Variable_Eps = .FALSE.
		vareps_min = 1e-20
		vareps_max = 1e-6

		DO n = 1,nmv
			ALLOCATE(a(imin:imax))
			ALLOCATE(b(imin:imax))
			ALLOCATE(c(imin:imax))
			ALLOCATE(r(imin:imax))

!			Left Biased Interpolation

!			Left Boundary (Non-compact)
			a(imin) = 0.0
			b(imin) = 1.0
			c(imin) = 0.0
			r(imin) =  (1.0_8/3.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) - (1.0_8/6.0_8)*Q_Cell(imin+2,n)
			a(imin+1) = 0.0
			b(imin+1) = 1.0
			c(imin+1) = 0.0
			r(imin+1) = -(1.0_8/6.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) + (1.0_8/3.0_8)*Q_Cell(imin+2,n)
!			Right Boundary (Non-compact)
			a(imax-1) = 0.0
			b(imax-1) = 1.0
			c(imax-1) = 0.0
			r(imax-1) = (1.0_8/3.0_8)*Q_Cell(imax-3,n) - (7.0_8/6.0_8)*Q_Cell(imax-2,n) + (11.0_8/6.0_8)*Q_Cell(imax-1,n)

!			Left Boundary (Compact)
!			a(imin) = 0.0
!			b(imin) = 2.0/3.0
!			c(imin) = 1.0/3.0
!			r(imin) = (1.0/6.0) * (Q_Cell(imin,n) + 5.0*Q_Cell(imin+1,n))
!			a(imin+1) = 3.0/10.0
!			b(imin+1) = 6.0/10.0
!			c(imin+1) = 1.0/10.0
!			r(imin+1) = (1.0/30.0) * (Q_Cell(imin,n) + 19.0*Q_Cell(imin+1,n) + 10.0*Q_Cell(imin+2,n))
!			Right Boundary (Compact)
!			a(imax-1) = 1.0/3.0
!			b(imax-1) = 2.0/3.0
!			c(imax-1) = 0.0
!			r(imax-1) = (1.0/6.0) * (5.0*Q_Cell(imax-1,n) + Q_Cell(imax,n))
!			Redundant
			a(imax) = 0.0
			b(imax) = 1.0
			c(imax) = 0.0
			r(imax) = 0.0
!			Interior
			DO i = imin+2, imax-2
!				Defining stencil points
				Q_m2 = Q_Cell(i-2, n)
				Q_m1 = Q_Cell(i-1, n)
				Q_c  = Q_Cell(i, n)
				Q_p1 = Q_Cell(i+1, n)
				Q_p2 = Q_Cell(i+2, n)
!				Candidate stencils
				f1 = (Q_m1 + 5.0_8*Q_c) / 6.0_8
				f2 = (5.0_8*Q_c + Q_p1) / 6.0_8
				f3 = (Q_c + 5.0_8*Q_p1) / 6.0_8
!				Smoothness Indicators
				b1 = (13.0_8/12.0_8) * (Q_m2-2.0_8*Q_m1+Q_c)*(Q_m2-2.0_8*Q_m1+Q_c) &
					+ (1.0_8/4.0_8) * (Q_m2-4.0_8*Q_m1+3.0_8*Q_c)*(Q_m2-4.0_8*Q_m1+3.0_8*Q_c)
				b2 = (13.0_8/12.0_8) * (Q_m1-2.0_8*Q_c+Q_p1) * (Q_m1-2.0_8*Q_c+Q_p1) &
					+ (1.0_8/4.0_8) * (Q_m1-Q_p1) * (Q_m1-Q_p1)
				b3 = (13.0_8/12.0_8) * (Q_c-2.0_8*Q_p1+Q_p2) * (Q_c-2.0_8*Q_p1+Q_p2) &
					+ (1.0_8/4.0_8) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				
				IF (Variable_Eps) THEN
					term = MIN(b1,b2,b3) / (vareps_min + MAX(b1,b2,b3) - MIN(b1,b2,b3))
					eps = vareps_max * MIN(1.0_8,term) + vareps_min
				ENDIF

				IF (Borges) THEN
!					Borges definition of tau
					tau = abs (b1 - b3)
				ELSEIF (Carpenter) THEN
!					Yamaleev-Carpenter definition of tau
					tau = (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2) &
					    * (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2)
				ELSE
					tau = 0.0_8
				ENDIF
!				Optimal Weights
				w1_opt = 0.2_8
				w2_opt = 0.5_8
				w3_opt = 0.3_8
!				WENO Limiting
				IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
					a1 = w1_opt * (1.0_8 + (tau/(eps+b1))**IS_p)
					a2 = w2_opt * (1.0_8 + (tau/(eps+b2))**IS_p)
					a3 = w3_opt * (1.0_8 + (tau/(eps+b3))**IS_p)
				ELSE
					a1 = w1_opt / ((eps+b1)**IS_p)
					a2 = w2_opt / ((eps+b2)**IS_p)
					a3 = w3_opt / ((eps+b3)**IS_p)
				ENDIF
!				Making the weights convex
				w1 = a1 / (a1 + a2 + a3)
				w2 = a2 / (a1 + a2 + a3)
				w3 = a3 / (a1 + a2 + a3)
				IF (Mapped) THEN
!					Mapping the WENO weights
					a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
					a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
					a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
				ENDIF

				r(i) = w1*f1 + w2*f2 + w3*f3
				a(i) = (2.0_8/3.0_8) * w1 + (1.0_8/3.0_8) * w2
				b(i) = (1.0_8/3.0_8) * w1 + (2.0_8/3.0_8) * (w2 + w3)
				c(i) = (1.0_8/3.0_8) * w3
			ENDDO !i-loop
			CALL TRIDIAGONAL_SOLVER(imin, imax, a, b, c, r, tridiag_flag)
			IF (tridiag_flag) THEN
				DO i = imin, imax-1
					Q_Face_Left(i,n) = r(i)
				ENDDO !i loop
			ELSE
!				Falling back on WENO5 if singular system for the
!				compact scheme
				Q_Face_Left(imin,n) =  (1.0_8/3.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) - (1.0_8/6.0_8)*Q_Cell(imin+2,n)
				Q_Face_Left(imin+1,n) = -(1.0_8/6.0_8)*Q_Cell(imin,n) + (5.0_8/6.0_8)*Q_Cell(imin+1,n) + (1.0_8/3.0_8)*Q_Cell(imin+2,n)
				Q_Face_Left(imax-1,n) = (1.0_8/3.0_8)*Q_Cell(imax-3,n) - (7.0_8/6.0_8)*Q_Cell(imax-2,n) + (11.0_8/6.0_8)*Q_Cell(imax-1,n)
!				Interior
				DO i = imin+2, imax-2
!					Defining stencil points
					Q_m2 = Q_Cell(i-2, n)
					Q_m1 = Q_Cell(i-1, n)
					Q_c  = Q_Cell(i, n)
					Q_p1 = Q_Cell(i+1, n)
					Q_p2 = Q_Cell(i+2, n)
!					Candidate Stencils
					f1 =  (1.0_8/3.0_8)*Q_c + (5.0_8/6.0_8)*Q_p1 - (1.0_8/6.0_8)*Q_p2
					f2 = -(1.0_8/6.0_8)*Q_m1 + (5.0_8/6.0_8)*Q_c + (1.0_8/3.0_8)*Q_p1
					f3 =  (1.0_8/3.0_8)*Q_m2 - (7.0_8/6.0_8)*Q_m1 + (11.0_8/6.0_8)*Q_c
!					Smoothness Indicators
					b1 = (13.0_8/12.0_8) * (Q_c-2.0_8*Q_p1+Q_p2)*(Q_c-2.0_8*Q_p1+Q_p2) &
						+ (1.0_8/4.0_8) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2)*(3.0_8*Q_c-4.0_8*Q_p1+Q_p2)
					b2 = (13.0_8/12.0_8) * (Q_m1-2.0_8*Q_c+Q_p1)*(Q_m1-2.0_8*Q_c+Q_p1) &
						+ (1.0_8/4.0_8) * (Q_m1-Q_p1)*(Q_m1-Q_p1)
					b3 = (13.0_8/12.0_8) * (Q_m2-2.0_8*Q_m1+Q_c)*(Q_m2-2.0_8*Q_m1+Q_c) &
						+ (1.0_8/4.0_8) * (Q_m2-4.0_8*Q_m1+3.0_8*Q_c)*(Q_m2-4.0_8*Q_m1+3.0_8*Q_c)
					IF (Borges) THEN
!						Borges definition of tau
						tau = abs (b1 - b3)
					ELSEIF (Carpenter) THEN
!						Yamaleev-Carpenter definition of tau
						tau = (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2) &
						    * (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2)
					ELSE
						tau = 0.0_8
					ENDIF
!					Optimal Weights
					w1_opt = 0.3
					w2_opt = 0.6
					w3_opt = 0.1
!					WENO limiting
					IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
						a1 = w1_opt * (1.0_8 + tau / (eps+b1))
						a2 = w2_opt * (1.0_8 + tau / (eps+b2))
						a3 = w3_opt * (1.0_8 + tau / (eps+b3))
					ELSE
						a1 = w1_opt / ((eps+b1)*(eps+b1))
						a2 = w2_opt / ((eps+b2)*(eps+b2))
						a3 = w3_opt / ((eps+b3)*(eps+b3))
					ENDIF
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
					IF (Mapped) THEN
!						Mapping the WENO weights
						a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
						a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
						a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!						Making the weights convex
						w1 = a1 / (a1 + a2 + a3)
						w2 = a2 / (a1 + a2 + a3)
						w3 = a3 / (a1 + a2 + a3)
					ENDIF
					Q_Face_Left(i,n) = w1*f1 + w2*f2 + w3*f3
				ENDDO !i loop
			ENDIF

			a = 0
			b = 0
			c = 0
			r = 0

!			Right Biased Interpolation

!			Left Boundary (Non-compact)
			a(imin+1) = 0.0
			b(imin+1) = 1.0
			c(imin+1) = 0.0
			r(imin+1) = (1.0_8/3.0_8)*Q_Cell(imin+3,n) - (7.0_8/6.0_8)*Q_Cell(imin+2,n) + (11.0_8/6.0_8)*Q_Cell(imin+1,n)
!			Right Boundary (Non-compact)
			a(imax) = 0.0
			b(imax) = 1.0
			c(imax) = 0.0
			r(imax) = (1.0_8/3.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) - (1.0_8/6.0_8)*Q_Cell(imax-2,n)
			a(imax-1) = 0.0
			b(imax-1) = 1.0
			c(imax-1) = 0.0
			r(imax-1) = -(1.0_8/6.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) + (1.0_8/3.0_8)*Q_Cell(imax-2,n)

!			Redundant (compact)
			a(imin) = 0.0
			b(imin) = 1.0
			c(imin) = 0.0
			r(imin) = 0.0
!			Left Boundary (compact)
!			a(imin+1) = 0.0
!			b(imin+1) = 2.0/3.0
!			c(imin+1) = 1.0/3.0
!			r(imin+1) = (1.0/6.0) * (Q_Cell(imin,n) + 5.0*Q_Cell(imin+1,n))
!			Right Boundary (Compact)
!			a(imax-1) = 1.0/10.0
!			b(imax-1) = 6.0/10.0
!			c(imax-1) = 3.0/10.0
!			r(imax-1) = (1.0/30.0) * (10.0*Q_Cell(imax-2,n) + 19.0*Q_Cell(imax-1,n) + Q_Cell(imax,n))
!			a(imax) = 1.0/3.0
!			b(imax) = 2.0/3.0
!			c(imax) = 0.0
!			r(imax) = (1.0/6.0) * (Q_Cell(imax,n) + 5.0*Q_Cell(imax-1,n))

!			Interior
			DO i = imin+2, imax-2
!				Defining stencil points
				Q_m1 = Q_Cell(i-2, n)
				Q_c = Q_Cell(i-1, n)
				Q_p1 = Q_Cell(i, n)
				Q_p2 = Q_Cell(i+1, n)
				Q_p3 = Q_Cell(i+2, n)
!				Candidate Stencils
				f1 = (Q_p2 + 5.0_8*Q_p1) / 6.0_8
				f2 = (5.0_8*Q_p1 + Q_c) / 6.0_8
				f3 = (Q_p1 + 5.0_8*Q_c) / 6.0_8
!				Smoothness Indicators
				b1 = (13.0_8/12.0_8) * (Q_p3-2.0_8*Q_p2+Q_p1) *(Q_p3-2.0_8*Q_p2+Q_p1) &
					+ (1.0_8/4.0_8) * (Q_p3-4.0_8*Q_p2+3.0_8*Q_p1) * (Q_p3-4.0_8*Q_p2+3.0_8*Q_p1)
				b2 = (13.0_8/12.0_8) * (Q_p2-2.0_8*Q_p1+Q_c) * (Q_p2-2.0_8*Q_p1+Q_c) &
					+ (1.0_8/4.0_8) * (Q_p2-Q_c) * (Q_p2-Q_c)
				b3 = (13.0_8/12.0_8) * (Q_p1-2.0_8*Q_c+Q_m1) * (Q_p1-2.0_8*Q_c+Q_m1) &
					+ (1.0_8/4.0_8) * (3.0_8*Q_p1-4.0_8*Q_c+Q_m1) *(3.0_8*Q_p1-4.0_8*Q_c+Q_m1)
				
				IF (Variable_Eps) THEN
					term = MIN(b1,b2,b3) / (vareps_min + MAX(b1,b2,b3) - MIN(b1,b2,b3))
					eps = vareps_max * MIN(1.0_8,term) + vareps_min
				ENDIF

				IF (Borges) THEN
!					Borges definition of tau
					tau = abs (b1 - b3)
				ELSEIF (Carpenter) THEN
!					Yamaleev-Carpenter definition of tau
					tau = (Q_p3-4.0_8*Q_p2+6.0_8*Q_p1-4.0_8*Q_c+Q_m1) &
					    * (Q_p3-4.0_8*Q_p2+6.0_8*Q_p1-4.0_8*Q_c+Q_m1)
				ELSE
					tau = 0.0_8
				ENDIF
!				Optimal Weights
				w1_opt = 0.2_8
				w2_opt = 0.5_8
				w3_opt = 0.3_8
!				WENO Limiting
				IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
					a1 = w1_opt * (1.0_8 + (tau/(eps+b1))**IS_p)
					a2 = w2_opt * (1.0_8 + (tau/(eps+b2))**IS_p)
					a3 = w3_opt * (1.0_8 + (tau/(eps+b3))**IS_p)
				ELSE
					a1 = w1_opt / ((eps+b1)**IS_p)
					a2 = w2_opt / ((eps+b2)**IS_p)
					a3 = w3_opt / ((eps+b3)**IS_p)
				ENDIF
!				Making the weights convex
				w1 = a1 / (a1 + a2 + a3)
				w2 = a2 / (a1 + a2 + a3)
				w3 = a3 / (a1 + a2 + a3)
				IF (Mapped) THEN
!					Mapping the WENO weights
					a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
					a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
					a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
				ENDIF

				r(i) = w1*f1 + w2*f2 + w3*f3
				a(i) = (1.0_8/3.0_8) * w3
				b(i) = (1.0_8/3.0_8) * w1 + (2.0_8/3.0_8) * (w2 + w3)
				c(i) = (2.0_8/3.0_8) * w1 + (1.0_8/3.0_8) * w2
			ENDDO !i-loop
			CALL TRIDIAGONAL_SOLVER(imin, imax, a, b, c, r, tridiag_flag)
			IF (tridiag_flag) THEN
				DO i = imin+1, imax
					Q_Face_Right(i,n) = r(i)
				ENDDO !i loop
			ELSE
!				Falling back on WENO5 if singular system for the
!				compact scheme
				Q_Face_Right(imin+1,n) = (1.0_8/3.0_8)*Q_Cell(imin+3,n) - (7.0_8/6.0_8)*Q_Cell(imin+2,n) + (11.0_8/6.0_8)*Q_Cell(imin+1,n)
				Q_Face_Right(imax,n) = (1.0_8/3.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) - (1.0_8/6.0_8)*Q_Cell(imax-2,n)
				Q_Face_Right(imax-1,n) = -(1.0_8/6.0_8)*Q_Cell(imax,n) + (5.0_8/6.0_8)*Q_Cell(imax-1,n) + (1.0_8/3.0_8)*Q_Cell(imax-2,n)
				DO i = imin+2, imax-2
!					Defining stencil points
					Q_m2 = Q_Cell(i+2, n)
					Q_m1 = Q_Cell(i+1, n)
					Q_c  = Q_Cell(i, n)
					Q_p1 = Q_Cell(i-1, n)
					Q_p2 = Q_Cell(i-2, n)
!					Candidate Stencils
					f1 =  (1.0_8/3.0_8)*Q_c + (5.0_8/6.0_8)*Q_p1 - (1.0_8/6.0_8)*Q_p2
					f2 = -(1.0_8/6.0_8)*Q_m1 + (5.0_8/6.0_8)*Q_c + (1.0_8/3.0_8)*Q_p1
					f3 =  (1.0_8/3.0_8)*Q_m2 - (7.0_8/6.0_8)*Q_m1 + (11.0_8/6.0_8)*Q_c
!					Smoothness indicators
					b1 = (13.0_8/12.0_8) * (Q_c-2.0_8*Q_p1+Q_p2)*(Q_c-2.0_8*Q_p1+Q_p2) &
						+ (1.0_8/4.0_8) * (3.0_8*Q_c-4.0_8*Q_p1+Q_p2)*(3.0_8*Q_c-4.0_8*Q_p1+Q_p2)
					b2 = (13.0_8/12.0_8) * (Q_m1-2.0_8*Q_c+Q_p1)*(Q_m1-2.0_8*Q_c+Q_p1) &
						+ (1.0_8/4.0_8) * (Q_m1-Q_p1)*(Q_m1-Q_p1)
					b3 = (13.0_8/12.0_8) * (Q_m2-2.0_8*Q_m1+Q_c)*(Q_m2-2.0_8*Q_m1+Q_c) &
						+ (1.0_8/4.0_8) * (Q_m2-4.0_8*Q_m1+3.0_8*Q_c)*(Q_m2-4.0_8*Q_m1+3.0_8*Q_c)
					IF (Borges) THEN
!						Borges definition of tau
						tau = abs (b1 - b3)
					ELSEIF (Carpenter) THEN
!						Yamaleev-Carpenter definition of tau
						tau = (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2) &
						    * (Q_m2-4.0_8*Q_m1+6.0_8*Q_c-4.0_8*Q_p1+Q_p2)
					ELSE
						tau = 0.0_8
					ENDIF
!					Optimal Weights
					w1_opt = 0.3
					w2_opt = 0.6
					w3_opt = 0.1
!					WENO limiting
					IF (Borges .OR. Carpenter .OR. No_Limiting) THEN
						a1 = w1_opt * (1.0_8 + tau / (eps+b1))
						a2 = w2_opt * (1.0_8 + tau / (eps+b2))
						a3 = w3_opt * (1.0_8 + tau / (eps+b3))
					ELSE
						a1 = w1_opt / ((eps+b1)*(eps+b1))
						a2 = w2_opt / ((eps+b2)*(eps+b2))
						a3 = w3_opt / ((eps+b3)*(eps+b3))
					ENDIF
!					Making the weights convex
					w1 = a1 / (a1 + a2 + a3)
					w2 = a2 / (a1 + a2 + a3)
					w3 = a3 / (a1 + a2 + a3)
					IF (Mapped) THEN
!						Mapping the WENO weights
						a1 = w1 * (w1_opt + w1_opt*w1_opt - 3.0_8*w1_opt*w1 + w1*w1) / (w1_opt*w1_opt + w1*(1.0_8-2.0_8*w1_opt))
						a2 = w2 * (w2_opt + w2_opt*w2_opt - 3.0_8*w2_opt*w2 + w2*w2) / (w2_opt*w2_opt + w2*(1.0_8-2.0_8*w2_opt))
						a3 = w3 * (w3_opt + w3_opt*w3_opt - 3.0_8*w3_opt*w3 + w3*w3) / (w3_opt*w3_opt + w3*(1.0_8-2.0_8*w3_opt))
!						Making the weights convex
						w1 = a1 / (a1 + a2 + a3)
						w2 = a2 / (a1 + a2 + a3)
						w3 = a3 / (a1 + a2 + a3)
					ENDIF
	
					Q_Face_Right(i,n) = w1*f1 + w2*f2 + w3*f3
				ENDDO !i loop
			ENDIF

			DEALLOCATE(a,b,c,r)

!			Overset correction
!			DO i = imin+2, imax-2
!				blank_fac  = iba(i)*iba(i-1)*iba(i+1)*iba(i-2)*iba(i+2)
!				Q_Face_Left(i,n) = Q_Face_Left(i,n)*blank_fac+(1.-blank_fac)*Q_Cell(i,n)
!				Q_Face_Right(i,n) = Q_Face_Right(i,n)*blank_fac+(1.-blank_fac)*Q_Cell(i,n)
!			ENDDO !i loop
		ENDDO	!n loop

	END SUBROUTINE CRWENO5_SCALAR

END MODULE WENO_ROUTINES
