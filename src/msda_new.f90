! --------------------------------------------------
SUBROUTINE msda(obj,nk,nvars,sigma,delta,pf,dfmax,pmax,nlam,flmin,ulam, &
     eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE

    ! ---------------- arguments ----------------
    DOUBLE PRECISION, PARAMETER :: big  = 9.9D30
    DOUBLE PRECISION, PARAMETER :: mfl  = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: tiny = 1.0D-14
    INTEGER,          PARAMETER :: mnlam = 6

    INTEGER, INTENT(IN)    :: nk, nvars, dfmax, pmax, nlam, verbose, maxit
    DOUBLE PRECISION       :: flmin
    DOUBLE PRECISION, INTENT(IN) :: sigma(nvars,nvars)
    DOUBLE PRECISION, INTENT(IN) :: delta(nk,nvars)
    DOUBLE PRECISION, INTENT(IN) :: ulam(nlam)
    DOUBLE PRECISION       :: pf(nvars)
    DOUBLE PRECISION, INTENT(IN) :: eps, sml

    INTEGER, INTENT(OUT)   :: nalam, npass, jerr
    INTEGER, INTENT(OUT)   :: m(pmax), ntheta(nlam)
    DOUBLE PRECISION, INTENT(OUT) :: theta(nk,pmax,nlam), alam(nlam), obj(nlam)

    ! ---------------- locals ----------------
    INTEGER :: i, j, k, l, jj
    INTEGER :: me, ni, mnl
    INTEGER :: outer_it, active_it
    INTEGER :: violation_found
    INTEGER :: converged_active
    INTEGER :: over_pmax

    INTEGER :: mm(nvars)               ! map var -> active-set position (0 if inactive)

    DOUBLE PRECISION :: al, alf
    DOUBLE PRECISION :: dif, dd, reldev, dev_old, dev_new
    DOUBLE PRECISION :: unorm, thresh
    DOUBLE PRECISION :: sigkk(nvars), invsigkk(nvars), pf2(nvars)

    DOUBLE PRECISION :: thetanew(nk,nvars)
    DOUBLE PRECISION :: thetatmp(nk)
    DOUBLE PRECISION :: d(nk)
    DOUBLE PRECISION :: r(nk,nvars)    ! residual: delta - theta*sigma

    ! ---------------- initialize outputs ----------------
    obj    = 0.0D0
    alam   = 0.0D0
    ntheta = 0
    theta  = 0.0D0
    m      = 0
    nalam  = 0
    npass  = 0
    jerr   = 0

    ! ---------------- basic checks ----------------
    IF (MAXVAL(pf) <= 0.0D0) THEN
        jerr = 10000
        RETURN
    END IF

    DO k = 1, nvars
        pf2(k) = MAX(0.0D0, pf(k))
    END DO

    DO k = 1, nvars
        sigkk(k) = sigma(k,k)
        IF (sigkk(k) <= tiny) THEN
            jerr = 9000 + k
            RETURN
        END IF
        invsigkk(k) = 1.0D0 / sigkk(k)
    END DO

    ! ---------------- lambda path setup ----------------
    mnl = MIN(mnlam, nlam)

    IF (flmin < 1.0D0) THEN
        flmin = MAX(mfl, flmin)
        IF (nlam > 1) THEN
            alf = flmin**(1.0D0 / DBLE(nlam - 1))
        ELSE
            alf = 1.0D0
        END IF
    ELSE
        alf = 1.0D0
    END IF

    ! ---------------- warm start ----------------
    thetanew = 0.0D0

    ! ==================================================
    ! lambda loop
    ! ==================================================
    DO l = 1, nlam

        ! ----- choose lambda -----
        IF (flmin >= 1.0D0) THEN
            al = ulam(l)
        ELSE
            IF (l == 1) THEN
                al = big
            ELSEIF (l == 2) THEN
                al = 0.0D0
                DO k = 1, nvars
                    IF (pf2(k) > 0.0D0) THEN
                        unorm = SQRT(DOT_PRODUCT(delta(:,k), delta(:,k)))
                        al = MAX(al, unorm / pf2(k))
                    END IF
                END DO
                al = al * alf
            ELSE
                al = al * alf
            END IF
        END IF

        ! ----- rebuild active set from warm start -----
        mm = 0
        m  = 0
        ni = 0
        over_pmax = 0

        DO k = 1, nvars
            unorm = SQRT(DOT_PRODUCT(thetanew(:,k), thetanew(:,k)))
            IF (unorm > tiny) THEN
                ni = ni + 1
                IF (ni > pmax) THEN
                    over_pmax = 1
                    EXIT
                END IF
                m(ni)  = k
                mm(k)  = ni
            ELSE
                thetanew(:,k) = 0.0D0
            END IF
        END DO

        IF (over_pmax == 1) THEN
            jerr = -10000 - l
            EXIT
        END IF

        ! ----- residual: r = delta - theta*sigma -----
        CALL compute_residual()

        ! ----- one full sweep over all variables -----
        dif = 0.0D0
        DO k = 1, nvars
            CALL update_block(k, dd, over_pmax)
            dif = MAX(dif, dd)
            npass = npass + 1
            IF (over_pmax == 1) EXIT
            IF (npass > maxit) THEN
                jerr = -l
                RETURN
            END IF
        END DO

        IF (over_pmax == 1) THEN
            jerr = -10000 - l
            EXIT
        END IF

        ! ----- active-set loop -----
        outer_it = 0
        DO
            outer_it = outer_it + 1

            dev_old = compute_obj()

            converged_active = 0
            active_it = 0

            ! active-only sweeps
            DO
                active_it = active_it + 1
                dif = 0.0D0

                j = 1
                DO WHILE (j <= ni)
                    k = m(j)
                    CALL update_block(k, dd, over_pmax)
                    dif = MAX(dif, dd)
                    npass = npass + 1

                    IF (over_pmax == 1) EXIT
                    IF (npass > maxit) THEN
                        jerr = -l
                        RETURN
                    END IF

                    ! if k was dropped from active set, do not advance j
                    IF (j <= ni) THEN
                        IF (m(j) == k) j = j + 1
                    END IF
                END DO

                IF (over_pmax == 1) EXIT

                dev_new = compute_obj()
                reldev = ABS(dev_new - dev_old) / MAX(1.0D0, ABS(dev_old))
                dev_old = dev_new

                IF (verbose == 1) THEN
                    CALL intpr('Current Lambda', -1, l, 1)
                    CALL dblepr('Active Sweep RelObj', -1, reldev, 1)
                END IF

                IF (dif < eps .OR. reldev < sml) THEN
                    converged_active = 1
                    EXIT
                END IF
            END DO

            IF (over_pmax == 1) THEN
                jerr = -10000 - l
                EXIT
            END IF

            ! KKT check on inactive variables
            violation_found = 0
            DO k = 1, nvars
                IF (mm(k) == 0 .AND. pf2(k) > 0.0D0) THEN
                    unorm = SQRT(DOT_PRODUCT(r(:,k), r(:,k)))
                    thresh = al * pf2(k)
                    IF (unorm > thresh * (1.0D0 + 1.0D-7)) THEN
                        CALL add_active(k, over_pmax)
                        violation_found = 1
                        IF (over_pmax == 1) EXIT
                    END IF
                END IF
            END DO

            IF (over_pmax == 1) THEN
                jerr = -10000 - l
                EXIT
            END IF

            IF (converged_active == 1 .AND. violation_found == 0) EXIT

            IF (outer_it >= maxit) THEN
                jerr = -l
                RETURN
            END IF
        END DO

        IF (jerr /= 0) EXIT

        ! ----- final save for this lambda -----
        theta(:,:,l) = 0.0D0
        IF (ni > 0) theta(:,1:ni,l) = thetanew(:,m(1:ni))

        me = 0
        DO j = 1, ni
            IF (SQRT(DOT_PRODUCT(theta(:,j,l), theta(:,j,l))) > tiny) me = me + 1
        END DO

        IF (me > dfmax) THEN
            jerr = -20000 - l
            EXIT
        END IF

        obj(l)    = compute_obj()
        ntheta(l) = ni
        alam(l)   = al
        nalam     = l

        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
    END DO

    RETURN

CONTAINS

    ! --------------------------------------------------
    SUBROUTINE compute_residual()
    ! r = delta - thetanew * sigma
    ! --------------------------------------------------
        INTEGER :: col
        r = delta
        DO col = 1, nvars
            IF (ABS(sigma(1,col)) >= 0.0D0) THEN
                DO j = 1, nvars
                    IF (ABS(sigma(j,col)) > tiny) THEN
                        r(:,col) = r(:,col) - sigma(j,col) * thetanew(:,j)
                    END IF
                END DO
            END IF
        END DO
    END SUBROUTINE compute_residual

    ! --------------------------------------------------
    SUBROUTINE add_active(k, over)
    ! --------------------------------------------------
        INTEGER, INTENT(IN)    :: k
        INTEGER, INTENT(INOUT) :: over
        IF (mm(k) == 0) THEN
            ni = ni + 1
            IF (ni > pmax) THEN
                over = 1
                RETURN
            END IF
            m(ni) = k
            mm(k) = ni
        END IF
    END SUBROUTINE add_active

    ! --------------------------------------------------
    SUBROUTINE drop_active(k)
    ! --------------------------------------------------
        INTEGER, INTENT(IN) :: k
        INTEGER :: pos, lastk

        pos = mm(k)
        IF (pos <= 0) RETURN

        lastk = m(ni)
        m(pos) = lastk
        mm(lastk) = pos

        m(ni) = 0
        mm(k) = 0
        ni = ni - 1
    END SUBROUTINE drop_active

    ! --------------------------------------------------
    SUBROUTINE update_block(k, dd, over)
    ! One block coordinate update for variable k
    ! --------------------------------------------------
        INTEGER, INTENT(IN)    :: k
        DOUBLE PRECISION, INTENT(OUT) :: dd
        INTEGER, INTENT(INOUT) :: over

        DOUBLE PRECISION :: v, lamk, norm_new

        dd = 0.0D0
        IF (pf2(k) <= 0.0D0) RETURN

        thetatmp = thetanew(:,k)

        ! u = theta_k + r_k / sigma_kk
        d = r(:,k) * invsigkk(k) + thetatmp

        unorm = SQRT(DOT_PRODUCT(d, d))
        lamk  = al * pf2(k) * invsigkk(k)

        IF (unorm > lamk) THEN
            v = 1.0D0 - lamk / unorm
            thetanew(:,k) = v * d
        ELSE
            thetanew(:,k) = 0.0D0
        END IF

        d  = thetanew(:,k) - thetatmp
        dd = MAXVAL(ABS(d))

        IF (dd > tiny) THEN
            ! rank-1 residual update:
            ! r(:,j) <- r(:,j) - sigma(j,k) * d
            DO jj = 1, nvars
                IF (ABS(sigma(jj,k)) > tiny) THEN
                    r(:,jj) = r(:,jj) - sigma(jj,k) * d
                END IF
            END DO

            norm_new = SQRT(DOT_PRODUCT(thetanew(:,k), thetanew(:,k)))

            IF (norm_new > tiny) THEN
                CALL add_active(k, over)
            ELSE
                thetanew(:,k) = 0.0D0
                IF (mm(k) /= 0) CALL drop_active(k)
            END IF
        END IF
    END SUBROUTINE update_block

    ! --------------------------------------------------
    DOUBLE PRECISION FUNCTION compute_obj()
    ! Objective:
    ! 0.5 * sum_i theta(i,:) * sigma * theta(i,:)^T
    ! - sum_i theta(i,:) * delta(i,:)^T
    ! + al * sum_k pf(k) * ||theta(:,k)||_2
    ! --------------------------------------------------
        DOUBLE PRECISION :: qf, pen, normk
        DOUBLE PRECISION :: tmp(nvars)
        INTEGER :: row, col

        qf = 0.0D0
        DO row = 1, nk
            tmp = MATMUL(sigma, thetanew(row,:))
            qf = qf + 0.5D0 * DOT_PRODUCT(thetanew(row,:), tmp) &
                    - DOT_PRODUCT(thetanew(row,:), delta(row,:))
        END DO

        pen = 0.0D0
        DO col = 1, nvars
            normk = SQRT(DOT_PRODUCT(thetanew(:,col), thetanew(:,col)))
            IF (normk > tiny) pen = pen + pf2(col) * normk
        END DO

        compute_obj = qf + al * pen
    END FUNCTION compute_obj

END SUBROUTINE msda
