! --------------------------------------------------
SUBROUTINE msda(obj,nk,nvars,sigma,delta,pf,dfmax,pmax,nlam,flmin,ulam,&
        eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
! Fixes vs original:
!   F1. Periodic exact recomputation of residual r      (prevents incremental drift)
!   F2. Periodic exact recomputation of objective dev   (prevents float accumulation)
!   F3. Combined convergence guard: require BOTH dif<eps AND rel-dev<sml to exit
!   F4. Absolute-value floor (1D-12) in all relative-change tests (guards dev~0)
!   F6. All floating-point literals promoted to double precision (D0 suffix)
!
! NOT changed:
!   Active-set warm-start logic is preserved exactly as in the original:
!   thetaold updated only over m(1:ni), maintaining monotone increasing support.
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - parameters - - -
    DOUBLE PRECISION, PARAMETER :: big    = 9.9D30
    DOUBLE PRECISION, PARAMETER :: mfl    = 1.0D-6
    DOUBLE PRECISION, PARAMETER :: abseps = 1.0D-12  ! floor for relative tests (F4)
    INTEGER,          PARAMETER :: mnlam  = 6
    INTEGER,          PARAMETER :: rsfreq = 10        ! resync period (F1,F2)
    ! - - - arg types - - -
    INTEGER::mnl
    INTEGER::nk
    INTEGER::nvars
    INTEGER::dfmax
    INTEGER::pmax
    INTEGER::nlam
    INTEGER::nalam
    INTEGER::npass
    INTEGER::jerr
    INTEGER::verbose
    INTEGER::maxit
    INTEGER::m(pmax)
    INTEGER::ntheta(nlam)
    DOUBLE PRECISION::flmin
    DOUBLE PRECISION::eps
    DOUBLE PRECISION::sml
    DOUBLE PRECISION::sigma(nvars,nvars)
    DOUBLE PRECISION::delta(nk,nvars)
    DOUBLE PRECISION::pf(nvars)
    DOUBLE PRECISION::ulam(nlam)
    DOUBLE PRECISION::theta(nk,pmax,nlam)
    DOUBLE PRECISION::alam(nlam)
    DOUBLE PRECISION::obj(nlam)
    ! - - - local declarations - - -
    INTEGER::mm(nvars)
    INTEGER::k
    INTEGER::j
    INTEGER::jj
    INTEGER::l
    INTEGER::vrg
    INTEGER::ni
    INTEGER::me
    INTEGER::psr
    DOUBLE PRECISION::dif
    DOUBLE PRECISION::v
    DOUBLE PRECISION::al
    DOUBLE PRECISION::alf
    DOUBLE PRECISION::unorm
    DOUBLE PRECISION::thetanew(nk,nvars)
    DOUBLE PRECISION::thetaold(nk,nvars)
    DOUBLE PRECISION::r(nk,nvars)
    DOUBLE PRECISION::ab(nk,nvars)
    DOUBLE PRECISION::d(nk)
    DOUBLE PRECISION::theta_sum(nk)
    DOUBLE PRECISION::thetatmp(nk)
    DOUBLE PRECISION::u(nk)
    DOUBLE PRECISION::loss_diff
    DOUBLE PRECISION::penalty_diff
    DOUBLE PRECISION::dev
    DOUBLE PRECISION::dev_tmp
    DOUBLE PRECISION::tmp1_new
    DOUBLE PRECISION::tmp2_new
    DOUBLE PRECISION::dev_new
    DOUBLE PRECISION::dev1_new
    DOUBLE PRECISION::dev2_new
    DOUBLE PRECISION::dev3_new
    DOUBLE PRECISION::tmp1_old
    DOUBLE PRECISION::tmp2_old
    DOUBLE PRECISION::dev_old
    DOUBLE PRECISION::dev1_old
    DOUBLE PRECISION::dev2_old
    DOUBLE PRECISION::dev3_old
    ! resync temporaries (F1,F2)
    INTEGER::rs_i
    INTEGER::rs_j
    DOUBLE PRECISION::rs_d1
    DOUBLE PRECISION::rs_d2
    DOUBLE PRECISION::rs_d3
! - - - begin - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    mnl = Min(mnlam, nlam)
    r = delta
    thetanew=0.0D0
    thetaold=0.0D0
    dev=0.0D0
    m=0
    mm=0
    npass=0
    ni=npass
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin = Max(mfl, flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    DO l=1,nlam

        dev = 0.0D0
        psr = 0

        IF(flmin>=1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al=0.0D0
                DO j = 1,nvars
                    IF(pf(j)>0.0D0) THEN
                        u = delta(:,j)
                        v = sqrt(dot_product(u,u))
                        al=max(al,v/pf(j))
                    ENDIF
                END DO
                al=al*alf
            ENDIF
        ENDIF
! --------- outer loop ----------------------------
        DO
            IF(ni>0) thetaold(:,m(1:ni))=thetanew(:,m(1:ni))
! --middle loop-------------------------------------
            DO
                npass=npass+1
                psr=psr+1
                dif=0.0D0
                dev_tmp = dev
                DO k=1,nvars
                    thetatmp=thetanew(:,k)
                    u = r(:,k)/sigma(k,k) + thetatmp
                    unorm = sqrt(dot_product(u,u))
                    v = unorm-al*pf(k)/sigma(k,k)
                    IF(v > 0.0D0) THEN
                        thetanew(:,k) = v*u/unorm
                    ELSE
                        thetanew(:,k)=0.0D0
                    ENDIF
                    d=thetanew(:,k)-thetatmp
                    theta_sum=thetanew(:,k)+thetatmp
                    IF(any(d/=0.0D0)) THEN
                        dif=max(dif,maxval(abs(d)))
                        loss_diff = sum(d*(0.5D0*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                        penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                        - sqrt(dot_product(thetatmp,thetatmp)))
                        dev = dev + loss_diff + penalty_diff
                        ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                        r=r-ab
                        IF(mm(k)==0) THEN
                            ni=ni+1
                            IF(ni>pmax) EXIT
                            mm(k)=ni
                            m(ni)=k
                        ENDIF
                    ENDIF
                ENDDO
                ! F1+F2: periodic exact resync of r and dev
                IF(psr >= rsfreq) THEN
                    psr = 0
                    r = delta
                    DO rs_j = 1,nvars
                        DO rs_i = 1,nvars
                            IF(any(thetanew(:,rs_i) /= 0.0D0)) THEN
                                r(:,rs_j) = r(:,rs_j) - sigma(rs_i,rs_j)*thetanew(:,rs_i)
                            ENDIF
                        ENDDO
                    ENDDO
                    rs_d1=0.0D0; rs_d2=0.0D0; rs_d3=0.0D0
                    DO rs_i = 1,nk
                        rs_d1 = rs_d1 + dot_product(matmul(thetanew(rs_i,:),sigma),thetanew(rs_i,:))
                        rs_d2 = rs_d2 + dot_product(thetanew(rs_i,:),delta(rs_i,:))
                    ENDDO
                    DO rs_j = 1,nvars
                        rs_d3 = rs_d3 + pf(rs_j)*sqrt(dot_product(thetanew(:,rs_j),thetanew(:,rs_j)))
                    ENDDO
                    dev = 0.5D0*rs_d1 - rs_d2 + al*rs_d3
                ENDIF
                ! F3+F4: combined, guarded convergence test
                IF(abs(dev_tmp) > abseps) THEN
                    IF(abs(dev-dev_tmp)/abs(dev_tmp) < sml .AND. dif < eps) EXIT
                ELSE
                    IF(dif < eps) EXIT
                ENDIF
                IF(ni>pmax) EXIT
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
                ENDIF
! --inner loop----------------------
                DO
                    npass=npass+1
                    psr=psr+1
                    dif=0.0D0
                    dev_tmp = dev
                    DO j=1,ni
                        k=m(j)
                        thetatmp=thetanew(:,k)
                        u = r(:,k)/sigma(k,k) + thetatmp
                        unorm = sqrt(dot_product(u,u))
                        v = unorm-al*pf(k)/sigma(k,k)
                        IF(v > 0.0D0) THEN
                            thetanew(:,k) = v*u/unorm
                        ELSE
                            thetanew(:,k)=0.0D0
                        ENDIF
                        d=thetanew(:,k)-thetatmp
                        theta_sum=thetanew(:,k)+thetatmp
                        IF(any(d/=0.0D0)) THEN
                            dif=max(dif,maxval(abs(d)))
                            loss_diff = sum(d*(0.5D0*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                            penalty_diff = al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                            - sqrt(dot_product(thetatmp,thetatmp)))
                            dev = dev + loss_diff + penalty_diff
                            ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                            r=r-ab
                        ENDIF
                    ENDDO
                    ! F1+F2: periodic exact resync inside inner loop
                    IF(psr >= rsfreq) THEN
                        psr = 0
                        r = delta
                        DO rs_j = 1,nvars
                            DO rs_i = 1,nvars
                                IF(any(thetanew(:,rs_i) /= 0.0D0)) THEN
                                    r(:,rs_j) = r(:,rs_j) - sigma(rs_i,rs_j)*thetanew(:,rs_i)
                                ENDIF
                            ENDDO
                        ENDDO
                        rs_d1=0.0D0; rs_d2=0.0D0; rs_d3=0.0D0
                        DO rs_i = 1,nk
                            rs_d1 = rs_d1 + dot_product(matmul(thetanew(rs_i,:),sigma),thetanew(rs_i,:))
                            rs_d2 = rs_d2 + dot_product(thetanew(rs_i,:),delta(rs_i,:))
                        ENDDO
                        DO rs_j = 1,nvars
                            rs_d3 = rs_d3 + pf(rs_j)*sqrt(dot_product(thetanew(:,rs_j),thetanew(:,rs_j)))
                        ENDDO
                        dev = 0.5D0*rs_d1 - rs_d2 + al*rs_d3
                    ENDIF
                    ! F3+F4: combined, guarded exit
                    IF(abs(dev_tmp) > abseps) THEN
                        IF(abs(dev-dev_tmp)/abs(dev_tmp) < sml .AND. dif < eps) EXIT
                    ELSE
                        IF(dif < eps) EXIT
                    ENDIF
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF
                ENDDO
            ENDDO
            IF(ni>pmax) EXIT
!--- this is the final check ------------------------
            vrg=1
            DO j=1,ni
                IF(maxval(abs(thetanew(:,m(j))-thetaold(:,m(j))))>=eps) THEN
                    vrg=0
                    EXIT
                ENDIF
            ENDDO
            IF(vrg==1) EXIT
            ! test deviance loop
            dev1_new = 0.0D0
            dev2_new = 0.0D0
            dev1_old = 0.0D0
            dev2_old = 0.0D0
            DO jj = 1,nk
                tmp1_new = dot_product(MATMUL(thetanew(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetanew(jj,m(1:ni)))
                tmp2_new = dot_product(thetanew(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_new = dev1_new + tmp1_new
                dev2_new = dev2_new + tmp2_new
                tmp1_old = dot_product(MATMUL(thetaold(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetaold(jj,m(1:ni)))
                tmp2_old = dot_product(thetaold(jj,m(1:ni)),delta(jj,m(1:ni)))
                dev1_old = dev1_old + tmp1_old
                dev2_old = dev2_old + tmp2_old
            ENDDO
            dev3_new = al * sum(pf(m(1:ni)) * sqrt(sum(thetanew(:,m(1:ni)) * thetanew(:,m(1:ni)), DIM = 1)))
            dev3_old = al * sum(pf(m(1:ni)) * sqrt(sum(thetaold(:,m(1:ni)) * thetaold(:,m(1:ni)), DIM = 1)))
            dev_new = 0.5D0 * dev1_new - dev2_new + dev3_new
            dev_old = 0.5D0 * dev1_old - dev2_old + dev3_old
            IF(verbose==1) THEN
                CALL intpr('Current Lambda',-1,l,1)
                CALL dblepr('Obj-func Jump',-1,abs(dev_new-dev_old)/(abs(dev_new)+abseps),1)
            ENDIF
            ! F4: guarded outer convergence check
            IF(abs(dev_new-dev_old)/(abs(dev_new)+abseps) < sml) EXIT
            ! test deviance loop end
        ENDDO
!--- final update variable save results------------
        IF(ni>pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni>0) theta(:,1:ni,l)=thetanew(:,m(1:ni))
        me = count(maxval(abs(theta(:,1:ni,l)),dim=1)/=0.0D0)
        IF(me>dfmax) THEN
            jerr=-20000-l
            EXIT
        ENDIF
        obj(l) = dev_new
        ntheta(l)=ni
        alam(l)=al
        nalam=l
        IF (l < mnl) CYCLE
        IF (flmin >= 1.0D0) CYCLE
    ENDDO
    RETURN
END SUBROUTINE msda
