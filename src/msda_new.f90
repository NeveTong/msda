! --------------------------------------------------
SUBROUTINE msda(obj,nk,nvars,sigma,delta,pf,dfmax,pmax,nlam,flmin,ulam,&
        eps,maxit,sml,verbose,nalam,theta,m,ntheta,alam,npass,jerr)
! --------------------------------------------------
    IMPLICIT NONE
    ! - - - arg types - - -
    DOUBLE PRECISION, PARAMETER :: big=9.9D30
    DOUBLE PRECISION, PARAMETER :: mfl=1.0D-6
    DOUBLE PRECISION, PARAMETER :: abs_eps=1.0D-12
    INTEGER, PARAMETER :: mnlam=6
    INTEGER, PARAMETER :: rsfreq=10
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
    DOUBLE PRECISION::dev_new
    DOUBLE PRECISION::dev_old
    DOUBLE PRECISION::d1n
    DOUBLE PRECISION::d2n
    DOUBLE PRECISION::d3n
    DOUBLE PRECISION::d1o
    DOUBLE PRECISION::d2o
    DOUBLE PRECISION::d3o
    INTEGER::i_rs
    INTEGER::j_rs
! - - - begin - - -
    IF(maxval(pf) <= 0.0D0) THEN
        jerr=10000
        RETURN
    ENDIF
    pf=max(0.0D0,pf)
! - - - some initial setup - - -
    mnl=min(mnlam,nlam)
    r=delta
    thetanew=0.0D0
    thetaold=0.0D0
    dev=0.0D0
    m=0
    mm=0
    npass=0
    ni=0
    theta=0.0D0
    obj=0.0D0
    alam=0.0D0
    ntheta=0
    nalam=0
    jerr=0
! --------- lambda loop ----------------------------
    IF(flmin < 1.0D0) THEN
        flmin=max(mfl,flmin)
        alf=flmin**(1.0D0/(nlam-1.0D0))
    ENDIF
    DO l=1,nlam

        dev=0.0D0

        IF(flmin >= 1.0D0) THEN
            al=ulam(l)
        ELSE
            IF(l > 2) THEN
                al=al*alf
            ELSE IF(l==1) THEN
                al=big
            ELSE IF(l==2) THEN
                al=0.0D0
                DO j=1,nvars
                    IF(pf(j) > 0.0D0) THEN
                        u=delta(:,j)
                        v=sqrt(dot_product(u,u))
                        al=max(al,v/pf(j))
                    ENDIF
                END DO
                al=al*alf
            ENDIF
        ENDIF

        psr=0

! --------- outer loop ----------------------------
        OUTER: DO
            ! F5: full warm-start so all active variables are captured
            thetaold=thetanew

! --middle loop-------------------------------------
            MIDDLE: DO
                npass=npass+1
                psr=psr+1
                dif=0.0D0
                dev_tmp=dev

                DO k=1,nvars
                    thetatmp=thetanew(:,k)
                    u=r(:,k)/sigma(k,k)+thetatmp
                    unorm=sqrt(dot_product(u,u))
                    v=unorm-al*pf(k)/sigma(k,k)
                    IF(v > 0.0D0) THEN
                        thetanew(:,k)=v*u/unorm
                    ELSE
                        thetanew(:,k)=0.0D0
                    ENDIF
                    d=thetanew(:,k)-thetatmp
                    theta_sum=thetanew(:,k)+thetatmp
                    IF(any(d /= 0.0D0)) THEN
                        dif=max(dif,maxval(abs(d)))
                        loss_diff=sum(d*(0.5D0*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                        penalty_diff=al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                            -sqrt(dot_product(thetatmp,thetatmp)))
                        dev=dev+loss_diff+penalty_diff
                        ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                        r=r-ab
                        IF(mm(k)==0) THEN
                            ni=ni+1
                            IF(ni > pmax) EXIT MIDDLE
                            mm(k)=ni
                            m(ni)=k
                        ENDIF
                    ENDIF
                ENDDO

                ! F1+F2: periodic exact resync of r and dev
                IF(psr >= rsfreq) THEN
                    psr=0
                    r=delta
                    DO j_rs=1,nvars
                        DO i_rs=1,nvars
                            IF(any(thetanew(:,i_rs) /= 0.0D0)) THEN
                                r(:,j_rs)=r(:,j_rs)-sigma(i_rs,j_rs)*thetanew(:,i_rs)
                            ENDIF
                        ENDDO
                    ENDDO
                    d1n=0.0D0; d2n=0.0D0; d3n=0.0D0
                    DO jj=1,nk
                        d1n=d1n+dot_product(matmul(thetanew(jj,:),sigma),thetanew(jj,:))
                        d2n=d2n+dot_product(thetanew(jj,:),delta(jj,:))
                    ENDDO
                    DO k=1,nvars
                        d3n=d3n+pf(k)*sqrt(dot_product(thetanew(:,k),thetanew(:,k)))
                    ENDDO
                    dev=0.5D0*d1n-d2n+al*d3n
                ENDIF

                ! F3+F4: combined convergence test with absolute floor
                IF(abs(dev_tmp) > abs_eps) THEN
                    IF(abs(dev-dev_tmp)/abs(dev_tmp) < sml .AND. dif < eps) EXIT MIDDLE
                ELSE
                    IF(dif < eps) EXIT MIDDLE
                ENDIF
                IF(ni > pmax) EXIT MIDDLE
                IF(npass > maxit) THEN
                    jerr=-l
                    RETURN
                ENDIF

! --inner loop----------------------
                INNER: DO
                    npass=npass+1
                    psr=psr+1
                    dif=0.0D0
                    dev_tmp=dev

                    DO j=1,ni
                        k=m(j)
                        thetatmp=thetanew(:,k)
                        u=r(:,k)/sigma(k,k)+thetatmp
                        unorm=sqrt(dot_product(u,u))
                        v=unorm-al*pf(k)/sigma(k,k)
                        IF(v > 0.0D0) THEN
                            thetanew(:,k)=v*u/unorm
                        ELSE
                            thetanew(:,k)=0.0D0
                        ENDIF
                        d=thetanew(:,k)-thetatmp
                        theta_sum=thetanew(:,k)+thetatmp
                        IF(any(d /= 0.0D0)) THEN
                            dif=max(dif,maxval(abs(d)))
                            loss_diff=sum(d*(0.5D0*theta_sum-r(:,k)-sigma(k,k)*thetatmp))
                            penalty_diff=al*pf(k)*(sqrt(dot_product(thetanew(:,k),thetanew(:,k))) &
                                -sqrt(dot_product(thetatmp,thetatmp)))
                            dev=dev+loss_diff+penalty_diff
                            ab=spread(sigma(:,k),dim=1,ncopies=nk)*spread(d,dim=2,ncopies=nvars)
                            r=r-ab
                        ENDIF
                    ENDDO

                    ! F1+F2: periodic exact resync inside inner loop
                    IF(psr >= rsfreq) THEN
                        psr=0
                        r=delta
                        DO j_rs=1,nvars
                            DO i_rs=1,nvars
                                IF(any(thetanew(:,i_rs) /= 0.0D0)) THEN
                                    r(:,j_rs)=r(:,j_rs)-sigma(i_rs,j_rs)*thetanew(:,i_rs)
                                ENDIF
                            ENDDO
                        ENDDO
                        d1n=0.0D0; d2n=0.0D0; d3n=0.0D0
                        DO jj=1,nk
                            d1n=d1n+dot_product(matmul(thetanew(jj,:),sigma),thetanew(jj,:))
                            d2n=d2n+dot_product(thetanew(jj,:),delta(jj,:))
                        ENDDO
                        DO k=1,nvars
                            d3n=d3n+pf(k)*sqrt(dot_product(thetanew(:,k),thetanew(:,k)))
                        ENDDO
                        dev=0.5D0*d1n-d2n+al*d3n
                    ENDIF

                    ! F3+F4: combined convergence test with absolute floor
                    IF(abs(dev_tmp) > abs_eps) THEN
                        IF(abs(dev-dev_tmp)/abs(dev_tmp) < sml .AND. dif < eps) EXIT INNER
                    ELSE
                        IF(dif < eps) EXIT INNER
                    ENDIF
                    IF(npass > maxit) THEN
                        jerr=-l
                        RETURN
                    ENDIF

                END DO INNER
            END DO MIDDLE

            IF(ni > pmax) EXIT OUTER

!--- this is the final check ------------------------
            vrg=1
            DO j=1,ni
                IF(maxval(abs(thetanew(:,m(j))-thetaold(:,m(j)))) >= eps) THEN
                    vrg=0
                    EXIT
                ENDIF
            ENDDO
            IF(vrg==1) EXIT OUTER

            ! compute objective exactly for outer convergence test
            d1n=0.0D0; d2n=0.0D0
            d1o=0.0D0; d2o=0.0D0
            DO jj=1,nk
                d1n=d1n+dot_product(matmul(thetanew(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetanew(jj,m(1:ni)))
                d2n=d2n+dot_product(thetanew(jj,m(1:ni)),delta(jj,m(1:ni)))
                d1o=d1o+dot_product(matmul(thetaold(jj,m(1:ni)),sigma(m(1:ni),m(1:ni))),thetaold(jj,m(1:ni)))
                d2o=d2o+dot_product(thetaold(jj,m(1:ni)),delta(jj,m(1:ni)))
            ENDDO
            d3n=al*sum(pf(m(1:ni))*sqrt(sum(thetanew(:,m(1:ni))*thetanew(:,m(1:ni)),DIM=1)))
            d3o=al*sum(pf(m(1:ni))*sqrt(sum(thetaold(:,m(1:ni))*thetaold(:,m(1:ni)),DIM=1)))
            dev_new=0.5D0*d1n-d2n+d3n
            dev_old=0.5D0*d1o-d2o+d3o

            IF(verbose==1) THEN
                CALL intpr('Current Lambda',-1,l,1)
                CALL dblepr('Obj-func Jump',-1,abs(dev_new-dev_old)/(abs(dev_new)+abs_eps),1)
            ENDIF

            ! F4: guarded outer convergence check
            IF(abs(dev_new-dev_old)/(abs(dev_new)+abs_eps) < sml) EXIT OUTER

        END DO OUTER

!--- final update variable save results------------
        IF(ni > pmax) THEN
            jerr=-10000-l
            EXIT
        ENDIF
        IF(ni > 0) theta(:,1:ni,l)=thetanew(:,m(1:ni))
        me=count(maxval(abs(theta(:,1:ni,l)),dim=1) /= 0.0D0)
        IF(me > dfmax) THEN
            jerr=-20000-l
            EXIT
        ENDIF

        ! recompute dev_new cleanly for storage
        d1n=0.0D0; d2n=0.0D0; d3n=0.0D0
        DO jj=1,nk
            d1n=d1n+dot_product(matmul(thetanew(jj,:),sigma),thetanew(jj,:))
            d2n=d2n+dot_product(thetanew(jj,:),delta(jj,:))
        ENDDO
        DO k=1,nvars
            d3n=d3n+pf(k)*sqrt(dot_product(thetanew(:,k),thetanew(:,k)))
        ENDDO
        dev_new=0.5D0*d1n-d2n+al*d3n

        obj(l)=dev_new
        ntheta(l)=ni
        alam(l)=al
        nalam=l
        IF(l < mnl) CYCLE
        IF(flmin >= 1.0D0) CYCLE
    ENDDO
    RETURN
END SUBROUTINE msda
