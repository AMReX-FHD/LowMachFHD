! This is a slightly modified version of fParallel/varden_ppm/bds_quadr_unlimited

module bds_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module

  implicit none
 
  private

  public :: bds
 
contains

      subroutine bds(mla,umac,s,s_update,force,dx,dt,start_comp,num_comp)
      ! modified for having the quadratic terms as well
      ! slxx and slyy are 2nd derivatives
      ! ave is the new constant for the polynomial

      type(ml_layout), intent(in   ) :: mla
      type(multifab) , intent(in   ) :: s(:)
      type(multifab) , intent(inout) :: s_update(:)
      type(multifab) , intent(in   ) :: force(:)
      type(multifab) , intent(in   ) :: umac(:,:)
      real(kind=dp_t), intent(in   ) :: dx(:,:),dt
      integer        , intent(in   ) :: start_comp, num_comp

      type(multifab) :: ave,slx,sly,slxy,slxx,slyy,sint,sc

      real(kind=dp_t), pointer :: uadvp(:,:,:,:)
      real(kind=dp_t), pointer :: vadvp(:,:,:,:)
      real(kind=dp_t), pointer ::   sop(:,:,:,:)
      real(kind=dp_t), pointer ::   snp(:,:,:,:)
      real(kind=dp_t), pointer ::    fp(:,:,:,:)
      real(kind=dp_t), pointer ::  avep(:,:,:,:)
      real(kind=dp_t), pointer ::  slxp(:,:,:,:)
      real(kind=dp_t), pointer ::  slyp(:,:,:,:)
      real(kind=dp_t), pointer :: slxyp(:,:,:,:)
      real(kind=dp_t), pointer :: slxxp(:,:,:,:)
      real(kind=dp_t), pointer :: slyyp(:,:,:,:)
      real(kind=dp_t), pointer ::   sip(:,:,:,:)
      real(kind=dp_t), pointer ::   scp(:,:,:,:)

      integer :: dm,ng,ng_u,ng_f,n,lev,i
      integer :: lo(2),hi(2)

      ! Only worry about one level
      lev = 1

      ! These all have one ghost cell and one component
      call multifab_build( ave,mla%la(lev),1,1)
      call multifab_build( slx,mla%la(lev),1,1)
      call multifab_build( sly,mla%la(lev),1,1)
      call multifab_build(slxy,mla%la(lev),1,1)
      call multifab_build(slxx,mla%la(lev),1,1)
      call multifab_build(slyy,mla%la(lev),1,1)

      ! These has two ghost cells and one component
      call multifab_build(sint,mla%la(lev),1,2)

      ! This has one ghost cell and four components
      call multifab_build(  sc,mla%la(lev),4,1)

      ng = s(1)%ng 
      ng_u = s_update(1)%ng
      ng_f = force(1)%ng
      dm = mla%dim

      do i = 1, nfabs(s(lev))
         uadvp => dataptr(umac(lev,1), i)
         vadvp => dataptr(umac(lev,2), i)

         sop   => dataptr(s(lev) , i)
         snp   => dataptr(s_update(lev), i)

         fp    => dataptr(force(lev), i)

         avep  => dataptr(ave , i)
         slxp  => dataptr(slx , i)
         slyp  => dataptr(sly , i)
         slxyp => dataptr(slxy, i)
         slxxp  => dataptr(slxx , i)
         slyyp  => dataptr(slyy , i)
         sip   => dataptr(sint, i)
         scp   => dataptr(sc  , i)
         lo =  lwb(get_box(s(lev), i))
         hi =  upb(get_box(s(lev), i))
 
         select case (dm)
            case (2)
             do n = start_comp, start_comp+num_comp-1
                 call bdsslope_2d(lo, hi, sop(:,:,1,n), ng, &
                                  avep(:,:,1,1), slxp(:,:,1,1), slyp(:,:,1,1), &
                                  slxyp(:,:,1,1), slxxp(:,:,1,1), slyyp(:,:,1,1), &
                                  sip(:,:,1,1), scp(:,:,1,:), dx(lev,:)) 

                 call  bdsconc_2d(lo, hi, sop(:,:,1,n), ng, snp(:,:,1,n), ng_u, &
                                  fp(:,:,1,n), ng_f, &
                                  avep(:,:,1,1), slxp(:,:,1,1), slyp(:,:,1,1), slxyp(:,:,1,1), &
                                  slxxp(:,:,1,1), slyyp(:,:,1,1), &
                                   sip(:,:,1,1),  scp(:,:,1,:), &
                                  uadvp(:,:,1,1), vadvp(:,:,1,1), dx(lev,:), dt)
              end do
            case (3)
               call parallel_abort("BDS advection not supported in 3D for now")  
         end select
      end do

      call multifab_destroy(ave)
      call multifab_destroy(slx)
      call multifab_destroy(sly)
      call multifab_destroy(slxy)
      call multifab_destroy(slxx)
      call multifab_destroy(slyy)
      call multifab_destroy(sint)
      call multifab_destroy(sc)

      end subroutine bds

      subroutine bdsslope_2d(lo,hi,s,ng,ave,slx,sly,slxy,slxx,slyy,sint,sc,dx)

      integer        , intent(in   ) :: lo(:), hi(:), ng
      real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:)
      real(kind=dp_t), intent(inout) ::  ave(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) ::  slx(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) ::  sly(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) :: slxy(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) :: slxx(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) :: slyy(lo(1)-1:,lo(2)-1:)
      real(kind=dp_t), intent(inout) :: sint(lo(1)-2:,lo(2)-2:)
      real(kind=dp_t), intent(inout) ::   sc(lo(1)-1:,lo(2)-1:,:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      real(kind=dp_t), allocatable :: diff(:,:)
      real(kind=dp_t), allocatable :: smin(:,:)
      real(kind=dp_t), allocatable :: smax(:,:)
      real(kind=dp_t), allocatable :: sumdif(:)
      real(kind=dp_t), allocatable :: sgndif(:)
      integer        , allocatable :: kdp(:)

      real(kind=dp_t) :: hx,hy,sumloc,redfac,redmax,div
      real(kind=dp_t) :: eps
      integer         :: inc1, inc2, inc3, inc4
      integer         :: i,j,k,ll,is,ie,js,je

      allocate(  diff(lo(1)-1:hi(1)+1,4))
      allocate(  smin(lo(1)-1:hi(1)+1,4))
      allocate(  smax(lo(1)-1:hi(1)+1,4))
      allocate(sumdif(lo(1)-1:hi(1)+1  ))
      allocate(sgndif(lo(1)-1:hi(1)+1  ))
      allocate(   kdp(lo(1)-1:hi(1)+1  ))

      hx = dx(1)
      hy = dx(2)
      is = lo(1)
      ie = hi(1)
      js = lo(2)
      je = hi(2)

      eps = 1.d-10

! estimate corner values and calculate slopes out of it

      do i = is-2,ie+1
        do j = js-2,je+1
          sint(i,j) = (s(i-1,j-1) + s(i-1,j+2) + s(i+2,j-1) + s(i+2,j+2) &
               - 7.d0*(s(i-1,j  ) + s(i-1,j+1) + s(i  ,j-1) + s(i+1,j-1) + & 
                       s(i  ,j+2) + s(i+1,j+2) + s(i+2,j  ) + s(i+2,j+1)) +  &
                49.d0*(s(i  ,j  ) + s(i+1,j  ) + s(i  ,j+1) + s(i+1,j+1)) ) / 144.d0
        enddo
      enddo

      do j = js-1,je+1
        do i = is-1,ie+1 

          slx(i,j) = 0.5d0*(sint(i  ,j) + sint(i  ,j-1) - &
                            sint(i-1,j) - sint(i-1,j-1) ) / hx
          sly(i,j) = 0.5d0*(sint(i  ,j) - sint(i  ,j-1) + &
                            sint(i-1,j) - sint(i-1,j-1) ) / hy
          slxy(i,j) = (sint(i,j  ) - sint(i  ,j-1) - &
                       sint(i-1,j) + sint(i-1,j-1) ) / (hx*hy)

        enddo
      enddo

! estimate 2nd derivatives at cell center (slxx and slyy are 1/2 of that) 
! and adjust constant to fit cell average

      do j = js-1,je+1
        do i = is-1,ie+1 

          slxx(i,j) = 0.5d0*( - s(i-2,j) + 12.d0*s(i-1,j) - 22.d0*s(i,j)  &
               + 12.d0*s(i+1,j) - s(i+2,j) ) / (8.d0*hx**2)
          slyy(i,j) = 0.5d0*( - s(i,j-2) + 12.d0*s(i,j-1) - 22.d0*s(i,j)  &
               + 12.d0*s(i,j+1) - s(i,j+2) ) / (8.d0*hy**2)

          ave(i,j) = s(i,j) - ( slxx(i,j)*hx**2 + slyy(i,j)*hy**2 ) / 12.d0

        enddo
      enddo
      
      deallocate(diff,smin,smax,sumdif,sgndif,kdp)
      return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! comment out these two lines to add in basic limiting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j = js-1,je+1

        do i = is-1,ie+1
          smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
          smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
          smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
          smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
          smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

          sc(i,j,4) = s(i,j) + 0.5d0*(hx*slx(i,j) + hy*sly(i,j))  &
                      + 0.25d0*hx*hy*slxy(i,j)
          sc(i,j,3) = s(i,j) + 0.5d0*(hx*slx(i,j) - hy*sly(i,j))  &
                      - 0.25d0*hx*hy*slxy(i,j)
          sc(i,j,2) = s(i,j) - 0.5d0*(hx*slx(i,j) - hy*sly(i,j)) &
                      - 0.25d0*hx*hy*slxy(i,j)
          sc(i,j,1) = s(i,j) - 0.5d0*(hx*slx(i,j) + hy*sly(i,j)) &
                      + 0.25d0*hx*hy*slxy(i,j)

          sc(i,j,4) = max(min(sc(i,j,4), smax(i,4)), smin(i,4))
          sc(i,j,3) = max(min(sc(i,j,3), smax(i,3)), smin(i,3))
          sc(i,j,2) = max(min(sc(i,j,2), smax(i,2)), smin(i,2))
          sc(i,j,1) = max(min(sc(i,j,1), smax(i,1)), smin(i,1))
 
         end do

        do ll = 1,3 
          do i = is-1,ie+1 
            sumloc = 0.25d0*(sc(i,j,4) + sc(i,j,3) +  &
                             sc(i,j,2) + sc(i,j,1))
            sumdif(i) = (sumloc - s(i,j))*4.d0
            sgndif(i) = sign(1.d0,sumdif(i))

            diff(i,4) = (sc(i,j,4) - s(i,j))*sgndif(i)
            diff(i,3) = (sc(i,j,3) - s(i,j))*sgndif(i)
            diff(i,2) = (sc(i,j,2) - s(i,j))*sgndif(i)
            diff(i,1) = (sc(i,j,1) - s(i,j))*sgndif(i)

            if (diff(i,1) .gt. eps) then
               inc1 = 1
            else
               inc1 = 0
            end if 

            if (diff(i,2) .gt. eps) then
               inc2 = 1
            else
               inc2 = 0
            end if 

            if (diff(i,3) .gt. eps) then
               inc3 = 1
            else
               inc3 = 0
            end if 

            if (diff(i,4) .gt. eps) then
               inc4 = 1
            else
               inc4 = 0
            end if 

            kdp(i) = inc1 + inc2 + inc3 + inc4

          enddo

          do k = 1,4 
            do i = is-1,ie+1 
              if (kdp(i).lt.1) then 
                 div = 1.d0
              else
                 div = dble(kdp(i))
              end if

              if (diff(i,k).gt.eps) then
                 redfac = sumdif(i)*sgndif(i)/div
                 kdp(i) = kdp(i)-1
              else
                 redfac = 0.d0
              end if

              if (sgndif(i) .gt. 0.d0) then
                 redmax = sc(i,j,k) - smin(i,k)
              else
                 redmax = smax(i,k) - sc(i,j,k)
              end if

              redfac = min(redfac,redmax)
              sumdif(i) = sumdif(i) - redfac*sgndif(i)
              sc(i,j,k) = sc(i,j,k) - redfac*sgndif(i)
            enddo
          enddo

        enddo

        do i = is-1,ie+1 
          slx(i,j) = 0.5d0*(sc(i,j,4) + sc(i,j,3) -  &
                            sc(i,j,1) - sc(i,j,2))/hx
          sly(i,j) = 0.5d0*(sc(i,j,4) + sc(i,j,2) -  &
                            sc(i,j,1) - sc(i,j,3))/hy
          slxy(i,j) = ( sc(i,j,1) + sc(i,j,4) &
                       -sc(i,j,2) - sc(i,j,3) ) / (hx*hy)
        enddo
      enddo

      deallocate(diff,smin,smax,sumdif,sgndif,kdp)

      end subroutine bdsslope_2d


      ! ***********************************************
      ! start of routine bdsconc_2d
      ! ***********************************************


      subroutine bdsconc_2d(lo,hi,s,ng,s_update,ng_u,force,ng_f, &
                            ave,slx,sly,slxy,slxx,slyy,sint,sc,uadv,vadv,dx,dt)

      integer        ,intent(in   ) :: lo(:), hi(:), ng, ng_u, ng_f
      real(kind=dp_t),intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:)
      real(kind=dp_t),intent(inout) ::   s_update(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t),intent(in   ) ::      force(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t),intent(in   ) ::  ave(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) ::  slx(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) ::  sly(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) :: slxy(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) :: slxx(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) :: slyy(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) :: sint(lo(1)- 2:,lo(2)- 2:)
      real(kind=dp_t),intent(in   ) ::   sc(lo(1)- 1:,lo(2)- 1:,:)
      real(kind=dp_t),intent(in   ) :: uadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) :: vadv(lo(1)- 1:,lo(2)- 1:)
      real(kind=dp_t),intent(in   ) ::   dx(:),dt

      real(kind=dp_t),allocatable ::   siphj(:,:)
      real(kind=dp_t),allocatable ::   sijph(:,:)
      real(kind=dp_t),allocatable ::    gamp(:)
      real(kind=dp_t),allocatable ::    gamm(:)
      real(kind=dp_t),allocatable ::      xm(:)
      real(kind=dp_t),allocatable ::      ym(:)
      real(kind=dp_t),allocatable ::       c(:,:)

      
      real(kind=dp_t) :: hx,hy,dt3rd,hxs,hys
      real(kind=dp_t) :: vtrans,stem,vaddif,vdif,vmult
      real(kind=dp_t) :: isign, jsign, force_local
      integer i,j,is,ie,js,je
      integer iup,jup

      real(kind=dp_t) :: uconv,vconv,divu
      real(kind=dp_t) :: u1,u2,v1,v2,uu,vv
      real(kind=dp_t) :: eps
      real(kind=dp_t), parameter :: two3rd = 2.d0/3.d0

      allocate(siphj(lo(1)  :hi(1)+1,lo(2):hi(2)  ))
      allocate(sijph(lo(1)  :hi(1)  ,lo(2):hi(2)+1))

      allocate( gamp(lo(1)-1:hi(1)+1))
      allocate( gamm(lo(1)-1:hi(1)+1))
      allocate(   xm(lo(1)-1:hi(1)+1))
      allocate(   ym(lo(1)-1:hi(1)+1))
      allocate(    c(lo(1)-1:hi(1)+1,4))

      eps = 1.d-8

      is = lo(1)
      js = lo(2)
      ie = hi(1)
      je = hi(2)
      hx = dx(1)
      hy = dx(2)

      dt3rd = dt / 3.d0

! note: the reconstructed function has now the form 
! ave + (x-xi)*slx + (y-yj)*sly + (x-xi)(y-yj)*slxy
! + (x-xi)^2*slxx + (y-yj)^2*slyy

      do j = js,je 

        ! ******************************* 
        ! calculate Gamma plus for flux F

        do i = is-1,ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j+1)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j
             jsign = 1.d0
             u2 = uadv(i+1,j)
          else 
             jup   = j+1
             jsign = -1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j+1) .gt. 0) then
                 u2 = uadv(i+1,j+1)
             end if
          end if
  
        vv = vadv(iup,j+1)

        hxs = hx*isign
        hys = hy*jsign

! quadrature rules:
! midpoint for const and linear terms
! average of midpoint of edges for bilin and quadr terms

        gamp(i) = ave(iup,jup)+     &
         (hxs*.5 - (u1+u2)*dt/3.d0)*slx(iup,jup) +   &
         (hys*.5 -    vv*dt/3.d0)*sly(iup,jup) +    &
         (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
         vv*(2.*u2+u1)*dt*dt)*slxy(iup,jup)/12.d0 + &
         ((hxs-dt*u1)**2 + (hxs-dt*u2)**2 + (hxs-dt*(u1+u2))**2)*slxx(iup,jup)/12.d0 + &
         (hys**2+2.*(hys-dt*vv)**2)*slyy(iup,jup)/12.d0

        enddo

        ! end of calculation of Gamma plus for flux F
        ! ****************************************


        ! *****************************************
        ! calculate Gamma minus for flux F

        do i = is-1,ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j-1
             jsign = 1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j-1) .gt. 0) then
                u2 = uadv(i+1,j-1)
             end if
          else 
             jup   = j
             jsign = -1.d0
             u2 = uadv(i+1,j)
          end if


        vv = vadv(iup,j)

        hxs = hx*isign
        hys = hy*jsign

        gamm(i) = ave(iup,jup)+     &
         (hxs*.5 - (u1+u2)*dt/3.d0)*slx(iup,jup) +   &
         (hys*.5 -    vv*dt/3.d0)*sly(iup,jup) +    &
         (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
         vv*(2.*u2+u1)*dt*dt)*slxy(iup,jup)/12.d0 + &
         ((hxs-dt*u1)**2 + (hxs-dt*u2)**2 + (hxs-dt*(u1+u2))**2)*slxx(iup,jup)/12.d0 + &
         (hys**2+2.*(hys-dt*vv)**2)*slyy(iup,jup)/12.d0


        enddo

        ! end of calculation of Gamma minus for flux F
        ! ****************************************


        ! *********************************
        ! calculate siphj

        do i = is-1, ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
             force_local = force(i,j)
          else
             iup   = i+1
             isign = -1.d0
             force_local = force(i+1,j)
          end if

! gamm and gamp are updated to handle quadratic terms
! stem is updated to handle quadratic terms (using a 1d 2-point Gaussian formula)
! this automatically should give `correct' volume integrals (in the non-div free case)

          hxs = hx*isign

          vdif = 0.5d0*dt*(vadv(iup,j+1)*gamp(i) -  &
                           vadv(iup,j)*gamm(i) ) / hy
          stem = ave(iup,j) + (hxs - uadv(i+1,j)*dt)*0.5d0*slx(iup,j) + &
               0.5d0*slxx(iup,j)*( (hxs/2. - uadv(i+1,j)*dt*(1.+sqrt(3.))/(2.*sqrt(3.)))**2 + &
                 (hxs/2. + uadv(i+1,j)*dt*(1.-sqrt(3.))/(2.*sqrt(3.)))**2   ) + &
               slyy(iup,j)*hy*hy/12.d0
          vaddif = stem*0.5d0*dt*(uadv(iup+1,j) - uadv(iup,j))/hx
          divu =  &
            (uadv(iup+1,j)-uadv(iup,j))/hx +  &
            (vadv(iup,j+1)-vadv(iup,j))/hy 
!          siphj(i+1,j) = stem - vdif - vaddif + 0.5d0*dt*(stem*divu + force_local)
          siphj(i+1,j) = stem - vdif - vaddif + 0.5d0*dt*force_local

        enddo
      enddo

      ! end of calculation of siphj
      ! *************************************

      do j = js-1,je 

        ! ********************************** 
        ! calculate Gamma plus for flux G


        do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i+1,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i
             isign = 1.d0
             v2 = vadv(i,j+1)
          else
             iup   = i+1
             isign = -1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i+1,j+1) .gt. 0) then
                v2 = vadv(i+1,j+1)
             end if
          end if


          uu = uadv(i+1,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamp(i) = ave(iup,jup)+ &
          (hys*.5 - (v1+v2)*dt/3.)*sly(iup,jup) +   &
          (hxs*.5 - uu*dt/3.)*slx(iup,jup) + &
          (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
          (2.*v2+v1)*uu*dt*dt)*slxy(iup,jup)/12.d0 + &
          (hxs**2 + 2.*(hxs-uu*dt)**2)*slxx(iup,jup)/12.d0 + &
          ((hys-v1*dt)**2+(hys-v2*dt)**2+(hys-(v1+v2)*dt)**2)*slyy(iup,jup)/12.d0

        enddo

        ! end of calculation of Gamma plus for flux G
        ! ****************************************


        ! *****************************************
        ! calculate Gamma minus for flux G

        do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i-1
             isign = 1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i-1,j+1) .gt. 0) then
                v2 = vadv(i-1,j+1)
             end if
          else
             iup   = i
             isign = -1.d0
             v2 = vadv(i,j+1)
          end if


          uu = uadv(i,jup)       

          hxs = hx*isign
          hys = hy*jsign

         gamm(i) = ave(iup,jup) +    &
          (hys*.5 - (v1+v2)*dt/3.)*sly(iup,jup) +    &
          (hxs*.5 - uu*dt/3.)*slx(iup,jup) +   &
          (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
          (2.*v2+v1)*uu*dt*dt)*slxy(iup,jup)/12.d0 + &
          (hxs**2 + 2.*(hxs-uu*dt)**2)*slxx(iup,jup)/12.d0 + &
          ((hys-v1*dt)**2+(hys-v2*dt)**2+(hys-(v1+v2)*dt)**2)*slyy(iup,jup)/12.d0
          
        enddo

        ! end of calculation of Gamma minus for flux G
        ! ****************************************


        ! *********************************
        ! calculate sijph

        do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
             force_local = force(i,j)
          else
             jup   = j+1
             jsign = -1.d0
             force_local = force(i,j+1)
          end if

          hys = hy*jsign

          vdif = 0.5d0*dt* &
            (uadv(i+1,jup)*gamp(i)-uadv(i,jup)*gamm(i))/hx
          stem = ave(i,jup) + (hys - vadv(i,j+1)*dt)*0.5d0*sly(i,jup) + &
               0.5d0*slyy(i,jup)*( (hys/2. - vadv(i,j+1)*dt*(1.+sqrt(3.))/(2.*sqrt(3.)))**2 + &
                 (hys/2. + vadv(i,j+1)*dt*(1.-sqrt(3.))/(2.*sqrt(3.)))**2   ) + &
               slxx(i,jup)*hx*hx/12.d0
          vaddif = stem*0.5d0*dt*(vadv(i,jup+1) - vadv(i,jup))/hy
          divu =  (uadv(i+1,jup)-uadv(i,jup))/hx +  &
                  (vadv(i,jup+1)-vadv(i,jup))/hy 
!         sijph(i,j+1) = stem - vdif - vaddif + 0.5d0*dt*(stem*divu + force_local)
          sijph(i,j+1) = stem - vdif - vaddif + 0.5d0*dt*force_local

        enddo
      enddo

      ! end of calculation of sijph
      ! *************************************

        do j = js,je 
        do i = is,ie 

          s_update(i,j) = s_update(i,j) -(  &
               (siphj(i+1,j)*uadv(i+1,j)-siphj(i,j)*uadv(i,j))/hx +  &
               (sijph(i,j+1)*vadv(i,j+1)-sijph(i,j)*vadv(i,j))/hy)

        enddo
        enddo

      deallocate(siphj,sijph,gamp,gamm,xm,ym,c)

      end subroutine bdsconc_2d

end module bds_module
