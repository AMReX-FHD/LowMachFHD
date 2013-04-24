module apply_precon_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use probin_gmres_module , only: mg_verbose, precon_type, stag_mg_verbosity
  use probin_common_module, only: visc_type
  use stag_mg_solver_module
  use macproject_module
  use div_and_grad_module
  use mac_applyop_module
  use convert_module
  use norm_inner_product_module
  use multifab_physbc_module

  implicit none

  private

  public :: apply_precon

contains

  ! This computes x = M^{-1} b using the approach in ./doc/PreconditionerNotes.tex
  subroutine apply_precon(mla,b_u,b_p,x_u,x_p,alpha,beta,gamma,theta,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: b_u(:,:)
    type(multifab) , intent(in   ) :: b_p(:)
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    type(multifab) , intent(in   ) :: alpha(:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: theta
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,dm
    real(kind=dp_t) :: mean_val_pres, mean_val_umac(mla%dim)

    type(multifab) ::           phi(mla%nlevel)
    type(multifab) ::       mac_rhs(mla%nlevel)
    type(multifab) ::      zero_fab(mla%nlevel)
    type(multifab) ::       x_p_tmp(mla%nlevel)
    type(multifab) :: alphainv_edge(mla%nlevel,mla%dim)
    type(multifab) ::       b_u_tmp(mla%nlevel,mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(phi(n)     ,mla%la(n),1,1)
       call multifab_build(mac_rhs(n) ,mla%la(n),1,0)
       call multifab_build(zero_fab(n),mla%la(n),1,0)
       call setval(zero_fab(n),0.d0,all=.true.)
       call multifab_build(x_p_tmp(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(alphainv_edge(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(b_u_tmp(n,i),mla%la(n),1,0,i)
       end do
    end do

    ! set the initial guess for Phi in the Poisson solve to 0
    ! set x_u = 0 as initial guess
    do n=1,nlevs
       call setval(phi(n),0.d0,all=.true.)
       do i=1,dm
          call setval(x_u(n,i),0.d0,all=.true.)
       end do
    end do

    ! 1 = projection preconditioner
    ! 2 = lower triangular preconditioner
    ! 3 = upper triangular preconditioner
    ! 4 = block diagonal preconditioner

    select case (abs(precon_type))

    case(1)  ! projection preconditioner

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! STEP 1: Solve for an intermediate state, x_u^star, using an implicit viscous solve
       !         x_u^star = A^{-1} b_u
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
        if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
           print*,""
           print*,"Begin viscous solve for x_u^star"
        end if

        ! x_u^star = A^{-1} b_u
        call stag_mg_solver(mla,alpha,beta,gamma,theta,x_u,b_u,dx,the_bc_tower)

        if (parallel_IOProcessor().and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"End viscous solve for x_u^star"
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Construct RHS for pressure Poisson problem
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! add set mac_rhs = D(x_u^star) to mac_rhs
        call compute_divu(mla,x_u,mac_rhs,dx)

        ! add b_p to mac_rhs
        do n=1,nlevs
          call multifab_plus_plus_c(mac_rhs(n),1,b_p(n),1,1,0)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Compute x_u and x_p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
          print*,""
          print*,"Begin projection"
        end if

        ! use multigrid to solve for Phi
        ! x_u^star is only passed in to get a norm for absolute residual criteria
        call macproject(mla,phi,x_u,alpha,mac_rhs,dx,the_bc_tower)

        if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
          print*,""
          print*,"End projection"
        end if

        ! compute alphainv_edge on faces by averaging and then inverting
        call average_cc_to_face_inv(nlevs,alpha,alphainv_edge,1,dm+2,1,the_bc_tower%bc_tower_array)

        ! x_u = x_u^star - (alpha I)^-1 grad Phi
        call subtract_weighted_gradp(mla,x_u,alphainv_edge,phi,dx)

        ! if precon_type = +1, or theta=0 then x_p = theta*Phi - c*beta*(mac_rhs)
        ! if precon_type = -1             then x_p = theta*Phi - c*beta*L_alpha Phi

        if ((precon_type .eq. 1) .or. (theta .eq. 0.d0)) then
          ! first set x_p = -mac_rhs 
          do n=1,nlevs
             call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)
             call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
          end do
        else
          ! first set x_p = -L_alpha Phi
          call mac_applyop(mla,x_p,phi,zero_fab,alphainv_edge,dx, &
                           the_bc_tower,bc_comp=dm+1,stencil_order=2)
        end if

        do n=1,nlevs

          if ( (abs(visc_type) .eq. 1) .or. (abs(visc_type) .eq. 2) ) then  
            ! multiply x_p by beta; x_p = -beta L_alpha Phi
            call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)

            if (abs(visc_type) .eq. 2) then
               ! multiply by c=2; x_p = -2*beta L_alpha Phi
               call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
            end if
 
          else if (abs(visc_type) .eq. 3) then

            ! multiply x_p by gamma, use mac_rhs a temparary to save x_p 
            call multifab_copy_c(mac_rhs(n),1,x_p(n),1,1,0)
            call multifab_mult_mult_c(mac_rhs(n),1,gamma(n),1,1,0)
            ! multiply x_p by beta; x_p = -beta L_alpha Phi
            call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
            ! multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
            call multifab_mult_mult_s_c(x_p(n),1,4.d0/3.d0,1,0)
            ! x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
            call multifab_plus_plus_c(x_p(n),1,mac_rhs(n),1,1,0)

          end if

          ! multiply Phi by theta
          call multifab_mult_mult_s_c(phi(n),1,theta,1,0)

          ! add theta*Phi to x_p
          call multifab_plus_plus_c(x_p(n),1,phi(n),1,1,0)

        end do

      case(2,5)  
        ! lower triangular, precon_type=-2 means using negative sign for Schur complement app

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve for x_u using an implicit viscous term
        !         A x_u = b_u
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
        if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
           print*,""
           print*,"Begin viscous solve for x_u"
        end if

        ! x_u = A^{-1} b_u
        call stag_mg_solver(mla,alpha,beta,gamma,theta,x_u,b_u,dx,the_bc_tower)

        if (parallel_IOProcessor().and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"End viscous solve for x_u"
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Solve a pressure Poisson problem for Phi
        !         L_alpha Phi = D x_u + b_p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! add set mac_rhs = D(x_u) to mac_rhs
        call compute_divu(mla,x_u,mac_rhs,dx)

        ! add b_p to mac_rhs
        do n=1,nlevs
          call multifab_plus_plus_c(mac_rhs(n),1,b_p(n),1,1,0)
        end do

        if (abs(theta) .gt. 0.d0) then 
           if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
             print*,""
             print*,"Begin projection"
           end if

           ! solves L_alpha Phi = mac_rhs
           ! x_u is only passed in to get a norm for absolute residual criteria
           call macproject(mla,phi,x_u,alpha,mac_rhs,dx,the_bc_tower)

           if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
             print*,""
             print*,"End projection"
           end if
        end if 
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Update x_p 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! x_p = theta*I *Phi - beta*mac_rhs 
        do n=1,nlevs
          ! beta part
          if ( (abs(visc_type) .eq. 1) .or. (abs(visc_type) .eq. 2) ) then 
             call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)
       
             ! multiply x_p by -1
             call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
      
             ! multiply x_p by beta; x_p = -beta L_alpha Phi
             call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
 
             if (abs(visc_type) .eq. 2) then
               ! multiply by c=2 for |viscous_type| = 2
               call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
             end if
 
          elseif (abs(visc_type) .eq. 3) then   
             ! beta part 
             call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)
        
             ! multiply x_p by -4/3
             call multifab_mult_mult_s_c(x_p(n),1,-4.d0/3.d0,1,0)
       
             ! multiply x_p by beta; x_p = -beta L_alpha Phi
             call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)

             ! gamma part
             ! x_p = theta*I *Phi - 4/3*beta*mac_rhs - gamma*mac_rhs 
             call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
             call multifab_mult_mult_c(mac_rhs(n),1,gamma(n),1,1,0)

             ! x_p = x_p + mac_rhs 
             call multifab_plus_plus_c(x_p(n),1,mac_rhs(n),1,1,0)

          end if

          if (abs(theta) .gt. 0.d0) then         
             ! multiply Phi by theta
             call multifab_mult_mult_s_c(phi(n),1,theta,1,0)   

             ! add theta*Phi to x_p
             call multifab_plus_plus_c(x_p(n),1,phi(n),1,1,0)
          end if   

          if (precon_type .eq. -2 .or. precon_type .eq. -5) then
            ! multiply x_p by -1, if precon_type=-2
            call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
          end if 
        end do

        if (abs(precon_type) .eq. 5) then

           ! compute = A^(-1)*(b_u-grad(x_p)) 

           ! we need gradients of x_p, and x_p doesn't necessarily have 
           ! a ghost cell so we create a temporary
           do n=1,nlevs
              call multifab_copy_c(x_p_tmp(n),1,x_p(n),1,1,0)
              call multifab_fill_boundary(x_p_tmp(n))
              call multifab_physbc(x_p_tmp(n),1,dm+1,1,the_bc_tower%bc_tower_array(n))
           end do

           ! contstruct b_u-grad(x_p)
           do n=1,nlevs
              do i=1,dm
                 call multifab_copy_c(b_u_tmp(n,i),1,b_u(n,i),1,1,0)
                 call setval(alphainv_edge(n,i),1.d0,all=.true.)
              end do
           end do
           call subtract_weighted_gradp(mla,b_u_tmp,alphainv_edge,x_p_tmp,dx)

           if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
              print*,""
              print*,"Begin viscous solve for x_u"
           end if

           ! compute = A^(-1)*(b_u-grad(x_p)) 
           call stag_mg_solver(mla,alpha,beta,gamma,theta,x_u,b_u_tmp,dx,the_bc_tower)

           if (parallel_IOProcessor().and. stag_mg_verbosity .ge. 1) then
              print*,""
              print*,"End viscous solve for x_u"
           end if

        end if
        
      case(3)  
        ! upper triangular, precon_type=-3 means using negative sign for Schur complement app

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve a pressure Poisson problem for Phi
        !         L_alpha Phi = b_p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do n=1,nlevs
          ! copy b_p to mac_rhs, mac_rhs is now the rhs for later use 
          call multifab_copy_c(mac_rhs(n),1,b_p(n),1,1,0)  
          if (precon_type .eq. -3) then         ! rhs times -1   
            call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
          end if 
        end do
  
        if (abs(theta) .gt. 0) then   
          if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
            print*,""
            print*,"Begin pressure Poisson solve"
          end if

          ! solves L_alpha Phi = mac_rhs
          ! x_u^star is only passed in to get a norm for absolute residual criteria
          call macproject(mla,phi,x_u,alpha,mac_rhs,dx,the_bc_tower)

          if (parallel_IOProcessor() .and. mg_verbose .ge. 1) then
            print*,""
            print*,"End pressure Poisson solve"
          end if
        end if
        
        do n=1,nlevs

          if ( (abs(visc_type) .eq. 1) .or. (abs(visc_type) .eq. 2) ) then
            ! x_p = beta*mac_rhs
            call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)      
            call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
          
            ! multiply by -c=-1 for |viscous_type| = 1
            call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
            if (abs(visc_type) .eq. 2) then
              ! multiply by -c=-2 for |viscous_type| = 2
              call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
            end if 

          elseif (abs(visc_type) .eq. 3) then
            ! x_p = 4/3*beta*mac_rhs
            call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)      
            call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
          
            ! multiply by -c=-4/3 for |viscous_type| = 3
            call multifab_mult_mult_s_c(x_p(n),1,-4.d0/3.d0,1,0)

            ! gamma part 
            call multifab_mult_mult_c(mac_rhs(n),1,gamma(n),1,1,0)
            call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
            ! x_p = x_p + gamma*mac_rhs
            call multifab_plus_plus_c(x_p(n),1,mac_rhs(n),1,1,0)
          end if
          
          if (abs(theta) .gt. 0) then 
             ! multiply phi by theta 
             call multifab_mult_mult_s_c(phi(n),1,theta,1,0)
             ! add phi to x_p                                           
             call multifab_plus_plus_c(x_p(n),1,phi(n),1,1,0)
          end if
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Compute the RHS for the viscous solve
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! we need gradients of x_p, and x_p doesn't necessarily have 
        ! a ghost cell so we create a temporary
        do n=1,nlevs
          call multifab_copy_c(x_p_tmp(n),1,x_p(n),1,1,0)
          call multifab_fill_boundary(x_p_tmp(n))
          call multifab_physbc(x_p_tmp(n),1,dm+1,1,the_bc_tower%bc_tower_array(n))
        end do

        ! contstruct b_u-grad(x_p)
        do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(b_u_tmp(n,i),1,b_u(n,i),1,1,0)
             call setval(alphainv_edge(n,i),1.d0,all=.true.)
          end do
        end do
        call subtract_weighted_gradp(mla,b_u_tmp,alphainv_edge,x_p_tmp,dx)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Solve for x_u using an implicit viscous term
        !         A x_u = (b_u-G x_p)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"Begin viscous solve for x_u"
        end if

        ! compute = A^(-1)*(b_u-grad(x_p)) 
        call stag_mg_solver(mla,alpha,beta,gamma,theta,x_u,b_u_tmp,dx,the_bc_tower)

        if (parallel_IOProcessor().and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"End viscous solve for x_u"
        end if

      case(4) 
        ! block diagonal

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve for x_u using an implicit viscous term
        !         A x_u = b_u
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (parallel_IOProcessor() .and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"Begin viscous solve for x_u"
        end if

        ! x_u = A^{-1} b_u
        call stag_mg_solver(mla,alpha,beta,gamma,theta,x_u,b_u,dx,the_bc_tower)

        if (parallel_IOProcessor().and. stag_mg_verbosity .ge. 1) then
          print*,""
          print*,"End viscous solve for x_u"
        end if

        ! x_p = -c*beta*b_p
        do n=1,nlevs
          ! x_p = beta*bp
          call multifab_copy_c(x_p(n),1,b_p(n),1,1,0)
          call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
          if ((abs(visc_type) .eq. 1) .and. (precon_type .eq. 4)) then
             ! multiply by -c=-1 for |viscous_type| = 1
             call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
          else if ((abs(visc_type) .eq. 2) .and. (precon_type .eq. 4)) then
             ! multiply by -c=-2 for |viscous_type| = 2
             call multifab_mult_mult_s_c(x_p(n),1,-2.d0,1,0)
          else if ((abs(visc_type) .eq. 2) .and. (precon_type .eq. -4)) then
             ! multiply by c=2 for precon_type = -4
             call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
          else if ((abs(visc_type) .eq. 3) .and. (precon_type .eq. 4)) then
             ! multiply by -c=-4/3 for |viscous_type| = 3
             call multifab_mult_mult_s_c(x_p(n),1,-4.d0/3.d0,1,0)
             ! gamma part: the sign is same as beta, use x_p_tmp as an immediate variable  
             call multifab_copy_c(x_p_tmp(n),1,b_p(n),1,1,0)
             call multifab_mult_mult_c(x_p_tmp(n),1,gamma(n),1,1,0)
             call multifab_mult_mult_s_c(x_p_tmp(n),1,-1.d0,1,0)
             call multifab_plus_plus_c(x_p(n),1,x_p_tmp(n),1,1,0)
             
          else if ((abs(visc_type) .eq. 3) .and. (precon_type .eq. -4)) then
             ! multiply by c=4/3 for precon_type = -4
             call multifab_mult_mult_s_c(x_p(n),1,4.d0/3.d0,1,0)
             ! gamma part: the sign is same as beta, use x_p_tmp as an immediate variable   
             call multifab_copy_c(x_p_tmp(n),1,b_p(n),1,1,0)
             call multifab_mult_mult_c(x_p_tmp(n),1,gamma(n),1,1,0)
             call multifab_plus_plus_c(x_p(n),1,x_p_tmp(n),1,1,0)             
           
          end if
        end do

     case default
        call bl_error('ERROR: Only precond_type 1 2 3 4 -2 -3 -4 supported')
     end select


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 4: Handle null-space issues in MG solvers
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! subtract off mean value: Single level only! No need for ghost cells
    call sum_umac_press(mla,x_p,x_u,mean_val_pres,mean_val_umac)
    
    ! The pressure Poisson problem is always singular:
    call multifab_sub_sub_s_c(x_p(1),1,mean_val_pres,1,0)
    
    ! The velocity problem is also singular for periodic systems with no identity piece
    if(all(mla%pmask(1:dm)) .and. (theta==0.0d0)) then
       do i=1,dm
          call multifab_sub_sub_s_c(x_u(1,i),1,mean_val_umac(i),1,0)
       end do
    end if 

    do n=1,nlevs
       call multifab_destroy(phi(n))
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(zero_fab(n))
       call multifab_destroy(x_p_tmp(n))
       do i=1,dm
          call multifab_destroy(alphainv_edge(n,i))
          call multifab_destroy(b_u_tmp(n,i))
       end do
    end do
    
  contains

    subroutine subtract_weighted_gradp(mla,x_u,alphainv_edge,phi,dx)

      type(ml_layout), intent(in   ) :: mla
      type(multifab ), intent(inout) :: x_u(:,:)
      type(multifab ), intent(in   ) :: alphainv_edge(:,:)
      type(multifab ), intent(in   ) :: phi(:)
      real(kind=dp_t), intent(in   ) :: dx(:,:)

      ! local
      integer :: i,dm,n,nlevs

      type(multifab) :: gradp(mla%nlevel,mla%dim)

      nlevs = mla%nlevel
      dm = mla%dim

      do n=1,nlevs
         do i=1,dm
            call multifab_build_edge(gradp(n,i),mla%la(n),1,1,i)
         end do
      end do

      call compute_gradp(mla,phi,gradp,dx)

      do n=1,nlevs
         do i=1,dm
            call multifab_mult_mult_c(gradp(n,i),1,alphainv_edge(n,i),1,1,0)
            call saxpy(x_u(n,i),-1.d0,gradp(n,i))
         end do
      end do

      do n=1,nlevs
         do i=1,dm
            call multifab_destroy(gradp(n,i))
         end do
      end do

    end subroutine subtract_weighted_gradp

  end subroutine apply_precon

end module apply_precon_module
