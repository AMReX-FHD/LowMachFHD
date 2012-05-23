PROGRAM TestLowMach_1
      USE Precision
      USE Random_Numbers
      USE Simple_Dense_Linear_Algebra
      IMPLICIT NONE

      INTEGER, PARAMETER :: shadow = 2
      INTEGER, DIMENSION (1) :: grid_shape, n_c, zeros = 0, ones = 1, n_solution_cells, cell_grid
      INTEGER, DIMENSION (1) :: r_shape
      INTEGER :: n_cells, dim, cell, face, step, i, j, k, n_steps, seed, n_vars, i_temp, use_fft, var, k_index, record, &
     & n_step_stats, n_stats, substep

      REAL (KIND=r_wp) :: r_temp, t, wall_speed, plug_speed, normalization, vel_scaling, pressure_scaling, vec_length, wall_dtemp
      REAL (KIND=r_wp) :: shear_visc, therm_cond, eq_density, eq_temperature, dt, cfl_factor, c_v, c_p, gamma, n_ppc, &
     & normalized_variance, dx, density, enthalpy, velocity, pressure, filter_width
      REAL (KIND=r_wp), DIMENSION (1) :: dr, position, r_temp_v
      REAL (KIND=r_wp) :: drho2, du2, dT2, drho2_av = 0, du2_av = 0, dT2_av = 0, dJ2, dJ2_av = 0, dE2, dE2_av = 0

      LOGICAL :: periodic (1), show_primary, render, dohydro = .TRUE.
      
      REAL (KIND=r_wp), DIMENSION (:, :), ALLOCATABLE :: conserved, fluxes, temp_arrays, facial_vals
      REAL (KIND=r_wp), DIMENSION (:), ALLOCATABLE :: densities, velocities, temperatures
      REAL (KIND=r_wp), DIMENSION (:), ALLOCATABLE :: adv_velocities, pressures, stochflux, stochfluxT
      INTEGER, DIMENSION (:), ALLOCATABLE :: grid_mask, rgrid_mask

      REAL (KIND=r_wp) :: fudge_factor = 0.1_r_wp
      LOGICAL :: diffadv = .FALSE., onestep = .TRUE., linearized = .FALSE., centered=.FALSE.
      
      periodic = .TRUE.

      WRITE (*,*) "Enter grid_shape, periodic, dr"
      READ (*,*) grid_shape, periodic (1:1), dr
      
      WRITE (*,*) "Enter eta, kappa, eq. density, temperature"
      READ (*,*) shear_visc, therm_cond, eq_density, eq_temperature
      
      WRITE (*,*) "Enter dt, n_steps, n_step_stats"
      READ (*,*) dt, n_steps, n_step_stats

      WRITE (*,*) "Enter variance, plug speed"
      READ (*,*) normalized_variance, plug_speed
      
      WRITE (*,*) "Enter linearized, diffadv, onestep, centered, fudge factor"
      READ (*,*) linearized, diffadv, onestep, centered, fudge_factor
      
      IF (diffadv) onestep = .TRUE.

      CALL UnpredictableSeeds()
      
      gamma = (2.0_r_wp+1) / 1
      c_v = 1.0_r_wp / (gamma-1.0_r_wp)
      c_p = gamma * c_v
      dx = dr (1)
      n_ppc = eq_density * PRODUCT (dr)
      
      WRITE (*,*) "CFLs: alpha=", Abs (plug_speed) * dt / dx, " beta=", dt * shear_visc / eq_density / dx ** 2, " gamma=", dt * &
     & therm_cond / eq_density / c_p / dx ** 2, " rs=", Abs (plug_speed) * dx * eq_density / shear_visc, Abs (plug_speed) * dx * &
     & c_p * eq_density / therm_cond
      
      n_vars = 1 + 2
      n_cells = PRODUCT (grid_shape)
      
      ALLOCATE (conserved(-1:n_vars, n_cells))
      ALLOCATE (fluxes(n_vars, 0:n_cells))
      ALLOCATE (densities(1-shadow:n_cells+shadow), velocities(1-shadow:n_cells+shadow), adv_velocities(1-shadow:n_cells+shadow), &
     & temperatures(1-shadow:n_cells+shadow), pressures(1-shadow:n_cells+shadow))
      ALLOCATE (temp_arrays(1-shadow:n_cells+shadow, 0:6), facial_vals(1:n_vars, 1-shadow:n_cells+shadow))
      ALLOCATE (stochflux(0:n_cells+1), stochfluxT(0:n_cells+1))

      ! Initialize
      CALL RandomNormal (conserved)
      r_temp = 0 ! Magnitude of initial perturbation
      conserved (1, 1:n_cells) = eq_density * (1+r_temp*conserved(1, 1:n_cells))
      conserved (2, 1:n_cells) = plug_speed * eq_density
      conserved (3, 1:n_cells) = eq_density * c_p * eq_temperature * (1+r_temp*conserved(3, 1:n_cells))
      CALL ConservedToPrimary ()
      CALL FixBCs ()

      ! Start iteration
      pressures = 0.0_r_wp
      t = 0
      n_stats = 0
      DO step = 1, n_steps
         CALL StochasticFluxes ()
         DO substep = 1, MERGE (1, 2, onestep)
            fluxes (1, 1:n_cells) = 0.0
            fluxes (2, 1:n_cells) = - stochflux (1:n_cells)
            fluxes (3, 1:n_cells) = - stochfluxT (1:n_cells)
            
            ! Initial guess for advective velocities
            
            ! Deviation from equation of state:
            IF (linearized) THEN
               DO cell = 1, n_cells
                  temp_arrays (cell, 0) = fudge_factor / gamma / dt * &
                 & ((densities(cell)-eq_density)*eq_temperature+eq_density*(temperatures(cell)-eq_temperature)) / eq_density / &
                 & eq_temperature
               END DO
            ELSE
               DO cell = 1, n_cells
                  temp_arrays (cell, 0) = fudge_factor * (densities(cell)*temperatures(cell)/eq_density/eq_temperature-1.0_r_wp) &
                 & / gamma / dt
               END DO
            END IF
            
            ! Advective velocities:
            IF (diffadv) THEN ! Keep them constant
               adv_velocities (1:n_cells) = plug_speed
            ELSE IF (onestep) THEN
               r_temp = 0.0_r_wp
               DO face = 1, n_cells
                  r_temp = r_temp + (dx*temp_arrays(face, 0))
                  adv_velocities (face) = plug_speed + r_temp + (gamma-1) / gamma / eq_density / eq_temperature * &
                 & (therm_cond*(temperatures(face+1)-temperatures(face))/dx+stochfluxT(face))
               END DO
            ELSE
               r_temp = 0.0_r_wp
               DO face = 1, n_cells
                  r_temp = r_temp + (dx*temp_arrays(face, 0))
                  IF (substep == 1) THEN
                     adv_velocities (face) = plug_speed + r_temp + (gamma-1) / gamma / eq_density / eq_temperature * &
                    & (therm_cond*(temperatures(face+1)-temperatures(face))/dx+stochfluxT(face))
                  ELSE
                     adv_velocities (face) = plug_speed + r_temp + (gamma-1) / gamma / eq_density / eq_temperature * &
                    & (0.5_r_wp*therm_cond*(temp_arrays(face+1, 4)+temperatures(face+1)-temp_arrays(face, &
                    & 4)-temperatures(face))/dx+stochfluxT(face))
                  END IF
               END DO
               IF (.FALSE.) WRITE (*,*) "TEST=", SUM (temp_arrays(1:n_cells, 0)), SUM (adv_velocities(1:n_cells)) / n_cells - &
              & plug_speed
            END IF
            
            DO cell = 1, shadow
               adv_velocities (n_cells+cell) = adv_velocities (cell)
            END DO
            DO cell = 1 - shadow, 0
               adv_velocities (cell) = adv_velocities (n_cells+cell)
            END DO
            
            ! Hyperbolic flux estimates
            IF (centered) THEN ! No upwinding
               CALL FacialInterpolation (u=densities, u_face=facial_vals(1, :))
               CALL FacialInterpolation (u=velocities, u_face=facial_vals(2, :))
               CALL FacialInterpolation (u=temperatures, u_face=facial_vals(3, :))
            ELSE ! Upwinding
               CALL SpatialSlopes (u=densities, du_dx=temp_arrays(1-shadow:n_cells+shadow, 1))
               CALL SpatialSlopes (u=velocities, du_dx=temp_arrays(1-shadow:n_cells+shadow, 2))
               CALL SpatialSlopes (u=temperatures, du_dx=temp_arrays(1-shadow:n_cells+shadow, 3))
               DO face = 1, n_cells
                  IF (linearized) THEN
                     IF (plug_speed >= 0.0_r_wp) THEN
                        r_temp = 0.5_r_wp * (dx-dt*plug_speed)
                        cell = face
                     ELSE
                        r_temp = - 0.5_r_wp * (dx+dt*plug_speed)
                        cell = face + 1
                     END IF
                     facial_vals (1, face) = densities (cell) + r_temp * temp_arrays (cell, 1) - 0.5_r_wp * dt * eq_density * &
                    & (adv_velocities(cell)-adv_velocities(cell-1)) / dx
                     facial_vals (2, face) = velocities (cell) + r_temp * temp_arrays (cell, 2) + 0.5_r_wp * dt / eq_density * (-&
                    & (pressures(cell)-pressures(cell-1))+shear_visc*(velocities(cell-1)-2*velocities(cell)+velocities(cell+&
                    & 1))/dx**2+(stochflux(cell)-stochflux(cell-1))/dx)
                     facial_vals (3, face) = temperatures (cell) + r_temp * temp_arrays (cell, 3) + 0.5_r_wp * dt / eq_density / &
                    & c_p * (therm_cond*(temperatures(cell-1)-2*temperatures(cell)+temperatures(cell+1))/dx**2+(stochfluxT(cell)-&
                    & stochfluxT(cell-1))/dx)
                  ELSE
                     IF (adv_velocities(face) >= 0.0_r_wp) THEN
                        r_temp = 0.5_r_wp * (dx-dt*adv_velocities(face))
                        cell = face
                     ELSE
                        r_temp = - 0.5_r_wp * (dx+dt*adv_velocities(face))
                        cell = face + 1
                     END IF
                     facial_vals (1, face) = densities (cell) + r_temp * temp_arrays (cell, 1) - 0.5_r_wp * dt * densities (cell) &
                    & * (adv_velocities(cell)-adv_velocities(cell-1)) / dx
                     facial_vals (2, face) = velocities (cell) + r_temp * temp_arrays (cell, 2) + 0.5_r_wp * dt / densities &
                    & (cell) * (-(pressures(cell)-pressures(cell-1))+shear_visc*(velocities(cell-1)-2*velocities(cell)+&
                    & velocities(cell+1))/dx**2+(stochflux(cell)-stochflux(cell-1))/dx)
                     facial_vals (3, face) = temperatures (cell) + r_temp * temp_arrays (cell, 3) + 0.5_r_wp * dt / densities &
                    & (cell) / c_p * (therm_cond*(temperatures(cell-1)-2*temperatures(cell)+temperatures(cell+1))/dx**2+&
                    & (stochfluxT(cell)-stochfluxT(cell-1))/dx)
                  END IF
               END DO
            END IF
            
            IF (linearized) THEN
               DO face = 1, n_cells
                  fluxes (1, face) = fluxes (1, face) + (facial_vals(1, face)*plug_speed+adv_velocities(face)*eq_density)
                  fluxes (2, face) = fluxes (2, face) + (eq_density*facial_vals(2, face)*plug_speed+pressures(face))
                  fluxes (3, face) = fluxes (3, face) + (c_p*eq_density*facial_vals(3, face)*plug_speed)
               END DO
            ELSE
               DO face = 1, n_cells
                  fluxes (1, face) = fluxes (1, face) + (facial_vals(1, face)*adv_velocities(face))
                  fluxes (2, face) = fluxes (2, face) + (facial_vals(1, face)*facial_vals(2, &
                 & face)*adv_velocities(face)+pressures(face))
                  fluxes (3, face) = fluxes (3, face) + (c_p*facial_vals(1, face)*facial_vals(3, face)*adv_velocities(face))
               END DO
            END IF
            
            ! Diffusive explicit fluxes
            DO cell = 1, n_cells
               fluxes (2, cell) = fluxes (2, cell) - (0.5_r_wp*shear_visc*(velocities(cell+1)-velocities(cell))/dx)
               fluxes (3, cell) = fluxes (3, cell) - (0.5_r_wp*therm_cond*(temperatures(cell+1)-temperatures(cell))/dx)
            END DO
            
            ! BC fix
            fluxes (:, 0) = fluxes (:, n_cells)
            
            ! New values
            PCStage: IF ((substep == 1) .AND. ( .NOT. onestep)) THEN
            
               ! Predictor stage has only temperature implicit solve
               DO cell = 1, n_cells
                  IF (linearized) THEN
                     density = eq_density
                  ELSE
                     density = conserved (1, cell) - dt * (fluxes(1, cell)-fluxes(1, cell-1)) / dx
                  END IF
                  enthalpy = conserved (3, cell) - dt * (fluxes(3, cell)-fluxes(3, cell-1)) / dx
                  temp_arrays (cell, 1) = enthalpy
                  temp_arrays (cell, 2) = density * c_p + dt * therm_cond / dx ** 2
                  temp_arrays (cell, 3) = - 0.5_r_wp * dt * therm_cond / dx ** 2
               END DO
               CALL PeriodicTridiagonalSolve (temp_arrays(1:n_cells, 3), temp_arrays(1:n_cells, 2), temp_arrays(1:n_cells, 3), &
              & temp_arrays(1:n_cells, 1), temp_arrays(1:n_cells, 4), n_cells)
               CALL FixTempBCs (4)
               
            ELSE PCStage
            
               DO cell = 1, n_cells
                  conserved (1:n_vars, cell) = conserved (1:n_vars, cell) - (dt*(fluxes(1:n_vars, cell)-fluxes(1:n_vars, &
                 & cell-1))/dx)
               END DO
               
               ! Velocity implicit solve
               DO cell = 1, n_cells
                  IF (linearized) THEN
                     density = eq_density
                  ELSE
                     density = conserved (1, cell)
                  END IF
                  temp_arrays (cell, 1) = conserved (3, cell)
                  temp_arrays (cell, 2) = density * c_p + dt * therm_cond / dx ** 2
                  temp_arrays (cell, 3) = - 0.5_r_wp * dt * therm_cond / dx ** 2
               END DO               
               CALL PeriodicTridiagonalSolve (temp_arrays(1:n_cells, 3), temp_arrays(1:n_cells, 2), temp_arrays(1:n_cells, 3), &
              & temp_arrays(1:n_cells, 1), temp_arrays(1:n_cells, 4), n_cells)
               temperatures (1:n_cells) = temp_arrays (1:n_cells, 4)
                              
               ! Temperature implicit solve
               DO cell = 1, n_cells
                  IF (linearized) THEN
                     density = eq_density
                  ELSE
                     density = conserved (1, cell)
                  END IF
                  temp_arrays (cell, 1) = conserved (2, cell)
                  temp_arrays (cell, 2) = density + dt * shear_visc / dx ** 2
                  temp_arrays (cell, 3) = - 0.5_r_wp * dt * shear_visc / dx ** 2
               END DO               
               CALL PeriodicTridiagonalSolve (temp_arrays(1:n_cells, 3), temp_arrays(1:n_cells, 2), temp_arrays(1:n_cells, 3), &
              & temp_arrays(1:n_cells, 1), temp_arrays(1:n_cells, 4), n_cells)
               velocities (1:n_cells) = temp_arrays (1:n_cells, 4)
               
               ! Density for the projection step:
               IF (.FALSE.) THEN ! This does not work because it violates conservation
                  densities (1:n_cells) = 0.5_r_wp * (densities(1:n_cells)+conserved(1, 1:n_cells))
               ELSE
                  densities (1:n_cells) = conserved (1, 1:n_cells)
               END IF
               
               CALL FixBCs ()
               
               ! Velocity projection
               DiffusionAdvection: IF ( .NOT. diffadv) THEN
               
                  ! Calculate the target velocity divergence
                  IF (linearized) THEN
                     DO cell = 1, n_cells
                        temp_arrays (cell, 1) = (gamma-1) / gamma / eq_density / eq_temperature * (therm_cond*(temperatures(cell-&
                       & 1)-2*temperatures(cell)+temperatures(cell+1))/dx**2+(stochfluxT(cell)-stochfluxT(cell-1))/dx) + &
                       & fudge_factor / gamma / dt * ((densities(cell)-eq_density)*eq_temperature+eq_density*(temperatures(cell)-&
                       & eq_temperature)) / eq_density / eq_temperature
                     END DO
                  ELSE
                     DO cell = 1, n_cells
                        temp_arrays (cell, 1) = (gamma-1) / gamma / eq_density / eq_temperature * (therm_cond*(temperatures(cell-&
                       & 1)-2*temperatures(cell)+temperatures(cell+1))/dx**2+(stochfluxT(cell)-stochfluxT(cell-1))/dx) + &
                       & fudge_factor * (densities(cell)*temperatures(cell)/eq_density/eq_temperature-1.0_r_wp) / gamma / dt
                     END DO
                  END IF
                  CALL FixTempBCs (1)
                  
                  ! We have to project it onto faces
                  DO face = 1, n_cells
                     temp_arrays (face, 2) = (temp_arrays(face, 1)+temp_arrays(face+1, 1)) / 2
                  END DO
                  
                  ! Pressure solve:
                  DO cell = 1, n_cells
                     temp_arrays (cell, 1) = ((velocities(cell+1)-velocities(cell))/dx-temp_arrays(cell, 2)) / dt
                     IF (linearized) THEN
                        temp_arrays (cell, 3) = 1.0_r_wp / eq_density / dx ** 2
                        temp_arrays (cell, 4) = - 2.0_r_wp / eq_density / dx ** 2
                        temp_arrays (cell, 5) = 1.0_r_wp / eq_density / dx ** 2
                     ELSE
                        temp_arrays (cell, 3) = 1.0_r_wp / densities (cell) / dx ** 2
                        temp_arrays (cell, 4) = - (1.0_r_wp/densities(cell+1)+1.0_r_wp/densities(cell)) / dx ** 2
                        temp_arrays (cell, 5) = 1.0_r_wp / densities (cell+1) / dx ** 2
                     END IF
                  END DO
                  CALL PeriodicTridiagonalSolve (a=temp_arrays(1:n_cells, 3), b=temp_arrays(1:n_cells, 4), &
                 & c=temp_arrays(1:n_cells, 5), r=temp_arrays(1:n_cells, 1), u=temp_arrays(1:n_cells, 6), n=n_cells)
                  CALL FixTempBCs ()
                  
                  ! Project:
                  DO cell = 1, n_cells
                     IF (linearized) THEN
                        density = eq_density
                     ELSE
                        density = densities (cell)
                     END IF
                     velocities (cell) = velocities (cell) - (dt*(temp_arrays(cell, 6)-temp_arrays(cell-1, 6))/dx/density)
                     pressures (cell) = pressures (cell) + (temp_arrays(cell, 6))
                  END DO
                  
                  IF (.TRUE.) THEN ! Test linear solver
                     DO cell = 1, n_cells
                        IF (linearized) THEN
                           density = eq_density
                        ELSE
                           density = densities (cell)
                        END IF
                        temp_arrays (cell, 3) = (temp_arrays(cell, 6)-temp_arrays(cell-1, 6)) / dx / density
                     END DO
                     CALL FixTempBCs (3)
                     DO cell = 1, n_cells
                        temp_arrays (cell, 4) = (temp_arrays(cell+1, 3)-temp_arrays(cell, 3)) / dx
                     END DO
                     DO face = 1, n_cells
                        r_temp = temp_arrays (face, 4) - temp_arrays (face, 1)
                        IF (Abs(r_temp) > 1.0E-3_r_wp) THEN
                           WRITE (*,*) face, " solver blew up remainder=", r_temp
                           STOP
                        ELSE IF (Abs(r_temp) > 1.0E-9_r_wp) THEN
                           WRITE (*,*) face, " mismatch=", r_temp
                        END IF
                     END DO
                  END IF
                  
                  IF (.TRUE.) THEN ! Test velocity divergence
                     velocities (n_cells+1) = velocities (1)
                     DO face = 1, n_cells
                        cell = face
                        r_temp = (velocities(cell+1)-velocities(cell)) / dx - temp_arrays (face, 2)
                        IF (Abs(r_temp) > 1.0E-3_r_wp) THEN
                           WRITE (*,*) face, " solver blew up div_mismatch=", r_temp
                        ELSE IF (Abs(r_temp) > 1.0E-9_r_wp) THEN
                           WRITE (*,*) face, " div error=", r_temp
                        END IF
                     END DO
                  END IF
                  
                  densities (1:n_cells) = conserved (1, 1:n_cells)
                  
               END IF DiffusionAdvection
            END IF PCStage
         END DO
         
         CALL PrimaryToConserved ()
         CALL FixBCs ()
         
         IF ((step/n_step_stats)*n_step_stats == step) THEN
            n_stats = n_stats + (1)
            IF (.TRUE.) THEN ! Verify conservation
               IF (.TRUE.) WRITE (*,*) "E/M/P/H=", REAL (SUM(conserved(0:n_vars, 1:n_cells), dim=2))
               IF (.FALSE.) WRITE (*,*) "<U^2>=", REAL (SUM((densities(1:n_cells)/eq_density-1)**2)) / n_cells, REAL &
              & (SUM((velocities(1:n_cells)-plug_speed)**2)) / n_cells, REAL (SUM((temperatures(1:n_cells)/eq_temperature-1)**2)) &
              & / n_cells
            END IF
         END IF
         IF (step >= n_steps) EXIT
         t = t + dt
      END DO
      
      WRITE (*,*) "Took ", step, " timesteps successfully"
      n_steps = step

CONTAINS
      SUBROUTINE ConservedToPrimary ()
         DO cell = 1, n_cells
            densities (cell) = conserved (1, cell)
            IF (linearized) THEN
               density = eq_density
            ELSE
               density = densities (cell)
            END IF
            velocities (cell) = conserved (2, cell) / density
            temperatures (cell) = conserved (3, cell) / density / c_p
         END DO
      END SUBROUTINE
      SUBROUTINE PrimaryToConserved ()
         DO cell = 1, n_cells
            conserved (1, cell) = densities (cell)
            IF (linearized) THEN
               density = eq_density
            ELSE
               density = densities (cell)
            END IF
            conserved (2, cell) = velocities (cell) * density
            conserved (3, cell) = c_p * temperatures (cell) * density
         END DO
      END SUBROUTINE
      SUBROUTINE FixBCs ()
         DO cell = 1, shadow
            densities (n_cells+cell) = densities (cell)
            velocities (n_cells+cell) = velocities (cell)
            temperatures (n_cells+cell) = temperatures (cell)
            pressures (n_cells+cell) = pressures (cell)
         END DO
         DO cell = 1 - shadow, 0
            densities (cell) = densities (n_cells+cell)
            velocities (cell) = velocities (n_cells+cell)
            temperatures (cell) = temperatures (n_cells+cell)
            pressures (cell) = pressures (n_cells+cell)
         END DO
      END SUBROUTINE
      SUBROUTINE FixTempBCs (var)
         INTEGER, INTENT (IN), OPTIONAL :: var
         IF (PRESENT(var)) THEN
            DO cell = 1, shadow
               temp_arrays (n_cells+cell, var) = temp_arrays (cell, var)
            END DO
            DO cell = 1 - shadow, 0
               temp_arrays (cell, var) = temp_arrays (n_cells+cell, var)
            END DO
         ELSE
            DO cell = 1, shadow
               temp_arrays (n_cells+cell, :) = temp_arrays (cell, :)
            END DO
            DO cell = 1 - shadow, 0
               temp_arrays (cell, :) = temp_arrays (n_cells+cell, :)
            END DO
         END IF
      END SUBROUTINE
      SUBROUTINE SpatialSlopes (u, du_dx)
         REAL (KIND=r_wp), DIMENSION (1-shadow:), INTENT (IN) :: u
         REAL (KIND=r_wp), DIMENSION (1-shadow:), INTENT (OUT) :: du_dx
         DO cell = 1, n_cells
            IF (.TRUE.) THEN
               du_dx (cell) = 2 * (u(cell+1)-u(cell-1)) / (3*dx) - (u(cell+2)-u(cell-2)) / (12*dx)
            ELSE
               du_dx (cell) = 3 * (u(cell+1)-u(cell-1)) / (4*dx) - (u(cell+2)-u(cell-2)) / (8*dx)
            END IF
         END DO
         DO cell = 1, shadow
            du_dx (n_cells+cell) = du_dx (cell)
         END DO
         DO cell = 1 - shadow, 0
            du_dx (cell) = du_dx (n_cells+cell)
         END DO
      END SUBROUTINE
      SUBROUTINE FacialInterpolation (u, u_face)
         REAL (KIND=r_wp), DIMENSION (1-shadow:), INTENT (IN) :: u
         REAL (KIND=r_wp), DIMENSION (1-shadow:), INTENT (OUT) :: u_face
         IF (.FALSE.) THEN
            DO cell = 1, n_cells
               u_face (cell) = (u(cell+1)+u(cell)) / 2
            END DO
         ELSE
            DO cell = 1, n_cells
               u_face (cell) = (7*(u(cell+1)+u(cell))-(u(cell+2)+u(cell-1))) / 12
            END DO
         END IF
         DO cell = 1, shadow
            u_face (n_cells+cell) = u_face (cell)
         END DO
         DO cell = 1 - shadow, 0
            u_face (cell) = u_face (n_cells+cell)
         END DO
      END SUBROUTINE
      SUBROUTINE StochasticFluxes ()
         CALL RandomNormal (stochflux(1:n_cells))
         CALL RandomNormal (stochfluxT(1:n_cells))
         stochflux (1:n_cells) = stochflux (1:n_cells) * (Sqrt(normalized_variance*2*shear_visc*eq_temperature/(dx*dt)))
         stochflux (0) = stochflux (n_cells)
         stochflux (n_cells+1) = stochflux (1)
         stochfluxT (1:n_cells) = stochfluxT (1:n_cells) * (Sqrt(normalized_variance*2*therm_cond*eq_temperature**2/(dx*dt)))
         stochfluxT (0) = stochfluxT (n_cells)
         stochfluxT (n_cells+1) = stochfluxT (1)
      END SUBROUTINE
END PROGRAM
