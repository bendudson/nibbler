;
; Nibbler - Drift-Kinetic 3D Monte-Carlo code
;
; Version history:
;
; Nov 2010 :  Ben Dudson <bd512@york.ac.uk>
;     * Initial version
;
; Evolve position in toroidal coordinates (R,Z,phi)
; and parallel velocity of particles.
;

; Cross-product of two vectors
FUNCTION cross, a, b
  RETURN, {r:(a.phi*b.z - b.phi*a.z), $
           z:(a.r*b.phi - b.r*a.phi), $
           phi:(a.z*b.r - b.z*a.r)}
END

FUNCTION dot, a, b
  RETURN, a.r*b.r + a.z*b.z + a.phi*b.phi
END

FUNCTION mul, vec, num
  RETURN, {r:(vec.r*num), z:(vec.z*num), phi:(vec.phi*num)}
END

FUNCTION add, a, b
  RETURN, {r:(a.r+b.r), z:(a.z+b.z), phi:(a.phi+b.phi)}
END

;Calculate equilibrium B field
FUNCTION getEqBfield, ri, zi
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d 
  
  r = INTERPOLATE(r1d, ri)
  drdi = INTERPOLATE(DERIV(r1d), ri)
  dzdi = INTERPOLATE(DERIV(z1d), zi)

  g = local_gradient(dctfpol, ri, zi)
  Bphi = g.f / r
  dBphidr = (g.dfdr/drdi)/r - Bphi/r
  dBphidz = (g.dfdz/dzdi)/r
  
  ;; For poloidal field derivatives, need second
  ;; derivatives of psi
     
  g = EvalCosP(dctpsi, x0=ri, y0=zi)
  ;; g contains [ psi, dpsidr, dpsidz, d2psidr2, d2pdidz2,
  ;; d2psidrdz ]
     
  Br = -g[2] / dzdi / r   ; Br = -(dpsi/dz) / R
  Bz = g[1] / drdi / r    ; Bz = (dpsi/dR) / R
  
  dBrdr = -( (g[5]/(drdi * dzdi)) + Br ) / r
  dBrdz = - g[4] / (dzdi^2) / r
     
  dBzdr = ( (g[3]/(drdi^2)) - Bz ) / r
  dBzdz = g[5] / (drdi * dzdi) / r
  
  RETURN, {Br:Br, Bz:Bz, Bphi:Bphi, $
           dBrdr:dBrdr, dBrdz:dBrdz, $
           dBzdr:dBzdr, dBzdz:dBzdz, $
           dBphidr:dBphidr, dBphidz:dBphidz}
END

PRO oplot_bvec, ri, zi
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d 
  
  r = INTERPOLATE(r1d, ri)
  drdi = INTERPOLATE(DERIV(r1d), ri)
  dzdi = INTERPOLATE(DERIV(z1d), zi)
  
  g = getEqBfield(ri, zi)
  
  Br = g.Br
  Bz = g.Bz
  
  B = SQRT(Br^2 + Bz^2)
  
  dr = Br / B
  dz = Bz / B
  
  r0 = INTERPOLATE(r1d, ri)
  z0 = INTERPOLATE(z1d, zi)
  
  OPLOT, [r0, r0+dr], [z0, z0+dz], thick=2, color=2
  
END

FUNCTION sign, x
  IF x GT 0. THEN RETURN, 1.
  RETURN, -1.
END

;; Calculate time derivatives
FUNCTION differential, t, y
  COMMON particle_com, N, mu, mass, charge
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d 
  
  ;; Extract position and velocity from y
  ;; ri and zi are indices into r and z arrays
  ri = y[0:(N-1)]
  zi = y[N:(2*N-1)]
  phi = y[(2*N):(3*N-1)]
  vpar = y[(3*N):(4*N-1)]
  
  ;; Get the equilibrium magnetic field at these locations
  
  B    = DBLARR(N)
  
  Br    = B
  Bz    = B
  Bphi  = B
  
  dBrdr   = B
  dBrdz   = B
  dBrdphi = B

  dBzdr   = B
  dBzdz   = B
  dBzdphi = B
  
  dBphidr   = B
  dBphidz   = B
  dBphidphi = B
  
  r = INTERPOLATE(r1d, ri)
  drdi = INTERPOLATE(DERIV(r1d), ri)
  dzdi = INTERPOLATE(DERIV(z1d), zi)
  
  FOR i=0, N-1 DO BEGIN
     eqB = getEqBfield(ri[i], zi[i])
     
     Br[i] = eqB.Br
     Bz[i] = eqB.Bz
     Bphi[i] = eqB.Bphi
     
     dBrdr[i] = eqB.dBrdr
     dBrdz[i] = eqB.dBrdz
     
     dBzdr[i] = eqB.dBzdr
     dBzdz[i] = eqB.dBzdz
     
     dBphidr[i] = eqB.dBphidr
     dBphidz[i] = eqB.dBphidz
  ENDFOR
  
  ;; Add the perturbed B from coils
  
  
  ;; Get magnitude of B 
  B = SQRT(Br^2 + Bz^2 + Bphi^2)
  
  ;; Gradient of |B| (for Grad-B drift)
  gradB = {r:(( Br*dBrdr + Bz*dBzdr + Bphi*dBphidr ) / B), $
           z:(( Br*dBrdz + Bz*dBzdz + Bphi*dBphidz ) / B), $
           phi:(( Br*dBrdphi + Bz*dBzdphi + Bphi*dBphidphi ) / B)}

  ;; Gradients of little b vector (for curvature drift)

  ;lbrdr   = (dBrdr - dBdr*Br/B)/B
  lbrdz   = (dBrdz - gradB.z*Br/B)/B
  lbrdphi = (dBrdphi - gradB.phi*Br/B)/B
  
  lbzdr   = (dBzdr - gradB.r*Bz/B)/B
  ;lbzdz   = (dBzdz - dBdz*Bz/B)/B
  lbzdphi = (dBzdphi - gradB.phi*Bz/B)/B
  
  lbphidr   = (dBphidr - gradB.r*Bphi/B)/B
  lbphidz   = (dBphidz - gradB.z*Bphi/B)/B
  ;lbphidphi = (dBphidphi - dBdphi*Bphi/B)/B
  
  bvec = {r:(Br/B), z:(Bz/B), phi:(Bphi/B)}
  curlb = {r:(lbzdphi/r - lbphidz), phi:(lbrdz - lbzdr), $
           z:((bvec.phi+r*lbphidr - lbrdphi)/r)}
  
  kappa = cross( curlb, bvec )  ; (b dot grad) b

  ;; Calculate perpendicular drift velocity
  invwcj = mass / (charge * B)  ; 1 / wcj
  
  vd = mul( cross( bvec, add( mul(kappa, vpar^2), mul(gradB, mu) )), invwcj )
  ;vd = {r:0.0, z:0.0, phi:0.0}
  
  ;; Add parallel velocity
  v = add( vd, mul(bvec, vpar) )
  
  ; Calculate parallel acceleration (mirror force)
  dvpardt = -mu * dot( bvec, gradB ) / mass
  
  ;cursor, a, b, /down

  return, [v.r / drdi, v.z / dzdi, v.phi / r, dvpardt]
END

PRO nibbler, Nparticles=Nparticles, shot=shot, electron=electron, temp=temp, psin=psin, $
  kpar=kpar, output=output
  COMMON particle_com, N, mu, mass, charge
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d

  IF NOT KEYWORD_SET(Nparticles) THEN Nparticles = 1 ; Number of particles
  IF NOT KEYWORD_SET(temp) THEN temp = 200 ; Temperature in eV
  IF NOT KEYWORD_SET(psin) THEN psin = 0.9 ; Starting psi location
  IF NOT KEYWORD_SET(kpar) THEN kpar = 0.3 ; Vpar^2 / V^2
  
  N = Nparticles
  
  IF KEYWORD_SET(shot) THEN BEGIN
     ;; Fetch MAST equilibrium from IDAM
     
     PRINT, "Sorry, code not written yet."
     RETURN
  ENDIF ELSE BEGIN
     ;; Read in neqdsk file to get Psi
     file = "g014220.00200"
     
     PRINT, "Reading G-EQDSK file"
     grid = read_neqdsk(file)
     
     nr = grid.nx               ; Number of points in R
     nz = grid.ny               ; Number of points in Z
     psi = grid.psi             ; Poloidal flux on [nr, nz] mesh
     
     fpol1d = grid.fpol         ; f on uniform psi grid
     
     r2d = grid.r
     z2d = grid.z
     
     r1d = REFORM(r2d[*,0])
     z1d = REFORM(z2d[0,*])
  ENDELSE
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; By this point, need
  ;; nr, nz
  ;; psi[nr,nz], fpol1d[psi], r1d[nr], z1d[nz]
  
  ;; Plot contour
  nlev = 100
  minf = MIN(psi)
  maxf = MAX(psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  safe_colors, /first
  CONTOUR, psi, r2d, z2d, levels=levels, color=1, /iso, xstyl=1, ysty=1, $
    yr=[-0.7,0.7], xr=[0.6,1.2]

  ;; Analyse the equilibrium
  aeq = analyse_equil( psi, r1d, z1d )
  
  ;; Overplot the separatrices, O-points
  oplot_critical, psi, r1d, z1d, aeq
  
  ;; Categorise points inside and outside core
  core = categorise(psi, aeq=aeq, psinorm=psinorm) ; 1 where in core, 0 otherwise
  
  ;; Interpolate f onto grid points
  npgrid = (FINDGEN(N_ELEMENTS(fpol1d))/(N_ELEMENTS(fpol1d)-1))
  fpol2d = DBLARR(nr, nz)
  FOR i=0, nr-1 DO BEGIN
    fpol2d[i,*] = INTERPOL(fpol1d, npgrid, psinorm[i,*], /spline)
  ENDFOR

  PRINT, "Calculating DCT of Psi..."
  DCT2Dslow, psi, dctpsi
  
  PRINT, "Calculating DCT of fpol..."
  DCT2Dslow, fpol2d, dctfpol
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Get starting point
  contour_lines, psinorm, findgen(nr), findgen(nz), levels=[psin], $
                 path_info=info, path_xy=xy
  ;; Get the primary O-point
  opt_ri = aeq.opt_ri[aeq.primary_opt]
  opt_zi = aeq.opt_zi[aeq.primary_opt]
  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ;; Find the surface closest to the o-point
      
     ind = closest_line(info, xy, opt_ri, opt_zi)
     info = info[ind]
  ENDIF ELSE info = info[0]
  rinds = xy[0,info.offset:(info.offset+info.n-1)]
  zinds = xy[1,info.offset:(info.offset+info.n-1)]
  
  OPLOT, INTERPOLATE(r1d, rinds), INTERPOLATE(z1d, zinds), color=4, thick=2
  
  ; Get starting position at outboard midplane
  ri0 = MAX(rinds, pos)
  zi0 = zinds[pos]
  
  ; Get B here. Should be minimum B
  eqB = getEqBfield(ri0, zi0)
  Bmin = SQRT(eqB.Br^2 + eqB.Bz^2 + eqB.Bphi^2)
  

  ; Get maximum B (to work out trapped/passing boundary)
  ri1 = MIN(rinds, pos)
  zi1 = zinds[pos]
  eqB = getEqBfield(ri1, zi1)
  Bmax = SQRT(eqB.Br^2 + eqB.Bz^2 + eqB.Bphi^2)
  
  trapbndry = 1 - Bmin / Bmax
  
  ;; Check if trapped
  PRINT, "Particle trapped if |kpar| = vpar^2/v^2 < ", trapbndry

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Generate initial spread of kinetic energies and magnetic moment
  
  ; Everything in SI units?
  
  ke = FLTARR(N)
  mu = FLTARR(N)
  
  IF KEYWORD_SET(electron) THEN BEGIN
     mass = 9.11e-31
     charge = -1.602e-19
  ENDIF ELSE BEGIN
     ; Deuterium ion
     mass = 2.*1.67e-27
     charge = 1.602e-19
  ENDELSE
  
  ;; Choose distribution of energy
  
  ke = ke + 1.5*1.602e-19*temp  ; Energy in J
  
  ;; Choose distribution of parallel velocity
  kepar = ABS(kpar)*ke ; Kinetic energy in parallel direction
  
  vpar = FLTARR(N) + SIGN(kpar)*SQRT( 2.*kepar / mass )
  
  ;; Get perpendicular K.E.
  keperp = ke - kepar
  
  ;; Magnetic moment is mu = keperp / B
  mu = keperp / Bmin
  
  ; Create starting vector
  
  y0 = [REPLICATE(ri0, N), REPLICATE(zi0, N), $ ; All start at outboard midplane
        !DPI*2.*DINDGEN(N)/DOUBLE(N), $         ; Distribute around toroidally
        vpar]

  ; To get typical timescale, use cyclotron timescale
  
  ctime = 2.*!PI * mass / (charge * Bmax)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Evolve
  
  dt = DOUBLE(ctime)
  y = y0
  time = DOUBLE(0.)
  FOR i=0, 1000 DO BEGIN
    dydt = differential(0., y)
    y = RK4(y, dydt, 0., dt, 'differential', /double)
    time = time + dt

    PLOTS, INTERPOLATE(r1d, y[0:(N-1)]), INTERPOLATE(z1d, y[N:(2*N-1)]), color=3, PSYM=3
    ;PRINT, y
  ENDFOR
  
  IF KEYWORD_SET(output) THEN BEGIN
    ; Dump the results to a file
    
    SAVE, psi, fpol2d, r2d, z2d, psinorm, aeq, $ ; Equilibrium data
      psin, rinds, zinds, $ ; Starting flux surface
      ri0, zi0, $ ; Starting location
      ke, mu, mass, charge, $ ; Particle quantities
      y0, time, y, $  ; Starting and end positions, time elapsed
      file=output
  ENDIF
  
  STOP
END
