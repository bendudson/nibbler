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
; See README file for useage instructions
;
; All quantities are in SI units except temperature in eV
;

; Cross-product of two vectors
FUNCTION cross, a, b
  RETURN, {r:(a.phi*b.z - b.phi*a.z), $
           z:(a.r*b.phi - b.r*a.phi), $
           phi:(a.z*b.r - b.z*a.r)}
END

; Dot-product of two vectors
FUNCTION dot, a, b
  RETURN, a.r*b.r + a.z*b.z + a.phi*b.phi
END

; Multiply a vector by a constant
FUNCTION mul, vec, num
  RETURN, {r:(vec.r*num), z:(vec.z*num), phi:(vec.phi*num)}
END

; Add two vectors together
FUNCTION add, a, b
  RETURN, {r:(a.r+b.r), z:(a.z+b.z), phi:(a.phi+b.phi)}
END

;Calculate equilibrium B field
FUNCTION getEqBfield, ri, zi, status=status
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d, Ar, Az, Aphi, includermp 

  status = 0
  nr = N_ELEMENTS(r1d)
  nz = N_ELEMENTS(z1d)
  IF (ri LT 0) OR (ri GT nr-1) OR (zi LT 0) OR (zi GT nz-1) THEN BEGIN
    ; Gone outside domain
    status = 1
    RETURN, 0
  ENDIF
  
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
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d, Ar, Az, Aphi, includermp
  
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

; Return the sign of a number (+1 or -1)
FUNCTION sign, x
  IF x GT 0. THEN RETURN, 1.
  RETURN, -1.
END

; Calculate 3D (r,z,phi) derivatives and second-derivatives
; given var[r,z,phi] and indices
FUNCTION get3DGradients, var, ri, zi, pi, dr, dz, dphi
  COMMON grad3d, calc, r, z, p, W, U, V
  s = SIZE(var, /dim)
  nr = s[0]
  nz = s[1]
  nphi = s[2]
  
  ri0 = ROUND(ri)
  zi0 = ROUND(zi)
  pi0 = ROUND(pi)
  
  IF N_ELEMENTS(calc) EQ 0 THEN BEGIN ; Calculate inverse matrix
     calc = 1
     
     ;; Coordinates of points in a cube (27 total)
     r = [ 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1, 0,-1, 1]
     z = [ 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1]
     p = [ 0, 0, 0, 0, 0, 0, 0, 0, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
     
     n = N_ELEMENTS(r)
     ;; Create the matrix to generate points from differentials
     A = TRANSPOSE([[fltarr(n)+1.], $ ; constant
                    [r], $            ; d/dr
                    [z], $            ; d/dz
                    [p], $            ; d/dphi
                    [0.5*r*r], $      ; d2/dr2
                    [0.5*z*z], $      ; d2/dz2
                    [0.5*p*p], $      ; d2/dp2
                    [r*z], $          ; d2/drdz
                    [r*p], $          ; d2/drdp
                    [z*p]])           ; d2/dzdp
  
     ;; Calculate SVD
     SVDC, A, W, U, V
  ENDIF
  
  n = N_ELEMENTS(r)
  ; Get the values at each point
  data = DBLARR(n)
  FOR i=0, n-1 DO data[i] = var[((ri0+r[i]) > 0) < nr, $
                                ((zi0+z[i]) > 0) < nz, $
                                (((pi0+p[i]) MOD nphi) + nphi) MOD nphi]
  

  res = SVSOL(U,W,V, data)
  
  RETURN, {f:res[0], ddr:res[1]/dr, ddz:res[2]/dz, ddp:res[3]/dphi, $
           d2dr2:res[4]/(dr*dr), d2dz2:res[5]/(dz*dz), d2dp2:res[6]/(dphi*dphi), $
           d2drdz:res[7]/(dr*dz), d2drdp:res[8]/(dr*dphi), d2dzdp:res[9]/(dz*dphi)}
END

; Add cartesian vectors
FUNCTION addCart, a, b
  ac = toCart(a)
  bc = toCart(b)

  RETURN, {x:(ac.x+bc.x), y:(ac.y+bc.y), z:(ac.z+bc.z)}
END

;; Calculate time derivatives
FUNCTION differential, t, y
  COMMON particle_com, N, mu, mass, charge
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d, Ar, Az, Aphi, includermp
  
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
  
  evolve = REPLICATE(DOUBLE(1.), N)
  FOR i=0, N-1 DO BEGIN
     eqB = getEqBfield(ri[i], zi[i], status=status)
     IF status THEN BEGIN
       evolve[i] = 0. ; Turn off evolution for this particle
     ENDIF ELSE BEGIN
     
       Br[i] = eqB.Br
       Bz[i] = eqB.Bz
       Bphi[i] = eqB.Bphi
       
       dBrdr[i] = eqB.dBrdr
       dBrdz[i] = eqB.dBrdz
       
       dBzdr[i] = eqB.dBzdr
       dBzdz[i] = eqB.dBzdz
       
       dBphidr[i] = eqB.dBphidr
       dBphidz[i] = eqB.dBphidz
     ENDELSE
  ENDFOR
  
  IF includermp THEN BEGIN
    ; Add the perturbed B from coils
    s = SIZE(Ar, /dim)
    IF N_ELEMENTS(s) NE 3 THEN STOP
    nphi = s[2] ; Number of points in phi
    dphi = !DPI*2. / DOUBLE(nphi)
    
    FOR i=0, N-1 DO BEGIN
      ; Calculate gradients at this location (Using SVD)
      IF evolve[i] THEN BEGIN
        Arg = get3DGradients(Ar,   ri[i], zi[i], phi[i]/dphi, drdi[i], dzdi[i], dphi)
        Azg = get3DGradients(Az,   ri[i], zi[i], phi[i]/dphi, drdi[i], dzdi[i], dphi)
        Apg = get3DGradients(Aphi, ri[i], zi[i], phi[i]/dphi, drdi[i], dzdi[i], dphi)
       
        ; Take curl of A
        Br[i] = Br[i] + Azg.ddp/r[i] - Apg.ddz
        Bz[i] = Bz[i] + ( Apg.f + r[i]*Apg.ddr - Arg.ddp ) / r[i]
        Bphi[i] = Bphi[i] + Arg.ddz - Azg.ddr
        
        ; Gradients of B
        
        dBrdr[i] = dBrdr[i] + Azg.d2drdp/r[i] - Azg.ddp/(r[i]*r[i]) - Apg.d2drdz
        dBrdz[i] = dBrdz[i] + Azg.d2dzdp/r[i] - Apg.d2dz2
        dBrdphi[i] = (Azg.d2dp2/r[i] - Apg.d2dzdp)/r[i]
        
        dBzdr[i] = dBzdr[i] - ( Apg.f + r[i]*Apg.ddr - Arg.ddp ) / (r[i]*r[i]) + $
          ( Apg.ddr + Apg.ddr + r[i]*Apg.d2dr2 - Arg.d2drdp ) / r[i]
        dBzdz[i] = dBzdz[i] + ( Apg.ddz + r[i]*Apg.d2drdz - Arg.d2dzdp ) / r[i]
        dBzdphi[i] = (Apg.ddp + r[i]*Apg.d2drdp - Arg.d2dp2)/(r[i]*r[i])
        
        dBphidr[i] = dBphidr[i] + Arg.d2drdz - Azg.d2dr2
        dBphidz[i] = dBphidz[i] + Arg.d2dz2 - Azg.d2drdz
        dBphidphi[i] = (Arg.d2dzdp - Azg.d2drdp)/r[i]
      ENDIF
    ENDFOR
  ENDIF
    
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

  return, [(v.r / drdi)*evolve, (v.z / dzdi)*evolve, $
           (v.phi / r)*evolve, dvpardt*evolve]
END

PRO nibbler, Nparticles=Nparticles, shot=shot, electron=electron, $
             temp=temp, psin=psin, $
             kpar=kpar, output=output, equil=equil, odd=odd, $
             current=current, rmp=rmp, $
             runfor=runfor, dtmul=dtmul
  COMMON particle_com, N, mu, mass, charge
  COMMON equil_com, dctpsi, dctfpol, r1d, z1d, Ar, Az, Aphi, includermp

  IF NOT KEYWORD_SET(Nparticles) THEN Nparticles = 1 ; Number of particles
  IF NOT KEYWORD_SET(temp) THEN temp = 200 ; Temperature in eV
  IF NOT KEYWORD_SET(psin) THEN psin = 0.9 ; Starting psi location
  IF NOT KEYWORD_SET(kpar) THEN kpar = 0.3 ; Vpar^2 / V^2
  IF NOT KEYWORD_SET(current) THEN current = 4e3 ; Coil current in Amps
  IF NOT KEYWORD_SET(runfor) THEN runfor=1e-3 ; Run time in seconds
  IF NOT KEYWORD_SET(dtmul) THEN dtmul = 1.  ; Multiply the timestep by this
  
  includermp = 0
  IF KEYWORD_SET(rmp) THEN includermp = 1
  
  N = Nparticles
  
  IF KEYWORD_SET(equil) THEN BEGIN
     ; Load a previously calculated equilibrium
     
     RESTORE, 'equil.idl'
     
     s = SIZE(psi, /dim)
     nr = s[0]
     nz = s[1]
  ENDIF ELSE BEGIN
     IF KEYWORD_SET(shot) THEN BEGIN
        ;; Fetch MAST equilibrium from IDAM
        
        PRINT, "Sorry, code not written yet."
        RETURN
     ENDIF ELSE BEGIN
        ;; Read in neqdsk file to get Psi
        file = "g014220.00200"
        
        PRINT, "Reading G-EQDSK file"
        grid = read_neqdsk(file)
        
        nr = grid.nx            ; Number of points in R
        nz = grid.ny            ; Number of points in Z
        psi = grid.psi          ; Poloidal flux on [nr, nz] mesh
        
        fpol1d = grid.fpol      ; f on uniform psi grid
        
        r2d = grid.r
        z2d = grid.z
        
        r1d = REFORM(r2d[*,0])
        z1d = REFORM(z2d[0,*])
     ENDELSE
     
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; By this point, need
     ;; nr, nz
     ;; psi[nr,nz], fpol1d[psi], r1d[nr], z1d[nz]

     ;; Analyse the equilibrium
     aeq = analyse_equil( psi, r1d, z1d )
     
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
     
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Generate RMP coil field
     
     PRINT, "Generating vector potential from coils..."
     
     ;;;; Define coil sets for MAST
     lower = {r1:1.311, z1:-0.791, r2:1.426, z2:-0.591, n:6, dphi:2.*!PI/12.}
     upper = {r1:1.311, z1:0.791, r2:1.426, z2:0.591, n:6, dphi:2.*!PI/12.}
     
     nphi = 64
     dz = 2.*!PI / FLOAT(nphi)
     
     rpos   = FLTARR(nr, nz, nphi)
     zpos   = rpos
     phipos = rpos
     FOR i=0, nphi-1 DO BEGIN
        rpos[*,*,i]   = r2d
        zpos[*,*,i]   = z2d
        phipos[*,*,i] = FLOAT(i)*dz
     ENDFOR
     
     ; Calculate for current of 1A (multiplied up after)
     I = 1.
     pos = {r:rpos, z:zpos, phi:phipos}
     A = AfromCoilSet(lower, I, pos, /fast)
     IF KEYWORD_SET(odd) THEN BEGIN
       A = addCart(A, AfromCoilSet(upper, I, pos, /fast))
     ENDIF ELSE BEGIN
       A = addCart(A, AfromCoilSet(upper, -I, pos, /fast))
     ENDELSE
     
     ; Convert to polar coordinates
     Ar   = A.x * COS(phipos) + A.y * SIN(phipos)
     Aphi = A.y * COS(phipos) - A.x * SIN(phipos)
     Az   = A.z
     
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Save file
     
     PRINT, "Saving equilibrium data to 'equil.idl'"
     
     SAVE, psi, psinorm, fpol1d, fpol2d, $  ; psi and fpol
       r1d, z1d, r2d, z2d, $   ; Radius, height in 1 and 2D
       dctpsi, dctfpol, $      ; DCT of psi and fpol
       aeq, $                  ; Analysis of equilibrium
       Ar, Az, Aphi, $         ; RMP field
       file="equil.idl"
  ENDELSE
  
  ; Multiply the RMP coil field
  Ar   = current * Ar
  Az   = current * Az
  Aphi = current * Aphi
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plotting
  
  ;; Plot contour
  nlev = 100
  minf = MIN(psi)
  maxf = MAX(psi)
  levels = findgen(nlev)*(maxf-minf)/FLOAT(nlev-1) + minf
  safe_colors, /first
  CONTOUR, psi, r2d, z2d, levels=levels, color=1, /iso, xstyl=1, ysty=1
  
  ;; Overplot the separatrices, O-points
  oplot_critical, psi, r1d, z1d, aeq
  
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
  
  dt = DOUBLE(ctime) * dtmul
  y = y0
  time = DOUBLE(0.)
  data = [[y0]]
  tarr = [0.0]
  nsteps = ROUND(runfor / dt)
  PRINT, "TIMESTEP = ", dt, " seconds"
  PRINT, "NUMBER OF STEPS: ", nsteps
  FOR i=0, nsteps DO BEGIN
    dydt = differential(0., y)
    y = RK4(y, dydt, 0., dt, 'differential', /double)
    time = time + dt

    PLOTS, INTERPOLATE(r1d, y[0:(N-1)]), INTERPOLATE(z1d, y[N:(2*N-1)]), color=3, PSYM=3
    ;PRINT, y
    IF i MOD 10 EQ 0 THEN BEGIN
      data = [[data], [y]]
      tarr = [tarr, time]
      WRITEU, -1, 13, "Progress: "+STR(100.*FLOAT(i)/nsteps)+"%"
    ENDIF
  ENDFOR
  
  IF KEYWORD_SET(output) THEN BEGIN
    ; Dump the results to a file
    
    SAVE, psi, fpol2d, r2d, z2d, psinorm, aeq, $ ; Equilibrium data
      psin, rinds, zinds, $ ; Starting flux surface
      ri0, zi0, $ ; Starting location
      ke, mu, mass, charge, $ ; Particle quantities
      tarr, data, $  ; time and 
      file=output
  ENDIF
  
  STOP
END
