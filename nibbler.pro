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
FUNCTION cross(a, b)
  RETURN {r:(a.phi*b.z - b.phi*a.z), $
          z:(a.r*b.phi - b.r*a.phi), $
          phi:(a.z*b.r - b.z*a.r)}
END

FUNCTION dot(a, b)
  RETURN, a.r*b.r + a.z*b.z + a.phi*b.phi
END

FUNCTION mul(vec, num)
  RETURN {r:(vec.r*num), z:(vec.z*num), phi:(vec.phi*num)}
END

FUNCTION add(a, b)
  RETURN {r:(a.r+b.r), z:(a.z+b.z), phi:(a.phi+b.phi)}
END

;; Calculate time derivatives
PRO deriv, t, y
  COMMON particle_com, N, ke, mu, mass, charge
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
  drdi = INTERPOLATE(DERIV(REFORM(r1d), ri)
  dzdi = INTERPOLATE(DERIV(REFORM(z1d), ri)
  
  FOR i=0, N-1 DO BEGIN
     g = local_gradient(dctfpol, ri[i], zi[i])
     Bphi[i] = g.f / r[i]
     dBphidr[i] = (g.dfdr/drdi[i])/r[i] - Bphi[i]/r[i]
     dBphidz[i] = (g.dfdz/dzdi[i])/r[i]
     
     ;; For poloidal field derivatives, need second
     ;; derivatives of psi
     
     g = EvalCosP(dctpsi, x0=ri[i], y0=zi[i])
     ;; g contains [ psi, dpsidr, dpsidz, d2psidr2, d2pdidz2,
     ;; d2psidrdz ]
     
     ;; NOTE: CHECK THESE SIGNS - A MINUS SOMEWHERE
     Br[i] = g[2] / dzdi
     Bz[i] = g[1] / drdi
     
     dBrdr[i] = g[5] / (drdi * dzdi)
     dBrdz[i] = g[4] / (dzdi^2)
     
     dBzdr[i] = g[3] / (drdi^2)
     dBzdz[i] = g[5] / (drdi * dzdi)
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
  lbrdz   = (dBrdz - dBdz*Br/B)/B
  lbrdphi = (dBrdphi - dBdphi*Br/B)/B
  
  lbzdr   = (dBzdr - dBdr*Bz/B)/B
  ;lbzdz   = (dBzdz - dBdz*Bz/B)/B
  lbzdphi = (dBzdphi - dBdphi*Bz/B)/B
  
  lbphidr   = (dBphidr - dBdr*Bphi/B)/B
  lbphidz   = (dBphidz - dBdz*Bphi/B)/B
  ;lbphidphi = (dBphidphi - dBdphi*Bphi/B)/B
  
  bvec = {r:(Br/B), z:(Bz/B), phi:(Bphi/B)}
  curlb = {r:(lbzdphi/r - lbphidz), phi:(lbrdz - lbzdr), $
           z:((bvec.phi+r*lbphidr - lbrdphi)/r)}
  
  kappa = cross_prod( curlb, bvec )  ; (b dot grad) b

  ;; Calculate perpendicular drift velocity
  invwcj = mass / (charge * B)  ; 1 / wcj
  
  vd = mul( cross_prod( bvec, add( mul(kappa, vpar^2), mul(gradB, mu) )), invwcj )
  
  ;; Add parallel velocity
  v = add( vd, mul(bvec, vpar) )
  
  ; Calculate parallel acceleration (mirror force)
  dvpardt = -mu * dot( bvec, gradB ) / mass
  
  return [v.r / drdi, v.z / dzdi, v.phi, dvpardt]
END

PRO nibbler, shot=shot
  COMMON particle_com, N, ke, mu
  COMMON equil_com, dctpsi, dctfpol, r, z

  N = 1
  
  IF KEYWORD_SET(shot) THEN BEGIN
     ;; Fetch MAST equilibrium from IDAM
     
     PRINT, "Sorry, code not written yet."
     RETURN
  ENDIF ELSE BEGIN
     ;; Read in neqdsk file to get Psi
     file = ""
     
     PRINT, "Reading G-EQDSK file"
     grid = read_neqdsk(file)
     
     nr = grid.nx               ; Number of points in R
     nz = grid.ny               ; Number of points in Z
     psi = grid.psi             ; Poloidal flux on [nr, nz] mesh
     
     ;; Categorise points inside and outside core
     core = categorise(psi, psinorm=psinorm)     ; 1 where in core, 0 otherwise
     
     ;; Now interpolate f onto grid points
     npgrid = (FINDGEN(N_ELEMENTS(g.fpol))/(N_ELEMENTS(g.fpol)-1))
     fpol = INTERPOL(grid.fpol, npgrid, psinorm, /spline)
     
  ENDELSE
  
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; By this point, need
  ;; nr, nz
  ;; psi[nr,nz], fpol[nr, nz]
  
  PRINT, "Calculating DCT of Psi..."
  DCT2Dslow, psi, dctpsi
  
  PRINT, "Calculating DCT of fpol..."
  DCT2Dslow, psi, dctfpol
  
  ;; Generate initial spread of kinetic energies and magnetic moment
  
  ke = FLTARR(N)
  mu = FLTARR(N)
  
  
END
