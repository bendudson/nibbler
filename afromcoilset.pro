; Calculate the Apar due to a given set of coils


; Convert polar to cartesian coordinates
FUNCTION polarToCart, vp
  RETURN, {x:vp.r*cos(vp.phi), y:vp.r*sin(vp.phi), z:vp.z}
END

; Convert cartesian to polar coordinates
FUNCTION cartToPolar, cp
  RETURN, {r:SQRT(cp.x^2 + cp.y^2), phi:atan(cp.y, cp.x), z:cp.z}
END

; Convert to cartesian
FUNCTION toCart, p
  CATCH, err
  IF err EQ 0 THEN BEGIN
    i = p.x
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; No 'x' component - convert 
    RETURN, polarToCart(p)
  ENDELSE
  ; Already in cartesian
  RETURN, p
END

; Convert to polar
FUNCTION toPolar, p
  CATCH, err
  IF err EQ 0 THEN BEGIN
    i = p.phi
  ENDIF ELSE BEGIN
    CATCH, /cancel
    ; No 'phi' component - convert 
    RETURN, cartToPolar(p)
  ENDELSE
  ; Already in polar
  RETURN, p
END

; Distance between two points
FUNCTION distance, p1, p2
  c1 = toCart(p1)
  c2 = toCart(p2)
  
  RETURN, SQRT((c1.x - c2.x)^2 + $
               (c1.y - c2.y)^2 + $
               (c1.z - c2.z)^2)
END

FUNCTION linediff, I, y
  COMMON linecom, c1, c2, c
  
  ; i between 0 and 1
  
  ; Current position
  ipos = {x:( i*c2.x + (1.-i)*c1.x ), $
          y:( i*c2.y + (1.-i)*c1.y ), $
          z:( i*c2.z + (1.-i)*c1.z )}
  
  d = distance(ipos, c)
  
  RETURN, [1./d]
END

; Wire from p1 to p2, carrying current 
; Get A at pos
FUNCTION AfromLine, p1, p2, current, pos, fast=fast 
  COMMON linecom, c1, c2, c
  
  ; Convert all coordinates to cartesian
  c1 = toCart(p1)
  c2 = toCart(p2)
  cin  = toCart(pos)
  
  len = distance(c1, c2) ; length of the wire
  
  Ivec = {x:current*(c2.x - c1.x)/len, $
          y:current*(c2.y - c1.y)/len, $
          z:current*(c2.z - c1.z)/len}
  
  IF KEYWORD_SET(fast) THEN BEGIN
    ; Use a fixed number of Simpson rule steps
    
    n = 2 ; Must be even
    h = len / FLOAT(n)
    FOR i=0, n DO BEGIN
      f = FLOAT(i) / FLOAT(n)
      ; Position along wire
      ipos = {x:( f*c2.x + (1.-f)*c1.x ), $
              y:( f*c2.y + (1.-f)*c1.y ), $
              z:( f*c2.z + (1.-f)*c1.z )}
      
      ; Distance
      d = distance(ipos, cin)
      
      IF i EQ 0 THEN BEGIN
        integral = 1. / d
      ENDIF ELSE IF i EQ n THEN BEGIN
        integral = integral + 1. / d
      ENDIF ELSE IF i MOD 2 EQ 1 THEN BEGIN
        integral = integral + 4. / d
      ENDIF ELSE BEGIN
        integral = integral + 2. / d
      ENDELSE
    ENDFOR
    integral = integral * (h / 3.) * 1.e-7
    result = {x:Ivec.x*integral, $
              y:Ivec.y*integral, $
              z:Ivec.z*integral}
  ENDIF ELSE BEGIN
    result = cin
    i = 0L
    REPEAT BEGIN
      c = {x:cin.x[i], y:cin.y[i], z:cin.z[i]}
      
      ; Integrate along the line
      a0 = [0.]
      a = LSODE(a0, 0., 1., 'linediff', lstat)
      a = a * 1.e-7 ; mu_0 / 4pi
      
      result.x[i] = Ivec.x*a[0]
      result.y[i] = Ivec.y*a[0]
      result.z[i] = Ivec.z*a[0]
      
      i = i + 1L
    ENDREP UNTIL i EQ N_ELEMENTS(cin.x)
  ENDELSE
  
  RETURN, result
END

FUNCTION AfromArc, p1, phi, current, pos, fast=fast 
  ; For now turn into a series of lines
  
  ps = p1
  pe = p1
  
  n = 10
  a = toCart(pos)
  
  dphi = phi / FLOAT(n)
  FOR i=0, n-1 DO BEGIN
    pe.phi = ps.phi + dphi
    
    a1 = AfromLine(ps, pe, current, pos, fast=fast)
    
    a.x = a.x + a1.x
    a.y = a.y + a1.y
    a.z = a.z + a1.z
    
    ps = pe
  ENDFOR
  
  RETURN, a
END

; Add cartesian vectors
FUNCTION addCart, a, b
  ac = toCart(a)
  bc = toCart(b)

  RETURN, {x:(ac.x+bc.x), y:(ac.y+bc.y), z:(ac.z+bc.z)}
END


; Add a coil, giving two corners
FUNCTION AfromCoil, p1, p2, current, pos, fast=fast
  
  ; Corners
  c0 = toPolar(p1)
  c2 = toPolar(p2)
  dphi = c2.phi - c0.phi
  
  c1 = {r:c0.r, z:c0.z, phi:c2.phi}
  c3 = {r:c2.r, z:c2.z, phi:c0.phi}
  
  A = AfromLine(c0, c1, current, pos, fast=fast)
  A = addCart(A, AfromLine(c1, c2, current, pos, fast=fast))
  A = addCart(A, AfromLine(c2, c3, current, pos, fast=fast))
  A = addCart(A, AfromLine(c3, c0, current, pos, fast=fast))
  
  RETURN, A
END

FUNCTION AfromCoilSet, set, current, pos, fast=fast
  
  shift = 2.*!PI / FLOAT(set.n)
  
  ; Go through the set of n coils
  FOR i=0, set.n-1 DO BEGIN
    c0 = {r:set.r1, z:set.z1, phi:(shift*i)}
    c1 = {r:set.r2, z:set.z2, phi:(c0.phi + set.dphi)}
    
    dA = AfromCoil(c0, c1, current*(-1.)^i, pos, fast=fast)
    IF i EQ 0 THEN A = dA ELSE A = addCart(A, dA)
  ENDFOR
  
  RETURN, A
END
