; Categorises points in psi[nx,ny]
;
; 1 = inside core
; 0 = outside core

FUNCTION categorise, psi, psinorm=psinorm, aeq=aeq
  s = SIZE(psi, /dim)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
     PRINT, "ERROR: psi passed to categorise needs to be 2D"
     RETURN, 0
  ENDIF
  nx = s[0]
  ny = s[1]
  
  IF NOT KEYWORD_SET(aeq) THEN BEGIN
     ;; Analyse the equilibrium
     r = (FINDGEN(nx)+1)/FLOAT(nx) ; create dummy r and z coordinates
     z = (FINDGEN(ny) / FLOAT(ny)) - 0.5
     aeq = analyse_equil( psi, r, z )
  ENDIF
  
  ;; Get the primary O-point
  opt_ri = aeq.opt_ri[aeq.primary_opt]
  opt_zi = aeq.opt_zi[aeq.primary_opt]
  opt_psi = aeq.opt_f[aeq.primary_opt]
  ;; Get innermost x-point
  xpt_psi = aeq.xpt_f[aeq.inner_sep]
  
  ;; Get a contour 
  level = 0.98*xpt_psi + 0.02*opt_psi
  contour_lines, psi, findgen(nx), findgen(ny), levels=[level], $
                 path_info=info, path_xy=xy

  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ;; Find the surface closest to the o-point
      
     ind = closest_line(info, xy, opt_ri, opt_zi)
     info = info[ind]
  ENDIF ELSE info = info[0]
  rinds = xy[0,info.offset:(info.offset+info.n-1)]
  zinds = xy[1,info.offset:(info.offset+info.n-1)]
  
  ; Now categorise all points 
  
  core = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
     FOR j=0, ny-1 DO BEGIN
        ;; Get crossing points for line to O-point and separatrix
        data = line_crossings([i,opt_ri], [j,opt_zi], 0, $
                              rinds, zinds, 1, ncross=ncross)
        IF ncross MOD 2 EQ 0 THEN core[i,j] = 1 ; inside
     ENDFOR
  ENDFOR
  
  ; Calculate normalised psi, 1 outside core
  psinorm = core * (psi - opt_psi) / (xpt_psi - opt_psi) + (1. - core)

  RETURN, core
END


