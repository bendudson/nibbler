; Categorises points in psi[nx,ny]
;
; 1 = inside core
; 0 = outside core

FUNCTION categorise, psi, psinorm=psinorm
  s = SIZE(psi, /dim)
  IF N_ELEMENTS(s) NE 2 THEN BEGIN
     PRINT, "ERROR: psi passed to categorise needs to be 2D"
     RETURN, 0
  ENDIF
  nx = s[0]
  ny = s[1]
  
  ;; Analyse the equilibrium
  
  r = (FINDGEN(nx)+1)/FLOAT(nx) ; create dummy r and z coordinates
  z = (FINDGEN(ny) / FLOAT(ny)) - 0.5
  a = analyse_equil( psi, r, z )
  
  ;; Get the primary O-point
  opt_ri = a.opt_ri[a.primary_opt]
  opt_zi = a.opt_zi[a.primary_opt]
  opt_psi = a.opt_f[a.primary_opt]
  ;; Get innermost x-point
  xpt_psi = a.xpt_f[a.inner_sep]
  
  ;; Get a contour 
  level = 0.98*xpt_psi + 0.02*opt_psi
  contour_lines, psi, findgen(nx), findgen(ny), levels=[level], $
                 path_info=info, path_xy=xy

  IF N_ELEMENTS(info) GT 1 THEN BEGIN
     ;; Find the surface closest to the o-point
      
     ind = closest_line(info, xy, opt_ri[primary_opt], opt_zi[primary_opt])
     info = info[ind]
  ENDIF ELSE info = info[0]
  STOP
  ; Now categorise all points 
  
  core = FLTARR(nx, ny)
  FOR i=0, nx-1 DO BEGIN
     FOR j=0, ny-1 DO BEGIN
        ;; Get crossing points for line to O-point and separatrix
        data = line_crossings([i,opt_ri], [j,opt_zi], 0, ncross=ncross)
        IF ncross MOD 2 EQ 0 THEN core[i,j] = 1 ; inside
     ENDFOR
  ENDFOR
  
  ; Calculate normalised psi, 1 outside core
  psinorm = core * (psi - opt_psi) / (xpt_psi - opt_psi) + (1. - core)

  RETURN, core
END


