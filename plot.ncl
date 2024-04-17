;************************************************
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

  load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"

;************************************************
;
;************************************************
begin

  hemi="NH"
  ;hemi="SH"

  ;minlat = 66.0
  minlat = 76.0
  ;minlat = 89.8
  ;minlat=88.0
  vmag = 0.75
  vthn = 0.008

;************************************************
; read in netCDF file
;************************************************

  ;fname="ice.mx025.bilin2d.nc"
  ;gr=addfile("/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/025/tripole.mx025.nc","r")
  ;varname="bilin2d"

  fname="ice.mx100.rgbilin2d.nc"
  gr=addfile("/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/100/tripole.mx100.nc","r")
  varname="rgbilin2d"

  wkstitle="native"
  ifile=addfile(fname,"r")
  plottitle=fname

  print(fname)

  ut = ifile->$varname$(0,:,:)
  vt = ifile->$varname$(1,:,:)
  print("ut min = "+min(ut)+" ut max = "+max(ut))
  print("vt min = "+min(vt)+" vt max = "+max(vt))

  ut@lon2d=gr->lonCt
  ut@lat2d=gr->latCt
  vt@lon2d=gr->lonCt
  vt@lat2d=gr->latCt

  ; ; use cice's internal angle variable
  ; urot = ut
  ; vrot = vt
  ; urot = cos(rotang)*ut - sin(rotang)*vt
  ; vrot = sin(rotang)*ut + cos(rotang)*vt

  ; urot = where(urot@lat2d .ge. minlat, urot, 0.0)
  ; vrot = where(vrot@lat2d .ge. minlat, vrot, 0.0)
  ; print("urot min = "+min(urot)+" urot max = "+max(urot))
  ; print("vrot min = "+min(vrot)+" vrot max = "+max(vrot))

  ;printVarSummary(urot)
  ;printVarSummary(vrot)

  ut = where(ut@lat2d .ge. minlat, ut, 0.0)
  vt = where(vt@lat2d .ge. minlat, vt, 0.0)

;************************************************
; add a poly line showing the polar seem along
; great circle path from (lon1,lat1) to (lon2,lat2)
;************************************************

   lon1=60.00
   lon2=-120.00
   lat1=60
   lat2=60

;************************************************
; define the map
;************************************************
  wkstype="X11"
  ;wkstype="png"
  wkstype@wkWidth=800
  wkstype@wkHeight=800
  wks= gsn_open_wks (wkstype,wkstitle)              ; send graphics to X11 window

  plot = new(2,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@vcRefMagnitudeF       = vmag              ; vec magnitude
  res@vcRefLengthF          = 0.10              ; size of ref vector
  ;res@vcMonoLineArrowColor  = True             ; multiple colors desired
  res@vcMinDistanceF        = vthn              ; thins arrows near pole
  ;res@vcLineArrowThicknessF = 1.5
  ;res@vcGlyphStyle          = "CurlyVector"
  ;res@vcPositionMode = "ArrowTail"

  res@gsnPolar   =  hemi                         ; specify the hemisphere
  if(hemi .eq. "NH")then
   res@mpMinLatF    =  minlat                          ; minimum lat to plot
   res@mpCenterLonF = -45
  end if
  if(hemi .eq. "SH")then
   res@mpMaxLatF    = -50                            ; minimum lat to plot
   res@mpCenterLonF = 0.0
  end if

;************************************************
; Resources for the polyline
;************************************************

  res@mpGreatCircleLinesOn = True
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF = 2.0            ; line thickness
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =  1            ; dashes
  ;pres@tfPolyDrawOrder = "postDraw"

; ***************************************************
; panel resources
; **************************************************

  resP                       = True                 ; modify the panel plot
  resP@gsnMaximize           = True                 ; maximize in frame
  resP@gsnFrame              = False                ; don't advance frame yet
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelLabelBar      = True                 ; turn on a common labelbar for the entire panel plot

  ;resP@gsnPanelMainFontHeightF = 0.015           ; new resource added in NCL V6.4.0
  resP@gsnPanelMainString = plottitle

;************************************************
;
;************************************************

  plot(0) = gsn_csm_vector_map_polar(wks,ut,vt,res)    ; create the plot
  ;plot(0) = gsn_csm_contour_map_polar(wks,urot,res)    ; create the plot

;************************************************
; make the panel
;************************************************

  ;gsn_panel(wks,plot,(/1,2/),resP)
  gsn_panel(wks,plot,(/1,1/),resP)
  frame(wks)

exit
end