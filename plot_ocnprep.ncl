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

;************************************************
; read in netCDF file
;************************************************

 ; these are EW, on Ct
  ;ew = 1
  ;fnames = (/"ice.mx025.bilin2d.nc", "ice.mx100.rgbilin2d.nc"/)
 ;ftypes = (/"native Ct", "regrid Ct"/)
  ;icelonname="lonCt"
  ;icelatname="latCt"

  ; these are IJ, on Bu
  ew = 0
  fnames = (/"ice.nc", "ice.mx100.nc"/)
  plottitle = (/"source","dest"/)
  icelonname="lonBu"
  icelatname="latBu"

  fsrc = addfile(fnames(0),"r")
  fdst = addfile(fnames(1),"r")
  print(fnames(0))
  print(fnames(1))
  gsrc = addfile("/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/025/tripole.mx025.nc","r")
  gdst = addfile("/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data/100/tripole.mx100.nc","r")

  ;usrc = fsrc->bilin2d(0,:,:)
  ;vsrc = fsrc->bilin2d(1,:,:)
  ;rotsrc = gsrc->anglet(:,:)

  usrc = fsrc->uvel(:,:)
  vsrc = fsrc->vvel(:,:)
  rotsrc = gsrc->angle(:,:)

  usrc@lon2d=gsrc->$icelonname$
  usrc@lat2d=gsrc->$icelatname$
  vsrc@lon2d=gsrc->$icelonname$
  vsrc@lat2d=gsrc->$icelatname$
  vthnsrc = 0.015

  ;udst = fdst->rgbilin2d(0,:,:)
  ;vdst = fdst->rgbilin2d(1,:,:)
 ;rotdst = gdst->anglet(:,:)

  udst = fdst->uvel(:,:)
  vdst = fdst->vvel(:,:)
  rotdst = gdst->angle(:,:)

  udst@lon2d=gdst->$icelonname$
  udst@lat2d=gdst->$icelatname$
  vdst@lon2d=gdst->$icelonname$
  vdst@lat2d=gdst->$icelatname$
  vthndst = 0.015

 ; if vectors are ij oriented
  if (ew .eq. 0) then
    urot = udst ; copy metadata
    vrot = vdst
    rotang = rotdst
    urot =  cos(rotang)*udst - sin(rotang)*vdst
    vrot =  sin(rotang)*udst + cos(rotang)*vdst
    udst = where(urot@lat2d .ge. minlat, urot, 0.0)
    vdst = where(vrot@lat2d .ge. minlat, vrot, 0.0)

    urot := usrc
    vrot := vsrc
    rotang := rotsrc
    urot =  cos(rotang)*usrc - sin(rotang)*vsrc
    vrot =  sin(rotang)*usrc + cos(rotang)*vsrc
    usrc = where(urot@lat2d .ge. minlat, urot, 0.0)
    vsrc = where(vrot@lat2d .ge. minlat, vrot, 0.0)
  end if


  ;urot = cos(rotang)*ut - sin(rotang)*vt
  ;vrot = sin(rotang)*ut + cos(rotang)*vt
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
  wks= gsn_open_wks (wkstype,"test")              ; send graphics to X11 window

  plot = new(2,graphic)

  res                      = True        ; plot mods desired
  res@gsnDraw              = False       ; do not draw until the end
  res@gsnFrame             = False       ; do not automatically advance 'frame'

  res@vcRefMagnitudeF       = vmag              ; vec magnitude
  res@vcRefLengthF          = 0.10              ; size of ref vector
  ;res@vcMonoLineArrowColor  = True             ; multiple colors desired
  res@vcMinDistanceF        = vthnsrc              ; thins arrows near pole
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
  pres@tfPolyDrawOrder = "postDraw"

; ***************************************************
; panel resources
; **************************************************

  resP                       = True                 ; modify the panel plot
  resP@gsnMaximize           = True                 ; maximize in frame
  resP@gsnFrame              = False                ; don't advance frame yet
  resP@gsnPanelYWhiteSpacePercent = 5
  resP@gsnPanelYWhiteSpacePercent = 5

  resP@gsnPanelMainFontHeightF = 0.015           ; new resource added in NCL V6.4.0
  resP@gsnPanelMainString = "ocnprepf90"
  resP@gsnPanelFigureStrings = (/plottitle(0),plottitle(1)/)

;************************************************
;
;************************************************
  res@vcMinDistanceF        = vthnsrc              ; thins arrows near pole
  plot(0) = gsn_csm_vector_map_polar(wks,usrc,vsrc,res)    ; create the plot
  res@vcMinDistanceF        = vthndst              ; thins arrows near pole
  plot(1) = gsn_csm_vector_map_polar(wks,udst,vdst,res)

;************************************************
; make the panel
;************************************************

  gsn_panel(wks,plot,(/2,1/),resP)
  frame(wks)

exit
end