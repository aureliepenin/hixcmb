pro load_my_colors, PLOT=PLOT

; 0 : black
; 1 : white
; 2 : red
; 3 : green
; 4 : blue
; 5 : purple
; 6 : yellow
; 7 : grey
; 8 : grey 2
; 9 : orange
; 10 : cyan
; 11 : marron
; 12 : olive
; 13 : bleu
; 14 : vert
; 15 : turquoise
; 16 : violet
; 17 : gold
; 18 : vert


TVLCT, [0,255,255,0,0,255,255,200,130,255,0,210,128,0,0,0,255,255,0], $
       [0,255,0,255,0,0,255,200,130,114,255,105,128,128,128,200,20,215,200], $
       [0,255,0,0,255,255,0,200,130,0,255,30,0,128,0,200,147,0,0]




IF keyword_set(PLOT) THEN BEGIN
    Ncol = 19

    SET_PLOT, 'ps'
    DEVICE, FILENAME='~/library/useful/colors_test.eps', /ENCAPSULATED, /COLOR, xs=Ncol, ys=5

    @symbols_plots_casey_eric
    USERSYM, cir, /FILL

    PLOT, [0], [0], xs=5, ys=5, xr=[-1,Ncol], yr=[-1,2], COLOR=0 

    FOR k=0, Ncol-1 DO BEGIN
        USERSYM, cir, /FILL
        OPLOT, [k], [0], PSYM=8, SYMSIZE=2.5, COLOR=k
        XYOUTS, k, 1, STRING(k,F='(I2)'), ALIGN=0.5, CHARSIZE=1.2, CHARTHICK=7, /DATA
    ENDFOR

    DEVICE, /CLOSE
    SET_PLOT, 'x'

ENDIF



end
