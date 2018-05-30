pro ncounts_lbin_all ;declares program

loadct,3 ;loads in color table
znames=['zem5','zem4','z001','z002','z003','z004','z006','z008','z010','z014','z020','z030','z040'] ;names of metallicity
mods=['0.1','0.2','0.3','0.4','0.5','0.6','0.7','0.8','0.9','1','1.1','1.2','1.3','1.5','1.5','1.7','2.1','2.3','2.5','3','4','5','6','7','8','9','10','20','30','40','50','60'] ;names of model? Doesn't appear elsewhere in code
modsN=[0.1,0.5,1,3,5,10,20,30] ;Same
!x.style=1
!y.style=1
!p.charsize=1.5
!p.charthick=3
!p.thick=3


;;;change this variable for search down to this initial mass
;lookdownto=15d0
lookdownto=0d0
;lookdownto=50d0


;Sets up time steps to run from 10^6.05 in steps of 0.1 dex
timebin=dblarr(51)
timebin(0)=10e0^(6.05)
for n=1,50 do timebin(n)=10e0^(6.05+n*0.1)-10e0^(6.05+0.1*(n-1))
print,timebin

for zn=0,12 do begin ;for each metallicity

   ;FIRST TREVOR EDIT: MAKING MODELS BE 170x100x51 = (10 star types * 17 star type bins)x100xnumber of time steps
   models = dblarr(170,100,51)

   xdum="abc"
   xdum2="abc"
   modelimf=0D
   modeltype=1
   openr,5,"/astro/store/gradscratch/tmp/tzdw/BPASS/BPASSv2.1_imf135_300/input_bpass_"+znames(zn)
   readf,5,totalimfnumber ;first line in file 5 is total stars formed
   xdum5=dblarr(2)

   while(eof(5) eq 0) do begin  ;;;will load in model location from input bpass file and then go to load in that model
      readf,5,xdum2 ;read total file name
      readf,5,xdum5 ;read IMF prob and model type (which is split in the next two lines)
      modelimf=xdum5(0)
      modeltype=xdum5(1)
;0 - merger,
;1 - normal primary,
;2 - normal secondary,
;3 - single model,
;4 - binary hmg to get black hole mass too
      mixedage=0D
      mixedimf=0D
      if(modeltype ge 1.9) then begin;i.e. 2->model is for the secondary or single, file has extra line here
         readf,5,xdum5 ;read extra line into variables for mixed imf/age (has to do with binary rejuvination?)
         mixedimf=xdum5(0)
         mixedage=xdum5(1)
         modelimf=modelimf-mixedimf ;individual model likelihood, must be positive. Enforced in next two lines
         if(modelimf lt 0d0) then modelimf=0d0
         if(mixedimf lt 0d0) then mixedimf=0d0
      endif
      initblackhole=0e0
      initperiod=0e0
      if(modeltype ge 3.9) then begin;i.e. 4
         readf,5,xdum5;initblackhole,initperiod ;!!so can put the black hole in at the end
      endif

      dummy=fltarr(73) ;Dummy is 73-entry long float array

      dum=strmid(xdum2,26,70) ;Get filename of model
      test1=FILE_TEST("/astro/store/gradscratch/tmp/tzdw/BPASS/all_mods/"+dum) ;Check if file exists.
                                ;   print,dum,test1
      if(test1 eq 1) then begin
         openr,6,"/astro/store/gradscratch/tmp/tzdw/BPASS/all_mods/"+dum ;Open model into filenumber 6
         print,dum ;Print model name
         lasttime=0e0
         lasttime2=mixedage
         count=-1
         itmax=0
         itmax2=0
         itlast=0
         dummy(10)=1e0
         itlast2=round(10e0*alog10(mixedage))-60
         initialmass=-100d0
         initialmasstest=0
         while(eof(6) eq 0 and dummy(1) le 100e9 and initialmasstest eq 0) do begin
            readf,6,dummy ;Read first timestep from model into dummy
            if(initialmass le -99d0) then begin
               initialmass=dummy(5) ;Set initial mass, which will be positive from here on out
               if(initialmass lt lookdownto) then initialmasstest=1 ;check to see if the initial mass checking is working
            endif
            count=count+1
            dt=dummy(1)-lasttime ;width of timestep
            if(dt gt 0d0) then begin;;;if actual timestep then go into loop

               it=round(10e0*alog10(dummy(1)))-60;;set up binning in time (which time bin are we in?)
               if(it le 0) then it=0
               if(it ge 50) then it=50
               itmax=it

               xh=dummy(10);;sets up some parameters
               T=dummy(3)
               gravity=alog10(6.67259d-8*1.9891d33*dummy(5)/((10d0^dummy(2))*6.9598d10)^2d0)
               COHE=(dummy(12)/3.0+dummy(14)/4.0)/dummy(11)
               logOfG= 3.676d+0 * T - 13.253d+0 ; using g_evol

               startype=0;;;;decide star type
               if(xh le 0.4 and T ge 4.45d0) then begin ;WR STAR
                  if(xh gt 1e-3) then begin
                     startype=8 ;for WNH
                  endif else begin
                     startype=9 ;!WN
                     if(COHE gt 0.03) then startype=10 ;!for WC
                  endelse
               endif else begin ;pre-WR star
                  startype=7   ;!M rsg
                  if(T ge 3.550) then startype=6 ;!K rsg
                  if(T ge 3.660) then startype=5 ;!YSG
                  if(T ge 3.9) then startype=4   ;!A bsg
                  if(T ge 4.041) then startype=3 ;!B bsg
                  if(T ge 4.48) then startype=1  ;!O bsg
                  if(gravity lt logOfG and T ge 4.519d+0) then startype=2; !Of star
               endelse



               x0=round(10e0*((dummy(4))))-20       ;log(L/Lsun)

               ;HERE IS THE BIGGEST CHANGE WHERE TREVOR IS EDITING THE CODE:
               L = dummy(4) ; get luminosity
               ;;; Basically we want every 0.1 dex of L above 3.5 and below 5.0 to add 10 to the startype so it bumps it over another set of columns
               lbin = 0
	             ;if(L ge 3.5 and L lt 5.0) then begin
	             ;   L_floor = L - 3.5 ;how much bigger than 3.5?
               ;   L_pointdex = L_floor / 0.1 ;how many 0.1's is that
               ;   lbin = CEIL(L_pointdex) ;e.g., 3.51 should be one L_bin over. 3.51-3.5 = 0.01. 0.01/0.1 = 0.1. CEIL(0.1) = 1
               ;endif else begin
               ;   if(L ge 5.0) then lbin = 16
               ;endelse
               ;uncomment the following line if you want to just reproduce the behavior from Eldridge+17
               if (L lt 4.9) then lbin = 1

	             deltacolumn = 10 * lbin ;10 types of stars, so we want to add 10 column numbers per L bin
               x1=startype-1+deltacolumn


               ;trevor changed the x1 checks to extend out to column 171. This does all the IMF adding
               if(x0 ge 0 and x0 le 100 and x1 ge 0 and x1 le 171 and initialmass gt lookdownto) then begin
                  if(it eq itlast) then begin
                     models(x1,x0,it)=models(x1,x0,it)+dt*(modelimf)
                  endif else begin
                     if((it-itlast) gt 1) then begin
                        dum1=(10d0^(6.05d0+0.1d0*itlast))
                        dum2=(10d0^(5.95d0+0.1d0*it))
                        models(x1,x0,itlast)=models(x1,x0,itlast)+(dum1-lasttime)*(modelimf)
                        models(x1,x0,it)=models(x1,x0,it)+(dummy(1)-dum2)*(modelimf)
                        for n1x=itlast+1,it-1 do begin
                           dum1=(10d0^(6.05d0+0.1d0*(n1x))-10d0^(5.95d0+0.1d0*(n1x)))
                           models(x1,x0,n1x)=models(x1,x0,n1x)+(dum1)*(modelimf)
                        endfor
                     endif else begin
                        dum1=10d0^(5.95d0+0.1d0*it)
                        models(x1,x0,itlast)=models(x1,x0,itlast)+(dum1-lasttime)*(modelimf)
                        models(x1,x0,it)=models(x1,x0,it)+(dummy(1)-dum1)*(modelimf)
                     endelse
                  endelse

                  if(modeltype ge 1.9) then begin ;i.e. 2
                     it2=round(10e0*alog10(dummy(1)+mixedage))-60
                     if(it2 le 0) then it2=0
                     if(it2 ge 50) then it2=50
                     itmax2=it2
                     if(it2 eq it2last) then begin
                        models(x1,x0,it2)=models(x1,x0,it2)+dt*(mixedimf)
                     endif else begin
                        if((it2-it2last) gt 1) then begin
                           dum1=(10d0^(6.05d0+0.1d0*it2last))
                           dum2=(10d0^(5.95d0+0.1d0*it2))
                           models(x1,x0,it2last)=models(x1,x0,it2last)+(dum1-lasttime2)*(mixedimf)
                           models(x1,x0,it2)=models(x1,x0,it2)+(dummy(1)+mixedage-dum2)*(mixedimf)
                           for n1=it2last+1,it2-1 do begin
                              dum1=(10d0^(6.05d0+0.1d0*(n1))-10d0^(5.95d0+0.1d0*(n1)))
                              models(x1,x0,n1)=models(x1,x0,n1)+(dum1)*(mixedimf)
                           endfor
                        endif else begin
                           dum1=10d0^(5.95d0+0.1d0*it2)
                           models(x1,x0,it2last)=models(x1,x0,it2last)+(dum1-lasttime2)*(mixedimf)
                           models(x1,x0,it2)=models(x1,x0,it2)+(dummy(1)+mixedage-dum1)*(mixedimf)
                        endelse
                     endelse
                  endif

               endif
               lasttime=dummy(1)
               itlast=round(10e0*alog10(dummy(1)))-60
               if(itlast le 0) then itlast=0
               if(itlast ge 50) then itlast=50
               it2last=round(10e0*alog10(dummy(1)+mixedage))-60
               if(it2last le 0) then it2last=0
               if(it2last ge 50) then it2last=50
               lasttime2=dummy(1)+mixedage
            endif
        endwhile


         close,6
      endif

   endwhile
   close,5

   openw,5,"~/ncounts_all_49"+znames(zn)+".dat"

   ;ANOTHER TREVOR EDIT: 170 'types' of star segregated by L plus time bin
   ;outputarray=dblarr(170)
   outputarray=dblarr(20)
   for n=0,50 do begin
      ;for m=0,169 do begin
      for m=0,19 do begin
         outputarray(m)=total(models(m,0:99,n))
      endfor
      printf,5,6.0+0.1*n,outputarray,FORMAT='(21E16.7)' ;FORMAT='(171E16.7)'
   endfor
   close,5




   SAVE, /VARIABLES, FILENAME = 'output-data-all.'+znames(zn)+'.idl.dat'

endfor

end
