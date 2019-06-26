pro plot1cmdmake

znames=['z002','z014','z0004']


timebin=dblarr(51)
timebin(0)=10e0^(6.05)
for n=1,50 do timebin(n)=10e0^(6.05+n*0.1)-10e0^(6.05+0.1*(n-1))
print,timebin

names=['rot','not']

for itype=1,0,-1 do begin

   for zn=0,2 do begin
      
      models=dblarr(10,100,51) ;;;this is the array which will store model data as goes through
      
      xdum="abc"
      xdum2="abc"
      modelimf=0D
      modeltype=1
      openr,5,"/astro/store/gradscratch/tmp/tzdw/Geneva/input/input_geneva_"+znames(zn)+"_"+names(itype)+"_imf135_300"
      readf,5,totalimfnumber
      xdum5=dblarr(2)
      
      while(eof(5) eq 0) do begin  ;;;will load in model location from input geneva file and then go to load in that model
         readf,5,xdum2
         readf,5,xdum5
         modelimf=xdum5(0)
         modeltype=xdum5(1)
;0 - merger,  
;1 - normal primary, 
;2 - normal secondary,
;3 - single model,    
;4 - binary hmg to get black hole mass too  
         mixedage=0D
         mixedimf=0D
         
         dummy=fltarr(44)
         
         dum=strmid(xdum2,26,70)
         test1=FILE_TEST("/astro/store/gradscratch/tmp/tzdw/"+xdum2)
         print,xdum2,test1
         if(test1 eq 1) then begin
            openr,6,"/astro/store/gradscratch/tmp/tzdw/"+xdum2
;            print,dum
            lasttime=0e0
            lasttime2=mixedage
            count=-1
            itmax=0
            itmax2=0
            itlast=0
            dummy(1)=1e0
            it2last=round(10e0*alog10(mixedage))-60
            if(it2last le 0) then it2last=0
            if(it2last ge 50) then it2last=50
            initialmass=-100d0
            initialmasstest=0
            while(eof(6) eq 0 and dummy(1) le 100e9) do begin ;dummy(1) is Geneva Age
               readf,6,dummy

               count=count+1
               dt=dummy(1)-lasttime
               if(dt gt 0d0) then begin ;;;if actual timestep then go into loop
                  
                  it=round(10e0*alog10(dummy(1)))-60 ;;set up binning in time
                  if(it le 0) then it=0
                  if(it ge 50) then it=50
                  itmax=it
                  
                  xh=dummy(5) ;;sets up some parameters. dummy(5) is X
                  T=dummy(4) ;; dummy(4) is Temperature
                  T_sol = 5.778d3
                  logT_Tsol = T - alog10(T_sol)
                  logR = 0.5d0*(dummy(3) - 4d0*logT_Tsol) ;dummy(3) is luminosity
                  gravity=alog10(6.67259d-8*1.9891d33*dummy(2)/((10d0^logR)*6.9598d10)^2d0);dummy(2) is mass
                  C = dummy(7) + dummy(8) ;Carbon mass fractions
                  O = dummy(10) + dummy(11) + dummy(12) ;Nitrogen mass fractions
                  He = dummy(6) ;Helium mass fraction
                  COHE=(C/3.0+O/4.0)/He
                  logOfG= 3.676d+0 * T - 13.253d+0 ; using g_evol
                  
                  startype=0 ;;;;decide star type
                  
                  ;;; Using Geneva-code criteria (Georgey et al. 2013, Section 5.4)
                  ; T_WR = 4.0
                  ; X_WR = 0.3
                  ; X_WN/WC = 1e-5
                  ; Use BPASS WN/WC criteria
                  ; There are no K/M, just RSG, but all T_RSG < 3.66 are RSG anyway so we'll stick with that
                  ; T_YSG < 3.8
                  ; No A or B, use BPASS
                  ; T_O = 4.5
                  if(xh le 0.3 and T ge 4.0d0) then begin ;WR STAR 
                     if(xh gt 1e-5) then begin
                        startype=8 ;for WNH
                     endif else begin
                        startype=9                     ;!WN
                        if(COHE gt 0.03) then startype=10 ;!for WC
                     endelse
                  endif else begin                                        ;pre-WR star
                     startype=7                                           ;!M rsg
                     if(T ge 3.550) then startype=6                       ;!K rsg
                     if(T ge 3.660) then startype=5                       ;!YSG
                     if(T ge 3.8) then startype=4                         ;!A bsg
                     if(T ge 4.041) then startype=3                       ;!B bsg
                     if(T ge 4.5) then startype=1                        ;!O bsg
                     if(gravity lt logOfG and T ge 4.519d+0) then startype=2 ; !Of star
                  endelse
                  
                  
                  
                  x0=round(10e0*((dummy(3)))) ;log(L/Lsun)
                  
                  x1=startype-1
                  
                  
                  
                  if(x0 ge 0 and x0 le 99 and x1 ge 0 and x1 le 9) then begin
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
                        it2last=it2
                        lasttime2=dummy(1)+mixedage                        
                     endif                    
                     lasttime=dummy(1)
                     itlast=it
                   endif
               endif
            endwhile
            
            
            close,6
         endif
         
      endwhile
      close,5

      for n=0,50 do begin
         models(0:9,0:99,n)=models(0:9,0:99,n)/timebin(n)
      endfor
      
      openw,5,"trevorfiles_"+"_"+znames(zn)+"_"+names(itype)+".dat"
      outputarray=dblarr(10)
      for n=0,50 do begin
         for m=0,9 do begin
            outputarray(m)=total(models(m,49:99,n))
         endfor
         printf,5,6.0+0.1*n,outputarray,FORMAT='(100E16.7)'
      endfor
      close,5
      
      
      
      
      SAVE, /VARIABLES, FILENAME = 'output-data-'+names(itype)+'-'+znames(zn)+'.idl.dat'
      
   endfor
endfor
end

