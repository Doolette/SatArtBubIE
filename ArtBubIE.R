#isobaric dissolution of bubbles in blood
#inner ear dissolved gas kinetics
#References
#R Core Team (2021). R: A language and environment for statistical computing. R Foundation for
#  Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
#  version 4.05 (2021-03-21)
#ArtBub functions based on Srinivasan RS, Gerth WA, and Powell MR. Mathematical models of diffusion-limited gas bubble dynamics in tissue. J Appl Physiol 86: 732-741, 1999.
#variable shell thickness based on Tikuisis P, Ward CA, and Venter RD. Bubble evolution in a stirred volume of liquid closed to mass transport. J Appl Phys 54: 9, 1983.
#InnerEarFF function based on Doolette DJ, and Mitchell SJ. Biophysical basis for inner ear decompression sickness. J Appl Physiol 94: 2145-2150, 2003.

library(deSolve)
#arterial bubble models, isobric ambient pressure
#constant shell thickness
Isobaric.ArtBub.3RWT.CT<-function(t, state, parameters){
  with(as.list(c(state, parameters)),{
  Pi<-Pamb-(PaO2+PaCO2+PH2O)+2*ST/r
  dr<-(sol*D*(1/h+1/r)*(Pt-Pi))/(Pamb-(PaO2+PaCO2+PH2O)+4*ST/(3*r))
  list(c(dr),Pi=Pi)})} #end with(as.list..., end function

  #variable shell thickness
Isobaric.ArtBub.3RWT.VT<-function(t, state, parameters){
  with(as.list(c(state, parameters)),{
  Pi<-Pamb-(PaO2+PaCO2+PH2O)+2*ST/r
  h=Z*(1-exp(-r/Z))
  dr<-(sol*D*(1/h+1/r)*(Pt-Pi))/(Pamb-(PaO2+PaCO2+PH2O)+4*ST/(3*r))
  list(c(dr),Pi=Pi, h=h)})} #end with(as.list..., end function

#examples
times<-seq(0,100,0.01)
state10<-c(r=10)
parametersAirSurf<-c(sol=0.0151,D=2900,ST=0.55268, Pamb=1, Pt=0.741178, PaO2=0.14672, PaCO2=0.0503, PH2O=0.0618)#D (micron^2/s);h (micron); ST (atm micron)
parametersHeliox250msw<-c(sol=0.0094,D=9300, ST=0.55268, Pamb=26, Pt=25.5382, PaO2=0.3490, PaCO2=0.0503, PH2O=0.0618)
#parameter pressures in atm

outAirSurfCT10<-ode(y=state10, times=times, func=Isobaric.ArtBub.3RWT.CT, parms=c(parametersAirSurf,h=1))
tail(outAirSurfCT10)
plot(outAirSurfCT10)
#bubble lifetime = 0.765s (h=1)

outAirSurfVT10<-ode(y=state10, times=times, func=Isobaric.ArtBub.3RWT.VT, parms=c(parametersAirSurf,Z=50))
tail(outAirSurfVT10)
plot(outAirSurfVT10)

outHeliox250mswCT10<-ode(y=state10, times=times, func=Isobaric.ArtBub.3RWT.CT, parms=c(parametersHeliox250msw,h=1))
tail(outHeliox250mswCT10)
plot(outHeliox250mswCT10)

outHeliox250mswVT10<-ode(y=state10, times=times, func=Isobaric.ArtBub.3RWT.VT, parms=c(parametersHeliox250msw,Z=50))
tail(outHeliox250mswVT10)
plot(outHeliox250mswVT10)


#Inner Ear Function using forcing functions for Part, must be named PartH.fun and PartN.fun
InnerEarFF<-function(t, yini, parameters){
  with(as.list(c(yini, parameters)),{
  #membraneous labyrinth
  dPvenN <- (Q*SbN*(PartN.fun(t) - PvenN) - PSpN*(PvenN - PperN) - PSeN*(PvenN - PendN))/(StN*Vven)
  dPvenH <- (Q*SbH*(PartH.fun(t) - PvenH) - PSpH*(PvenH - PperH) - PSeH*(PvenH - PendH))/(StH*Vven)
  #perilymph
  dPperN <- (PSpN*(PvenN - PperN))/(SbN*Vper)#  - PSrwN*(PperN - PmeN)
  dPperH <- (PSpH*(PvenH - PperH))/(SbH*Vper)#  - PSrwH*(PperH - PmeH)
  #endolymph
  dPendN <- (PSeN*(PvenN - PendN))/(SbN*Vend)
  dPendH <- (PSeH*(PvenH - PendH))/(SbH*Vend)
  list(c(dPvenH,dPperH,dPendH,dPvenN,dPperN,dPendN),Pven=PvenN+PvenH,Pper=PperN+PperH,Pend=PendN+PendH)
  })} #end with(as.list..., end function

# example forcing functions for PartH, PartN
t.s<-cumsum(c(0,4.9,0.1,165/60,35)*60)#convert minutes to s and cumsum
Pamb<-c(31.24,31.24,31.24,26.25,26.25)
FHe<-c(0.9871959,0.9871959,0.9847619,0.9847619,0.9847619)
PartH<- FHe*(Pamb-0.0618)-0.0503#PH2O=0.0618, PaCO2=0.0503
PartH.fun<-approxfun(x=t.s, y=PartH, method="linear", rule=2)
PartN.fun<-approxfun(x=c(0,2565), y=c(0,0), method="linear", rule=2)#PrtN=0

times<-seq(0,2565,1)
ie.fixed<-c(Q=3.6E-4, StN=0.015, SbN=0.015, StH=0.01, SbH=0.01, PSeN=1.3E-5*0.015*3/0.025, PSpN=1.3E-5*0.015*3/0.033, PSeH=3.94E-5*0.01*3/0.025, PSpH=3.94E-5*0.01*3/0.033, Vven=0.070, Vend=0.038, Vper=0.166)
yini<-c(PvenH=30.72869,PperH=30.72869,PendH=30.72869,PvenN=0,PperN=0,PendN=0)
Upward.ode<-ode(y=yini, times=times, func=InnerEarFF, parms=ie.fixed)
