
### Combined aerosol and close contact modules - base model
### Code updated 10.25.2021

#### Setup ####

####################################################################

rm(list=ls()) #clear all variables (this is good to have at the top of the script to clear any preexisting variables in your environment that might interfere)


# Open Packages
library(mvtnorm)
library(mc2d) 
library(dplyr)
set.seed(12345)

#10001 simulations - variability
ndvar(10001)


# Source functions --------------------------------------------------------

# Make sure it is in the same working directory as current code
# Or add file path before file name
source("C:\\Users\\jsoboli\\Dropbox\\Emory 2020\\Covid Grant\\SARS-CoV-2 Baseline Model Analysis Check for Manuscript 5.13.21\\base_model_functions.R")


# Set controls/interventions -------------------------------------------------------

#Master controls for both aerosol and close contact modules
m.Event <-"cough"
m.room.exchange <- mcstoc(runif, type="V",min=2, max =2)
m.Humidity<-"high"

#Intervention controls for both aerosol and close contact modules
i.Clean1<-0
i.Clean2<-0
i.Clean3<-0
i.Clean4<-0
i.Clean5<-0
i.Clean6<-0
i.Clean7<-0
i.Clean8<-0
i.Clean9<-0
i.Clean10<-0
i.Clean11<-0
i.Clean12<-0
i.Clean13<-0
i.Clean14<-0
i.Clean15<-0
i.Clean16<-0
i.index.mask<-"none"
i.susceptible.mask<-"none"
i.HW<-"no"
i.Glove<-"no"
i.surface.clean.eff <- mcstoc(runif, type = "V", min = 0.999, max = 0.999)


# Aerosol function call ---------------------------------------------------

aero.dose <- Aerofunc(Event = m.Event, room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8,Clean9=i.Clean9,Clean10=i.Clean10,Clean11=i.Clean11,Clean12=i.Clean12,Clean13=i.Clean13,Clean14=i.Clean14,Clean15=i.Clean15,Clean16=i.Clean16, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)


aero.dose.clean <- select(aero.dose,
                          starts_with("a"),
                          starts_with("f"))
#View(aero.dose.clean)



# Close contact function call ---------------------------------------------

dose50601m <- Dosefunc(Event = m.Event, Volume.Fraction="50-60", Distance=1.0, Vol.Frac.Dist.Name = "50601m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose601001m <- Dosefunc(Event = m.Event, Volume.Fraction="60-100", Distance=1.0, Vol.Frac.Dist.Name = "601001m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose1007501m <- Dosefunc(Event = m.Event, Volume.Fraction="100+", Distance=1.0, Vol.Frac.Dist.Name = "1007501m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose50602m <- Dosefunc(Event = m.Event, Volume.Fraction="50-60", Distance=2.0, Vol.Frac.Dist.Name = "50602m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose601002m <- Dosefunc(Event = m.Event, Volume.Fraction="60-100", Distance=2.0, Vol.Frac.Dist.Name = "601002m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)
dose50603m <- Dosefunc(Event = m.Event, Volume.Fraction="50-60", Distance=3.0, Vol.Frac.Dist.Name = "50603m", room.exchange = m.room.exchange, Clean1=i.Clean1, Clean2=i.Clean2, Clean3=i.Clean3, Clean4=i.Clean4, Clean5=i.Clean5, Clean6=i.Clean6, Clean7=i.Clean7, Clean8=i.Clean8, sc.eff.p = i.surface.clean.eff, Humidity= m.Humidity, Infected.Mask = i.index.mask, HW = i.HW, Glove = i.Glove, Susceptible.Mask = i.susceptible.mask)


dose5060_1m <- select(dose50601m,
                      starts_with("a"),
                      starts_with("f"))

dose60100_1m <- select(dose601001m,
                       starts_with("a"),
                       starts_with("f"))


dose100750_1m <- select(dose1007501m,
                        starts_with("a"),
                        starts_with("f"))


dose5060_2m <- select(dose50602m,
                      starts_with("a"),
                      starts_with("f"))


dose60100_2m <- select(dose601002m,
                       starts_with("a"),
                       starts_with("f"))


dose5060_3m <- select(dose50603m,
                      starts_with("a"),
                      starts_with("f"))


# Combine aerosol and close contact doses ---------------------------------

#adding aerosol and aerosol fomite to close contact doses
dose1m <-cbind (dose5060_1m, dose60100_1m, dose100750_1m, aero.dose.clean)
dose2m <-cbind (dose5060_2m, dose60100_2m, aero.dose.clean)
dose3m <-cbind (dose5060_3m, aero.dose.clean)


dose1m <- mutate( dose1m, 
                  a1m1h = a50601m1h + a601001m1h + a1007501m1h + aero1h,
                  a1m2h = a50601m2h + a601001m2h + a1007501m2h + aero2h,
                  a1m3h = a50601m3h + a601001m3h + a1007501m3h + aero3h,
                  a1m4h = a50601m4h + a601001m4h + a1007501m4h + aero4h,
                  a1m5h = a50601m5h + a601001m5h + a1007501m5h + aero5h,
                  a1m6h = a50601m6h + a601001m6h + a1007501m6h + aero6h,
                  a1m7h = a50601m7h + a601001m7h + a1007501m7h + aero7h,
                  a1m8h = a50601m8h + a601001m8h + a1007501m8h + aero8h,
                  f1m1h = f50601m1h + f601001m1h + f1007501m1h + f1h,
                  f1m2h = f50601m2h + f601001m2h + f1007501m2h + f2h,
                  f1m3h = f50601m3h + f601001m3h + f1007501m3h + f3h,
                  f1m4h = f50601m4h + f601001m4h + f1007501m4h + f4h,
                  f1m5h = f50601m5h + f601001m5h + f1007501m5h + f5h,
                  f1m6h = f50601m6h + f601001m6h + f1007501m6h + f6h,
                  f1m7h = f50601m7h + f601001m7h + f1007501m7h + f7h,
                  f1m8h = f50601m8h + f601001m8h + f1007501m8h + f8h,
                  af1m1h = af50601m1h + af601001m1h + af1007501m1h + aerof1h,
                  af1m2h = af50601m2h + af601001m2h + af1007501m2h + aerof2h,
                  af1m3h = af50601m3h + af601001m3h + af1007501m3h + aerof3h,
                  af1m4h = af50601m4h + af601001m4h + af1007501m4h + aerof4h,
                  af1m5h = af50601m5h + af601001m5h + af1007501m5h + aerof5h,
                  af1m6h = af50601m6h + af601001m6h + af1007501m6h + aerof6h,
                  af1m7h = af50601m7h + af601001m7h + af1007501m7h + aerof7h,
                  af1m8h = af50601m8h + af601001m8h + af1007501m8h + aerof8h)
dose1m_risk <- select(dose1m, a1m1h:af1m8h)

dose1m_contribution <-mutate(dose1m,
                  droplet1m1h = a50601m1h + a601001m1h + a1007501m1h,
                  droplet1m2h = a50601m2h + a601001m2h + a1007501m2h,
                  droplet1m3h = a50601m3h + a601001m3h + a1007501m3h,
                  droplet1m4h = a50601m4h + a601001m4h + a1007501m4h,
                  droplet1m5h = a50601m5h + a601001m5h + a1007501m5h,
                  droplet1m6h = a50601m6h + a601001m6h + a1007501m6h,
                  droplet1m7h = a50601m7h + a601001m7h + a1007501m7h,
                  droplet1m8h = a50601m8h + a601001m8h + a1007501m8h)

dose1m_contrib_adf <- select(dose1m_contribution, droplet1m1h:droplet1m8h, aero1h:aero8h, f1m1h:f1m8h, af1m1h:af1m8h)

dose2m <- mutate( dose2m, 
                  a2m1h = a50602m1h + a601002m1h + aero1h ,
                  a2m2h = a50602m2h + a601002m2h + aero2h,
                  a2m3h = a50602m3h + a601002m3h + aero3h,
                  a2m4h = a50602m4h + a601002m4h + aero4h,
                  a2m5h = a50602m5h + a601002m5h + aero5h,
                  a2m6h = a50602m6h + a601002m6h + aero6h,
                  a2m7h = a50602m7h + a601002m7h + aero7h,
                  a2m8h = a50602m8h + a601002m8h + aero8h,
                  f2m1h = f50602m1h + f601002m1h + f1h,
                  f2m2h = f50602m2h + f601002m2h + f2h,
                  f2m3h = f50602m3h + f601002m3h + f3h,
                  f2m4h = f50602m4h + f601002m4h + f4h,
                  f2m5h = f50602m5h + f601002m5h + f5h,
                  f2m6h = f50602m6h + f601002m6h + f6h,
                  f2m7h = f50602m7h + f601002m7h + f7h,
                  f2m8h = f50602m8h + f601002m8h + f8h,
                  af2m1h = af50602m1h + af601002m1h + aerof1h,
                  af2m2h = af50602m2h + af601002m2h + aerof2h,
                  af2m3h = af50602m3h + af601002m3h + aerof3h,
                  af2m4h = af50602m4h + af601002m4h + aerof4h,
                  af2m5h = af50602m5h + af601002m5h + aerof5h,
                  af2m6h = af50602m6h + af601002m6h + aerof6h,
                  af2m7h = af50602m7h + af601002m7h + aerof7h,
                  af2m8h = af50602m8h + af601002m8h + aerof8h )

dose2m_risk <- select(dose2m, a2m1h:af2m8h)

dose2m_contribution <-mutate(dose2m,
                             droplet2m1h = a50602m1h + a601002m1h,
                             droplet2m2h = a50602m2h + a601002m2h,
                             droplet2m3h = a50602m3h + a601002m3h,
                             droplet2m4h = a50602m4h + a601002m4h,
                             droplet2m5h = a50602m5h + a601002m5h,
                             droplet2m6h = a50602m6h + a601002m6h,
                             droplet2m7h = a50602m7h + a601002m7h,
                             droplet2m8h = a50602m8h + a601002m8h)

dose2m_contrib_adf <- select(dose2m_contribution, droplet2m1h:droplet2m8h, aero1h:aero8h, f2m1h:f2m8h, af2m1h:af2m8h)

dose3m <- mutate( dose3m, 
                  a3m1h = a50603m1h + aero1h,
                  a3m2h = a50603m2h + aero2h,
                  a3m3h = a50603m3h + aero3h,
                  a3m4h = a50603m4h + aero4h,
                  a3m5h = a50603m5h + aero5h,
                  a3m6h = a50603m6h + aero6h,
                  a3m7h = a50603m7h + aero7h,
                  a3m8h = a50603m8h + aero8h,
                  f3m1h = f50603m1h + f1h,
                  f3m2h = f50603m2h + f2h,
                  f3m3h = f50603m3h + f3h,
                  f3m4h = f50603m4h + f4h,
                  f3m5h = f50603m5h + f5h,
                  f3m6h = f50603m6h + f6h,
                  f3m7h = f50603m7h + f7h,
                  f3m8h = f50603m8h + f8h,
                  af3m1h = af50603m1h + aerof1h,
                  af3m2h = af50603m2h + aerof2h,
                  af3m3h = af50603m3h + aerof3h,
                  af3m4h = af50603m4h + aerof4h,
                  af3m5h = af50603m5h + aerof5h,
                  af3m6h = af50603m6h + aerof6h,
                  af3m7h = af50603m7h + aerof7h,
                  af3m8h = af50603m8h + aerof8h)

dose3m_risk <- select(dose3m, a3m1h:af3m8h)

dose3m_contribution <-mutate(dose3m,
                             droplet3m1h = a50603m1h,
                             droplet3m2h = a50603m2h,
                             droplet3m3h = a50603m3h,
                             droplet3m4h = a50603m4h,
                             droplet3m5h = a50603m5h,
                             droplet3m6h = a50603m6h,
                             droplet3m7h = a50603m7h,
                             droplet3m8h = a50603m8h)

dose3m_contrib_adf <- select(dose3m_contribution, droplet3m1h:droplet3m8h, aero1h:aero8h, f3m1h:f3m8h, af3m1h:af3m8h)

# Calculate risk ----------------------------------------------------------

#dose-response parameter updated to Julian et al 2020 preprint
krisk=0.00680

#vaccination for susceptible worker
fullvaccine <- mcstoc(runif, type="V",min=0.01, max =0.23)
partialvaccine <- mcstoc(runif, type="V",min=0.26, max =0.48)
fullvaccine <- unmc (fullvaccine,drop = TRUE)
fullvaccine <- as.data.frame(fullvaccine)
partialvaccine<-unmc (partialvaccine,drop = TRUE)
partialvaccine<-as.data.frame(partialvaccine)

# Pull combined doses through dose response for data frame output of risk

#aerosol module risk (>3m)
riskaero.df = 1-exp(-krisk*aero.dose.clean)
riskaero.quant<-as.data.frame(t(apply(riskaero.df, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
riskaero.mean<- as.data.frame(colMeans(riskaero.df))
riskaero.comb.stats<- cbind(riskaero.quant, riskaero.mean)
#View(riskaero.comb.stats)
write.csv(riskaero.comb.stats, "C:\\Users\\jsoboli\\Desktop\\aerorisk.csv", row.names=TRUE)

#close contact 1m risk
risk1m.df = 1-exp(-krisk*dose1m_risk)
risk1mvaxx <-cbind(risk1m.df, fullvaccine, partialvaccine)
risk1mvaxxfull <-mutate(risk1mvaxx,
                        
                      fullvaxxaf1m1h = fullvaccine * af1m1h,
                      fullvaxxaf1m2h = fullvaccine * af1m2h,
                      fullvaxxaf1m3h = fullvaccine * af1m3h,
                      fullvaxxaf1m4h = fullvaccine * af1m4h,
                      fullvaxxaf1m5h = fullvaccine * af1m5h,
                      fullvaxxaf1m6h = fullvaccine * af1m6h,
                      fullvaxxaf1m7h = fullvaccine * af1m7h,
                      fullvaxxaf1m8h = fullvaccine * af1m8h,
                      partialvaxxaf1m1h = partialvaccine * af1m1h,
                      partialvaxxaf1m2h = partialvaccine * af1m2h,
                      partialvaxxaf1m3h = partialvaccine * af1m3h,
                      partialvaxxaf1m4h = partialvaccine * af1m4h,
                      partialvaxxaf1m5h = partialvaccine * af1m5h,
                      partialvaxxaf1m6h = partialvaccine * af1m6h,
                      partialvaxxaf1m7h = partialvaccine * af1m7h,
                      partialvaxxaf1m8h = partialvaccine * af1m8h)

risk1m.quant<-as.data.frame(t(apply(risk1mvaxxfull, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
risk1m.mean<- as.data.frame(colMeans(risk1mvaxxfull))
risk1m.comb.stats<- cbind(risk1m.quant, risk1m.mean)
#View(risk1m.comb.stats)
write.csv(risk1m.comb.stats, "C:\\Users\\jsoboli\\Desktop\\risk1m.csv", row.names=TRUE)


#close contact 1m dose contribution
dose1m_contrib_adf.quant <-as.data.frame(t(apply(dose1m_contrib_adf, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
#View(dose1m_contrib_adf.quant)
write.csv(dose1m_contrib_adf.quant, "C:\\Users\\jsoboli\\Desktop\\dose1mcontrib.csv", row.names=TRUE)

#close contact 2m risk
risk2m.df = 1-exp(-krisk*dose2m_risk)
risk2mvaxx <-cbind(risk2m.df, fullvaccine, partialvaccine)
risk2mvaxxfull <-mutate(risk2mvaxx,
                        
                        fullvaxxaf2m1h = fullvaccine * af2m1h,
                        fullvaxxaf2m2h = fullvaccine * af2m2h,
                        fullvaxxaf2m3h = fullvaccine * af2m3h,
                        fullvaxxaf2m4h = fullvaccine * af2m4h,
                        fullvaxxaf2m5h = fullvaccine * af2m5h,
                        fullvaxxaf2m6h = fullvaccine * af2m6h,
                        fullvaxxaf2m7h = fullvaccine * af2m7h,
                        fullvaxxaf2m8h = fullvaccine * af2m8h,
                        partialvaxxaf2m1h = partialvaccine * af2m1h,
                        partialvaxxaf2m2h = partialvaccine * af2m2h,
                        partialvaxxaf2m3h = partialvaccine * af2m3h,
                        partialvaxxaf2m4h = partialvaccine * af2m4h,
                        partialvaxxaf2m5h = partialvaccine * af2m5h,
                        partialvaxxaf2m6h = partialvaccine * af2m6h,
                        partialvaxxaf2m7h = partialvaccine * af2m7h,
                        partialvaxxaf2m8h = partialvaccine * af2m8h)

risk2m.quant<-as.data.frame(t(apply(risk2mvaxxfull, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
risk2m.mean<- as.data.frame(colMeans(risk2mvaxxfull))
risk2m.comb.stats<- cbind(risk2m.quant, risk2m.mean)
write.csv(risk2m.comb.stats, "C:\\Users\\jsoboli\\Desktop\\risk2m.csv", row.names=TRUE)

#close contact 2m dose contribution
dose2m_contrib_adf.quant <-as.data.frame(t(apply(dose2m_contrib_adf, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
#View(dose2m_contrib_adf.quant)
write.csv(dose2m_contrib_adf.quant, "C:\\Users\\jsoboli\\Desktop\\dose2mcontrib.csv", row.names=TRUE)

#close contact 3m risk
risk3m.df = 1-exp(-krisk*dose3m_risk)
risk3mvaxx <-cbind(risk3m.df, fullvaccine, partialvaccine)
risk3mvaxxfull <-mutate(risk3mvaxx,
                        
                        fullvaxxaf3m1h = fullvaccine * af3m1h,
                        fullvaxxaf3m2h = fullvaccine * af3m2h,
                        fullvaxxaf3m3h = fullvaccine * af3m3h,
                        fullvaxxaf3m4h = fullvaccine * af3m4h,
                        fullvaxxaf3m5h = fullvaccine * af3m5h,
                        fullvaxxaf3m6h = fullvaccine * af3m6h,
                        fullvaxxaf3m7h = fullvaccine * af3m7h,
                        fullvaxxaf3m8h = fullvaccine * af3m8h,
                        partialvaxxaf3m1h = partialvaccine * af3m1h,
                        partialvaxxaf3m2h = partialvaccine * af3m2h,
                        partialvaxxaf3m3h = partialvaccine * af3m3h,
                        partialvaxxaf3m4h = partialvaccine * af3m4h,
                        partialvaxxaf3m5h = partialvaccine * af3m5h,
                        partialvaxxaf3m6h = partialvaccine * af3m6h,
                        partialvaxxaf3m7h = partialvaccine * af3m7h,
                        partialvaxxaf3m8h = partialvaccine * af3m8h)

risk3m.quant<-as.data.frame(t(apply(risk3mvaxxfull, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
risk3m.mean<- as.data.frame(colMeans(risk3mvaxxfull))
risk3m.comb.stats<- cbind(risk3m.quant, risk3m.mean)
write.csv(risk3m.comb.stats, "C:\\Users\\jsoboli\\Desktop\\risk3m.csv", row.names=TRUE)

#close contact 3m dose contribution
dose3m_contrib_adf.quant <-as.data.frame(t(apply(dose3m_contrib_adf, 2, quantile, probs=c(0.05, 0.25, 0.5, 0.75, 0.95))))
#View(dose3m_contrib_adf.quant)
write.csv(dose3m_contrib_adf.quant, "C:\\Users\\jsoboli\\Desktop\\dose3mcontrib.csv", row.names=TRUE)











