import numpy as np



#a class to describe TolTEC cluster sensitvity
class TolTEC:

    #instantiate
    def __init__(self,atmFactor):
        self.atmFactor = atmFactor
        self.FOV_arcmin2 = np.pi*2.**2
        self.ynorm = 0.5e-4
        self.vnorm = 500.      #km/s
        self.taunorm = 0.005
        self.MS1p1 = 13.
    
    #time to map given depth in mJy at 270GHz (1.1mm)
    def time_S_hours(self, depth_mJy, area_deg2):
        return area_deg2/(self.MS1p1/self.atmFactor)/(depth_mJy**2)

    #time to map given y with 5sigma/pixel at 150GHz (2mm)
    def time_tSZ_mins(self,y,area_arcmin2,sigmaPerPix):
        return 138.*self.atmFactor/7.*\
            area_arcmin2/(2.*self.FOV_arcmin2)*\
            (self.ynorm/y)**2*\
            (sigmaPerPix/5.)**2

    #time to map given v with 2sigma/pixel at 220GHz (1.4mm)
    def time_kSZ_hours(self,v,area_arcmin2,taue,sigmaPerPix):
        return 115.*self.atmFactor/7.*\
            area_arcmin2/(2.*self.FOV_arcmin2)*\
            (self.taunorm/taue)**2 *\
            (self.vnorm/v)**2 *\
            (sigmaPerPix/2.)**2
    



        
