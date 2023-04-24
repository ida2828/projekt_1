import numpy as np



class Transformacje:
    def __init__(self, model: str = 'wgs84'):
        
        if model == "wgs84":
            self.a = 6378137.0 
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        else:
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b)/self.a
        self.e2 = 2 * self.flattening - self.flattening ** 2
    
    def transform_XYZ2BLH(self,X,Y,Z):
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z/(p*(1 - self.e2)))
        while True:
            N = self.a/ np.sqrt(1-self.e2 * np.sin(f)**2)
            h = (p / np.cos(f)) - N
            fp = f
            f = np.arctan(Z/(p*(1 - self.e2 * (N / (N + h)))))
            if np.abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)
        return(f,l,h)
    
    def transform_flh2XYZ (self, f, l, h):
            N = self.a/ np.sqrt(1-self.e2 * np.sin(f)**2)
            X=(N+h)*np.cos(f)*np.cos(l)
            Y=(N+h)*np.cos(f)*np.sin(l)
            Z=(N*(1-self.e2)+h)*np.sin(f)
            return(X, Y, Z)
    
    
transformator_wgs84 = Transformacje("wgs84")
print(transformator_wgs84.transform_XYZ2BLH(3782550,1084730,5002940))

transformator_grs80 = Transformacje("grs80")
print(transformator_grs80.transform_XYZ2BLH(3782550,1084730,5002940))






