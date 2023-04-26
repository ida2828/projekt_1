import numpy as np
from math import sin, cos, sqrt, atan, atan2, degrees, radians

def dms(x):
    sig = ' '
    if x < 0:
        sig = '-'
        x = abs(x)
    x = x * 180/np.pi
    d = int(x)
    m = int((x - d) * 60)
    s = (x - d - m/60) * 3600
    print(sig, "%3d %2d %7.5f" %(d, abs(m), abs(s)))

class Transformacje:
    def __init__(self, model: str = 'wgs84'):
        
        if model == "wgs84":
            self.a = 6378137.0 
            self.b = 6356752.31424518
        elif model == "grs80":
            self.a = 6378137.0
            self.b = 6356752.31414036
        elif model == "Krasowskiego":
            self.a = 6378245.0
            self.b = 6357397.155
        else:
            raise NotImplementedError(f"{model} model niezdefiniowany")
        self.flattening = (self.a - self.b)/self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
    
    def transform_XYZ2BLH(self,X,Y,Z, output = 'dec_degree'):
        p = np.sqrt(X**2 + Y**2)
        f = np.arctan(Z/(p*(1 - self.ecc2)))
        while True:
            N = self.a/ np.sqrt(1-self.ecc2 * np.sin(f)**2)
            h = (p / np.cos(f)) - N
            fp = f
            f = np.arctan(Z/(p*(1 - self.ecc2 * (N / (N + h)))))
            if np.abs(fp - f) < (0.000001/206265):
                break
        l = np.arctan2(Y,X)
        if output == "dec_degree": #stopnie dziesiętne
            return degrees(f), degrees(l), h 
        elif output == "dms": #stopnie minuty sekundy
            f = self.dms(degrees(f))
            l = self.dms(degrees(l))
            return f"{f[0]:02d}:{f[1]:02d}:{f[2]:.2f}", f"{l[0]:02d}:{l[1]:02d}:{l[2]:.2f}", f"{h:.3f}"
        else:
            raise NotImplementedError(f"{output} - format niezdefiniowany")
            
    
    def transform_BLH2XYZ (self, f, l, h,output = "XYZ"):
        N = self.a/ np.sqrt(1-self.ecc2 * np.sin(f)**2)
        X=(N+h)*np.cos(f)*np.cos(l)
        Y=(N+h)*np.cos(f)*np.sin(l)
        Z=(N*(1-self.ecc2)+h)*np.sin(f)
        if output == "XYZ":
            return(X,Y,Z)
        else:
            raise NotImplementedError(f"{output} - format niezdefiniowany")
            
    
     
    def transform_XYZ2neu(self, xa, ya, za, xb, yb, zb, phi, lam, h, output = "NEU,phi_stopnie,lam_stopnie"):
        N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
        X0 = (N + h) * np.cos(phi) * np.cos(lam)
        Y0 = (N + h) * np.cos(phi) * np.sin(lam)
        Z0 = ((1 - self.ecc2) * N + h) * np.sin(phi)
        
            # obliczenie współrzędnych XYZ
        dxyz = np.array([xb,yb,zb]) - np.array([xa,ya,za])
        X, Y, Z = dxyz - np.array([X0, Y0, Z0])
        
            # Macierz obrotu
        sin_phi = np.sin(phi)
        cos_phi = np.cos(phi)
        sin_lam = np.sin(lam)
        cos_lam = np.cos(lam)
        
        R = np.array([[-sin_lam,           cos_lam,           0],
                          [-sin_phi*cos_lam, -sin_phi*sin_lam,  cos_phi],
                          [ cos_phi*cos_lam,  cos_phi*sin_lam,  sin_phi]])
        
            # Obliczenie współrzędnych NEU
        NEU = np.dot(R, np.array([X, Y, Z]))
        
         
        N = NEU[0]
        E = NEU[1]
        U = -NEU[2]
        
        
        phi_stopnie = np.degrees(phi)
        lam_stopnie = np.degrees(lam)
        if output == "NEU,phi_stopnie,lam_stopnie":
            return N,E,U,phi_stopnie,lam_stopnie
        else:
            raise NotImplementedError(f"{output} - format niezdefiniowany")
        
        
    def transform_u2000(self, f, l, l0, s, output = "x00, y00, xgk, ygk"):
        m = 0.999923
        N=self.a/np.sqrt(1-self.ecc2*(np.sin(f))**2)
        t = np.tan(f)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(f))**2
       
        l= np.deg2rad(l)
        l0 =np.deg2rad(l0)
        d_l = l - l0

        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
    

        sig = self.a * ((A0*f) - (A2*np.sin(2*f)) + (A4*np.sin(4*f)) - (A6*np.sin(6*f)))
    
        xgk = sig + ((d_l**2)/2) * N *np.sin(f) * np.cos(f) * (1 + ((d_l**2)/12) * ((np.cos(f))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((d_l**4)/360) * ((np.cos(f))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        ygk = d_l * (N*np.cos(f)) * (1 + ((((d_l**2)/6) * (np.cos(f))**2) * (1-t**2+n2)) +  (((d_l**4)/(120)) * (np.cos(f)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
    
        x00 =m * xgk
        y00 =m * ygk + (s*1000000) + 500000
        if output == "x00, y00, xgk, ygk""": 
            return(x00, y00, xgk, ygk)
        else:
          raise NotImplementedError(f"{output} - format niezdefiniowany")  
    
        
    
    def transform_u1992(self,f, l, output = "x92, y92, xgk, ygk"):
        m = 0.9993
        N = self.a/(np.sqrt(1-self.ecc2 * np.sin(f)**2))
        t = np.tan(f)
        e_2 = self.ecc2/(1-self.ecc2)
        n2 = e_2 * (np.cos(f))**2
    
        l0 = np.deg2rad(19)
        d_l = l - l0
    
        A0 = 1 - (self.ecc2/4) - ((3*(self.ecc2**2))/64) - ((5*(self.ecc2**3))/256)   
        A2 = (3/8) * (self.ecc2 + ((self.ecc2**2)/4) + ((15 * (self.ecc2**3))/128))
        A4 = (15/256) * (self.ecc2**2 + ((3*(self.ecc2**3))/4))
        A6 = (35 * (self.ecc2**3))/3072 
    
        sigma = self.a * ((A0*f) - (A2*np.sin(2*f)) + (A4*np.sin(4*f)) - (A6*np.sin(6*f)))
    
        xgk = sigma + ((d_l**2)/2) * N *np.sin(f) * np.cos(f) * (1 + ((d_l**2)/12) * ((np.cos(f))**2) * (5 - t**2 + 9*n2 + 4*(n2**2)) + ((d_l**4)/360) * ((np.cos(f))**4) * (61 - (58*(t**2)) + (t**4) + (270*n2) - (330 * n2 *(t**2))))
        ygk = d_l * (N*np.cos(f)) * (1 + ((((d_l**2)/6) * (np.cos(f))**2) * (1-t**2+n2)) +  (((d_l**4)/(120)) * (np.cos(f)**4)) * (5 - (18 * (t**2)) + (t**4) + (14*n2) - (58*n2*(t**2))))
        
        x92 = m*xgk - 5300000
        y92 = m*ygk + 500000
        if output == "x92, y92, xgk, ygk":
            return(x92, y92, xgk, ygk)
        else:
          raise NotImplementedError(f"{output} - format niezdefiniowany")  


if __name__=='__main_':
    
    import argparse
    #XYZ2BLH
    parser = argparse.ArgumentParser(description='Transformacje współrzędnych')
    parser.add_argument('X', type=float, help='Współrzędna x')
    parser.add_argument('Y', type=float, help='Współrzędna y')
    parser.add_argument('Z', type=float, help='Współrzędna z')
    args = parser.parse_args()
    #BLH2XYZ
    parser.add_argument('f', type=float, help='Wartosc fi')
    parser.add_argument('l', type=float, help='Wartosc lambda')
    parser.add_argument('h', type=float, help='Wartosc h')
    args = parser.parse_args()
    #XYZ2NEU
    parser.add_argument('xa', type=float, help='X punktu A')
    parser.add_argument('ya', type=float, help='Y punktu A')
    parser.add_argument('za', type=float, help='Z punktu A')
    parser.add_argument('xb', type=float, help='X punktu B')
    parser.add_argument('yb', type=float, help='Y punktu B')
    parser.add_argument('zb', type=float, help='Z punktu B')
    parser.add_argument('phi', type=float, help='Wartosc fi')
    parser.add_argument('lam', type=float, help='Wartosc lambda')
    parser.add_argument('h', type=float, help='Wartosc h')
    args = parser.parse_args()
    #fl2u2000
    parser.add_argument('fi', type=float, help='Wartosc fi')
    parser.add_argument('lambda', type=float, help='Wartosc lambda')
    parser.add_argument('l0', type=float, help='Numer pasa')
    parser.add_argument('s', type=float, help='Numer odpowiadający numerowi pasa')
    args = parser.parse_args()
    #fl2u1992
    parser.add_argument('fi', type=float, help='Wartosc fi')
    parser.add_argument('lambda', type=float, help='Wartosc lambda')
    
    with open('plikprzykladowedane.py', 'r') as f:
        lines = f.readlines()

        X = []
        Y = []
        Z = []
        
        FI=[]
        LAM=[]
        H=[]
        
        N=[]
        E=[]
        U=[]
        
        X2000=[]
        Y2000=[]
        
        
        X1992=[]
        Y1992=[]
        
        
        for line in lines:
            if "print" in line:
                continue
            elif "," in line:
               rozdzielone_wsp=line.split(',')
               X.append(rozdzielone_wsp[0])
               Y.append(rozdzielone_wsp[1])
               Z.append(rozdzielone_wsp[2])
            elif "-" in line:   
               rozdzielone_wsp=line.split('-')
               FI.append(rozdzielone_wsp[0])
               LAM.append(rozdzielone_wsp[1])
               H.append(rozdzielone_wsp[2])
            elif ";" in line:
               rozdzielone_wsp=line.split(';')
               N.append(rozdzielone_wsp[0])
               E.append(rozdzielone_wsp[1])
               U.append(rozdzielone_wsp[2])
            elif "--" in line:
               rozdzielone_wsp=line.split('--')
               X2000.append(rozdzielone_wsp[0])
               Y2000.append(rozdzielone_wsp[1])
            else:
               rozdzielone_wsp=line.split(' ')
               X1992.append(rozdzielone_wsp[0])
               Y1992.append(rozdzielone_wsp[1])
        
        for (x,y,z) in zip(X,Y,Z):
            transformator_wgs84 = Transformacje(model = "wgs84")
            f, l, h = transformator_wgs84.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))
            
            transformator_grs80 = Transformacje(model = "grs80")
            f, l, h = transformator_grs80.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))

            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            f, l, h = transformator_Krasowskiego.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            X,Y,Z = transformator_wgs84.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_grs80 = Transformacje(model = "grs80")
            X,Y,Z = transformator_grs80.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            X,Y,Z = transformator_Krasowskiego.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            N, E, U, phi_stopnie, lam_stopnie = transformator_wgs84.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_grs80 = Transformacje(model = "grs80")
            N, E, U, phi_stopnie, lam_stopnie = transformator_grs80.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            N, E, U, phi_stopnie, lam_stopnie = transformator_Krasowskiego.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            x00, y00,xgk,ygk = transformator_wgs84.transform_u2000(f, l, l0, s)
            print(x00, y00,xgk,ygk)
            
            transformator_grs80 = Transformacje(model = "grs80")
            x00, y00,xgk,ygk = transformator_grs80.transform_u2000(f, l, l0,s)
            print(x00, y00,xgk,ygk)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            x00, y00,xgk,ygk = transformator_Krasowskiego.transform_u2000(f, l, l0,s)
            print(x00, y00,xgk,ygk)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            x92, y92,xgk,ygk = transformator_wgs84.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
            
            transformator_grs80 = Transformacje(model = "grs80")
            x92, y92,xgk,ygk = transformator_grs80.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")

            x92, y92,xgk,ygk = transformator_Krasowskiego.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
             
      
       
        
         







