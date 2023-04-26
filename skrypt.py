import numpy as np

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
            raise NotImplementedError(f"{model} model not implemented")
        self.flattening = (self.a - self.b)/self.a
        self.ecc2 = 2 * self.flattening - self.flattening ** 2
    
    def transform_XYZ2BLH(self,X,Y,Z):
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
        return(f,l,h)
    
    def transform_BLH2XYZ (self, f, l, h):
            N = self.a/ np.sqrt(1-self.ecc2 * np.sin(f)**2)
            X=(N+h)*np.cos(f)*np.cos(l)
            Y=(N+h)*np.cos(f)*np.sin(l)
            Z=(N*(1-self.ecc2)+h)*np.sin(f)
            return(X, Y, Z)
    
     
    def transform_XYZ2neu(self, xa, ya, za, xb, yb, zb, phi, lam, h):
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
        
        return (N, E, U, phi_stopnie, lam_stopnie)
        

    def transform_u2000(self, f, l, l0, s):
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
        return(x00, y00,xgk,ygk) 
    
    def transform_u1992(self,f, l):
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
        return(x92, y92, xgk, ygk)


if __name__=='__main_':
    
    import argparse
    parser = argparse.ArgumentParser(description='Transformacje współrzędnych')
    parser.add_argument('X', type=float, help='Współrzędna x')
    parser.add_argument('Y', type=float, help='Współrzędna y')
    parser.add_argument('Z', type=float, help='Współrzędna z')
    args = parser.parse_args()
    
    parser.add_argument('f', type=float, help='Wartosc fi')
    parser.add_argument('l', type=float, help='Wartosc lambda')
    parser.add_argument('h', type=float, help='Wartosc h')
    args = parser.parse_args()
    
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
    
    parser.add_argument('fi', type=float, help='Wartosc fi')
    parser.add_argument('lambda', type=float, help='Wartosc lambda')
    parser.add_argument('h', type=float, help='Wartosc h')
    parser.add_argument('l0', type=float, help='Numer pasa')
    parser.add_argument('s', type=float, help='Numer odpowiadający numerowi pasa')
    args = parser.parse_args()
    
   


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
            else:
                rozdzielone_wsp=line.split(',')
                X.append(rozdzielone_wsp[0])
                Y.append(rozdzielone_wsp[1])
                Z.append(rozdzielone_wsp[2])
        
        for (x,y,z) in zip(X,Y,Z):
            transformator_wgs84 = Transformacje(model = "wgs84")
            X = 3664940.500; Y = 1409153.590; Z = 5009571.170
            f, l, h = transformator_wgs84.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))
            
            transformator_grs80 = Transformacje(model = "grs80")
            X = 3664940.500; Y = 1409153.590; Z = 5009571.170
            f, l, h = transformator_grs80.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))

            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            X = 3664940.500; Y = 1409153.590; Z = 5009571.170
            f, l, h = transformator_Krasowskiego.transform_XYZ2BLH(X,Y,Z)
            print(f, l, h)
            (dms(f))
            (dms(l))
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            f = 3664940.500; l = 1409153.590; h = 5009571.170
            X,Y,Z = transformator_wgs84.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_grs80 = Transformacje(model = "grs80")
            f = 3664940.500; l = 1409153.590; h = 5009571.170
            X,Y,Z = transformator_grs80.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            f = 3664940.500; l = 1409153.590; h = 5009571.170
            X,Y,Z = transformator_Krasowskiego.transform_BLH2XYZ(f,l,h)
            print(X,Y,Z)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            xa = 15445; ya = 1544; za = 45541; xb = 451; yb = 154; zb = 154; phi = 17; lam= 179; h = 100
            N, E, U, phi_stopnie, lam_stopnie = transformator_wgs84.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_grs80 = Transformacje(model = "grs80")
            xa = 1; ya = 1; za = 1; xb = 1; yb = 1; zb = 1; phi = 1; lam= 1; h = 1
            N, E, U, phi_stopnie, lam_stopnie = transformator_grs80.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            xa = 1; ya = 1; za = 1; xb = 1; yb = 1; zb = 1; phi = 1; lam= 1; h = 1
            N, E, U, phi_stopnie, lam_stopnie = transformator_Krasowskiego.transform_XYZ2neu(xa, ya, za, xb, yb, zb, phi, lam, h)
            print(N, E, U, phi_stopnie, lam_stopnie)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            f = 9; l = 16; l0= 15; s =5
            x00, y00,xgk,ygk = transformator_wgs84.transform_u2000(f, l, l0, s)
            print(x00, y00,xgk,ygk)
            
            transformator_grs80 = Transformacje(model = "grs80")
            f = 9; l = 16; l0= 15; s=5
            x00, y00,xgk,ygk = transformator_grs80.transform_u2000(f, l, l0,s )
            print(x00, y00,xgk,ygk)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            f = 9; l = 14; l0= 15; s=5
            x00, y00,xgk,ygk = transformator_Krasowskiego.transform_u2000(f, l, l0,s)
            print(x00, y00,xgk,ygk)
            
            transformator_wgs84 = Transformacje(model = "wgs84")
            f = 9; l = 14
            x92, y92,xgk,ygk = transformator_wgs84.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
            
            transformator_grs80 = Transformacje(model = "grs80")
            f = 9; l = 14
            x92, y92,xgk,ygk = transformator_grs80.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
            
            transformator_Krasowskiego = Transformacje(model = "Krasowskiego")
            f = 9; l = 14
            x92, y92,xgk,ygk = transformator_Krasowskiego.transform_u1992(f, l)
            print(x92, y92, xgk, ygk)
      
       
        
         







