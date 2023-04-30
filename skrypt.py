import numpy as np
from math import sin, cos, sqrt, atan, atan2, degrees, radians
import argparse


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

print('test')

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
    
    def XYZ2BLH(self,X,Y,Z): 
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
    
    
    def BLH2XYZ (self, f, l, h,output = "XYZ"):
        N = self.a/ np.sqrt(1-self.ecc2 * np.sin(f)**2)
        X=(N+h)*np.cos(f)*np.cos(l)
        Y=(N+h)*np.cos(f)*np.sin(l)
        Z=(N*(1-self.ecc2)+h)*np.sin(f)
        if output == "XYZ":
            return(X,Y,Z)
        else:
            raise NotImplementedError(f"{output} - format niezdefiniowany")
            
    
     
    def XYZ2neu(self, xa, ya, za, xb, yb, zb, phi, lam, h, output = "NEU,phi_stopnie,lam_stopnie"):
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
       
    def u2000(self, phi, lam, h):
            if lam>13.5 and lam <16.5:
                s = 5
                l0 = 15.00000
            elif lam>16.5 and lam <19.5:
                s = 6
                l0 = 18.00000
            elif lam>19.5 and lam <22.5:
                s = 7
                l0 = 21.00000
            elif lam>22.5 and lam <25.5:
                s = 8
                l0 = 24.00000
            
            l0=np.deg2rad(l0)
            phi=np.deg2rad(phi)
            lam=np.deg2rad(lam)
            
            a2 = self.a**2
            b2 = a2 * (1 - self.ecc2)
            e_2 = (a2 - b2)/b2
            dl = lam - l0
            dl2 = dl**2
            dl4 = dl**4
            t = np.tan(phi)
            t2 = t**2
            t4 = t**4
            n2 = e_2 * (np.cos(phi)**2)
            n4 = n2 ** 2
            N = self.a / np.sqrt(1 - self.ecc2 * np.sin(phi)**2)
            e4 = self.ecc2**2
            e6 = self.ecc2**3
            A0 = 1 - (self.ecc2/4) - ((3*e4)/64) - ((5*e6)/256)
            A2 = (3/8) * (self.ecc2 + e4/4 + (15*e6)/128)
            A4 = (15/256) * (e4 + (3*e6)/4)
            A6 = (35*e6)/3072
            sigma = self.a * ((A0 * phi) - A2 * np.sin(2*phi) + A4 * np.sin(4*phi) - A6 * np.sin(6*phi))
            xgk = sigma + ((dl**2)/2) * N * np.sin(phi) * np.cos(phi) * (1 + ((dl**2)/12)*(np.cos(phi)**2)*(5 - t2 + 9 * n2 + 4 * n4) + (dl4/360) * (np.cos(phi)**4)*(61 - (58 * t2) + t4 + (270 * n2) - (330 * n2 * t2)))
            ygk = dl * N * np.cos(phi) * (1 + (dl2/6) * (np.cos(phi)**2) * (1 - t2 + n2) + (dl4/120) * (np.cos(phi)**4) * (5 - (18 * t2) + t4 + (14 * n2) - 58 * n2 * t2))
            
            x00 = xgk * 0.999923
            y00 = ygk * 0.999923 + s * 1000000 + 500000
            return(x00, y00)
        
    
    def u1992(self,f, l, output = "x92, y92, xgk, ygk"):
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
            return(x92, y92)
        else:
          raise NotImplementedError(f"{output} - format niezdefiniowany")  


if __name__=='__main__':
    

    parser = argparse.ArgumentParser(description='Transformacje współrzędnych')
    
    parser.add_argument('-f', type=str, help='Scieżka do pliku z danymi do transformacji')
    parser.add_argument('-t', type=str, help='rodzaj transformacji')
    parser.add_argument('-m', type=str, help='model elipsoidy', choices=['wgs84', 'grs80', 'Krasowskiego'])
    parser.add_argument('-fk', type=str, help='Ścieżka do pliku do którego program zapisuje dane')
    args = parser.parse_args() 
    print(args.f, args.t, args.m, args.fk)
    
    
#XYZ2BLH   
    if args.t == 'XYZ2BLH':
        if args.m == 'wgs84':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                

                X=[]
                Y=[]
                Z=[]
                
                for line in lines:
                        rozdzielone_wsp=line.split(',')
                        X.append(float(rozdzielone_wsp[0]))
                        Y.append(float(rozdzielone_wsp[1]))
                        Z.append(float(rozdzielone_wsp[2]))
        
                for (x,y,z) in zip(X,Y,Z):
                    XYZ_wgs84 =Transformacje(model = "wgs84")
                    f,l,h = XYZ_wgs84.XYZ2BLH(float(x), float(y), float(z))
                    wyniki.write("{}, {}, {}\n".format(f, l, h))
               
        elif args.m =='grs80':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                X=[]
                Y=[]
                Z=[]
                for line in lines:
                        rozdzielone_wsp=line.split(',')
                        X.append(rozdzielone_wsp[0])
                        Y.append(rozdzielone_wsp[1])
                        Z.append(rozdzielone_wsp[2])
                
                for (x,y,z) in zip(X,Y,Z):
                    XYZ_grs80=Transformacje(model = "grs80")
                    f,l,h= XYZ_grs80.XYZ2BLH(float(x), float(y), float(z))
                            
                    wyniki.write("{}, {}, {}\n".format(f, l, h))
                        
        elif args.m == 'Krasowskiego':
            with open(args.f, 'r') as file:
                    lines = file.readlines()
                    wyniki = open(args.fk, 'w')
                    
                    X=[]
                    Y=[]
                    Z=[]
                    for line in lines:
                            rozdzielone_wsp=line.split(',')
                            X.append(rozdzielone_wsp[0])
                            Y.append(rozdzielone_wsp[1])
                            Z.append(rozdzielone_wsp[2])
            
                    for (x,y,z) in zip(X,Y,Z):
                        #utworzy obiekt o tych współrzędnych
                        XYZ_Krasowskiego=Transformacje(model = "Krasowskiego")
                        f,l,h = XYZ_Krasowskiego.XYZ2BLH(float(x), float(y), float(z))
                        wyniki.write("{}, {}, {}\n".format(f, l, h))
                
                    
                    
#BLH2XYZ
    elif args.t == 'BLH2XYZ':
            if args.m == 'wgs84':
                with open(args.f, 'r') as file:
                    lines = file.readlines()
                    wyniki = open(args.fk, 'w')
                    
                    B=[]
                    L=[]
                    H=[]
                    for line in lines:

                            rozdzielone_wsp = line.split(',')
                            B.append(rozdzielone_wsp[0])
                            L.append(rozdzielone_wsp[1])
                            H.append(rozdzielone_wsp[2])
                            
                            
                    for (b,l,h) in zip(B,L,H):
                        BLH_wgs84=Transformacje(model = "wgs84")
                        X,Y,Z = BLH_wgs84.BLH2XYZ(float(b), float(l), float(h))
                        wyniki.write("{}, {}, {}\n".format(X,Y,Z))
            
            elif args.m =='grs80':
                with open(args.f, 'r') as file:
                    lines = file.readlines()
                    wyniki = open(args.fk, 'w')
                    
                    B=[]
                    L=[]
                    H=[]
                    
                    for line in lines:
                    
                            rozdzielone_wsp = line.split(',')
                            B.append(rozdzielone_wsp[0])
                            L.append(rozdzielone_wsp[1])
                            H.append(rozdzielone_wsp[2])
                           
                            
                    for (b,l,h) in zip(B,L,H):
                        BLH_grs80=Transformacje(model = "grs80")
                        X,Y,Z = BLH_grs80.BLH2XYZ(float(b), float(l), float(h))
                        wyniki.write("{}, {}, {}\n".format(X,Y,Z))
                        
            elif args.m == 'Krasowskiego':
                with open(args.f, 'r') as file:
                        lines = file.readlines()
                        wyniki = open(args.fk, 'w')
                        
                        B=[]
                        L=[]
                        H=[]
                        for line in lines:
                            
                            rozdzielone_wsp = line.split(',')
                            B.append(rozdzielone_wsp[0])
                            L.append(rozdzielone_wsp[1])
                            H.append(rozdzielone_wsp[2])
                            
                                
                                
                        for (b,l,h) in zip(B,L,H):
                            BLH_Krasowskiego=Transformacje(model = "Krasowskiego")
                            X,Y,Z = BLH_Krasowskiego.BLH2XYZ(float(b), float(l), float(h))
                            wyniki.write("{}, {}, {}\n".format(X,Y,Z))
                    
#XYZ2neu                
    elif args.t == 'XYZ2neu':
        if args.m == 'wgs84':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                xa=[]
                ya=[]
                za=[]
                xb=[]
                yb=[]
                zb=[]
                phi=[]
                lam=[]
                h=[]
                
                
                for line in lines:
                    
                        rozdzielone_wsp=line.split(',')
                        xa.append(float(rozdzielone_wsp[0]))
                        ya.append(float(rozdzielone_wsp[1]))
                        za.append(float(rozdzielone_wsp[2]))
                        xb.append(float(rozdzielone_wsp[3]))
                        yb.append(float(rozdzielone_wsp[4]))
                        zb.append(float(rozdzielone_wsp[5]))
                        phi.append(float(rozdzielone_wsp[6]))
                        lam.append(float(rozdzielone_wsp[7]))
                        h.append(float(rozdzielone_wsp[8]))
        
                for (XA,YA,ZA, XB, YB, ZB, PHI, LAM, H) in zip(xa,ya,za,xb,yb,zb,phi,lam,h):
                    neu_wgs84 =Transformacje(model = "wgs84")
                    N, E, U, phi_stopnie, lam_stopnie = neu_wgs84.XYZ2neu(float(XA), float(YA), float(ZA), float(XB), float(YB), float(ZB),float(PHI), float(LAM), float(H))
                    wyniki.write("{}, {}, {}\n".format(N, E, U, phi_stopnie, lam_stopnie))
                    
                
        elif args.m =='grs80':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                xa=[]
                ya=[]
                za=[]
                xb=[]
                yb=[]
                zb=[]
                phi=[]
                lam=[]
                h=[]
                
                
                for line in lines:
                    
                        rozdzielone_wsp=line.split(',')
                        xa.append(float(rozdzielone_wsp[0]))
                        ya.append(float(rozdzielone_wsp[1]))
                        za.append(float(rozdzielone_wsp[2]))
                        xb.append(float(rozdzielone_wsp[3]))
                        yb.append(float(rozdzielone_wsp[4]))
                        zb.append(float(rozdzielone_wsp[5]))
                        phi.append(float(rozdzielone_wsp[6]))
                        lam.append(float(rozdzielone_wsp[7]))
                        h.append(float(rozdzielone_wsp[8]))
        
                for (XA,YA,ZA, XB, YB, ZB, PHI, LAM, H) in zip(xa,ya,za,xb,yb,zb,phi,lam,h):
                    neu_grs80 =Transformacje(model = "grs80")
                    N, E, U, phi_stopnie, lam_stopnie = neu_grs80.XYZ2neu(float(XA), float(YA), float(ZA), float(XB), float(YB), float(ZB),float(PHI), float(LAM), float(H))
                    wyniki.write("{}, {}, {}\n".format(N, E, U, phi_stopnie, lam_stopnie))
                
        elif args.m == 'Krasowskiego':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                xa=[]
                ya=[]
                za=[]
                xb=[]
                yb=[]
                zb=[]
                phi=[]
                lam=[]
                h=[]
                
                
                for line in lines:
                    
                        rozdzielone_wsp=line.split(',')
                        xa.append(float(rozdzielone_wsp[0]))
                        ya.append(float(rozdzielone_wsp[1]))
                        za.append(float(rozdzielone_wsp[2]))
                        xb.append(float(rozdzielone_wsp[3]))
                        yb.append(float(rozdzielone_wsp[4]))
                        zb.append(float(rozdzielone_wsp[5]))
                        phi.append(float(rozdzielone_wsp[6]))
                        lam.append(float(rozdzielone_wsp[7]))
                        h.append(float(rozdzielone_wsp[8]))
        
                for (XA,YA,ZA, XB, YB, ZB, PHI, LAM, H) in zip(xa,ya,za,xb,yb,zb,phi,lam,h):
                    neu_Krasowskiego =Transformacje(model = "Krasowskiego")
                    N, E, U, phi_stopnie, lam_stopnie = neu_Krasowskiego.XYZ2neu(float(XA), float(YA), float(ZA), float(XB), float(YB), float(ZB),float(PHI), float(LAM), float(H))
                    wyniki.write("{}, {}, {}\n".format(N, E, U, phi_stopnie, lam_stopnie))

#flh2 2000
    elif args.t == 'u2000':
        if args.m == 'wgs84':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                phi=[]
                lam=[]
                h=[]
                
                for line in lines:
                    
                        rozdzielone_wsp=line.split(',')
                        phi.append(float(rozdzielone_wsp[0]))
                        lam.append(float(rozdzielone_wsp[1]))
                        h.append(float(rozdzielone_wsp[2]))
        
                for (PHI, LAM, H) in zip(phi,lam,h):
                    flh_wgs84 =Transformacje(model = "wgs84")
                    x00,y00 = flh_wgs84.u2000(float(PHI), float(LAM), float(H))
                    wyniki.write("{}, {}\n".format(x00, y00))
                    
                
                
        elif args.m =='grs80':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                phi=[]
                lam=[]
                h=[]
                
                for line in lines:
                   
                        rozdzielone_wsp=line.split(',')
                        phi.append(float(rozdzielone_wsp[0]))
                        lam.append(float(rozdzielone_wsp[1]))
                        h.append(float(rozdzielone_wsp[2]))
                
                for (PHI, LAM, H) in zip(phi, lam, h):
                    flh_grs80=Transformacje(model = "grs80")
                    x00, y00= flh_grs80.u2000(float(PHI), float(LAM), float(H))
                    wyniki.write("{}, {}\n".format(x00, y00))
                
        elif args.m == 'Krasowskiego':
            with open(args.f, 'r') as file:
                    lines = file.readlines()
                    wyniki = open(args.fk, 'w')
                    
                    phi=[]
                    lam=[]
                    h=[]
                    for line in lines:
                        
                            rozdzielone_wsp=line.split(',')
                            phi.append(rozdzielone_wsp[0])
                            lam.append(rozdzielone_wsp[1])
                            h.append(rozdzielone_wsp[2])
            
                    for (PHI,LAM,H) in zip(phi, lam, h):
                        #utworzy obiekt o tych współrzędnych
                        flh_Krasowskiego=Transformacje(model = "Krasowskiego")
                        x00, y00 = flh_Krasowskiego.u2000(float(PHI), float(LAM), float(H))
                        wyniki.write("{}, {}\n".format(x00, y00))
                        
#flh2 1992
    elif args.t =='u1992':
        if args.m =='wgs84':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
             
                B=[]
                L=[]
                
                for line in lines:
                    
                    rozdzielone_wsp = line.split(',')
                    B.append(rozdzielone_wsp[0])
                    L.append(rozdzielone_wsp[1])
    
                for (b,l) in zip(B,L):
                        BLH_wgs84=Transformacje(model = "wgs84")
                        X1992, Y1992 = BLH_wgs84.u1992(float(b), float(l))
                        wyniki.write("{}, {}\n".format(X1992, Y1992))
                 
                    
        elif args.m =='grs80':
            with open(args.f, 'r') as file:
                 lines = file.readlines()
                 wyniki = open(args.fk, 'w')
                 
                 B=[]
                 L=[]
                 
                 for line in lines:
                     
                     rozdzielone_wsp = line.split(',')
                     B.append(rozdzielone_wsp[0])
                     L.append(rozdzielone_wsp[1])
                     
                         
                 for (b,l) in zip(B,L):
                     #utworzy obiekt o tych współrzędnych
                     BLH_grs80=Transformacje(model = "grs80")
                     X1992, Y1992 = BLH_grs80.u1992(float(b), float(l))
                     wyniki.write("{}, {}\n".format(X1992, Y1992))
                     
                 
        elif args.m =='Krasowski':
            with open(args.f, 'r') as file:
                lines = file.readlines()
                wyniki = open(args.fk, 'w')
                
                B=[]
                L=[]
                
                for line in lines:
                    
                    rozdzielone_wsp = line.split(',')
                    B.append(rozdzielone_wsp[0])
                    L.append(rozdzielone_wsp[1])
                   
                        
                for (b,l) in zip(B,L):
                    #utworzy obiekt o tych współrzędnych
                    BLH_Krasowski=Transformacje(model = "Krasowski")
                    X1992, Y1992 = BLH_Krasowski.u1992(float(b), float(l))
                    wyniki.write("{},{}\n".format(X1992, Y1992))


       
        
        






