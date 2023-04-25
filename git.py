# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:09:13 2023

@author: 48531
"""
from math import *
import numpy as np
import argparse




class Transformations:
    def __init__(self, model):
       # self.X = X
        #self.Y = Y
        #self.Z = Z
        """
        Parametry elipsoid:
        a - duża półos elipsoidy 
        e2 - kwadrat mimosrodu elipsoidy
            + WGS84:
            + GRS80:
            + Elipsoida Krasowskiego:
        """
        if model == 'WGS84':
            self.a = 6378137.000
            self.e2 = 0.00669438002290
        elif model == 'GRS80':
            self.a = 6378137.000
            self.e2 = 0.00669438002290
        elif model == 'Elipsoida Krasowskiego':
            self.a = 6378245.000
            self.e2 = 0.00669342162296
        else:
            raise AttributeError("Jesli napisałes małymi literami, spróbuj napisać wielkimi")
            #Pomocnicze funkcje
            
            
            
        """
        Poniższe funkcje są funkcjami pomocniczymi dla obliczeń transformacji
        """
    def Npu(self, fi):     #promien krzywizny w I wertykale
        N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
        return(N)

    def Sigma(self, fi):
        A0 = 1 - (self.e2/4) - (3*(self.e2)**2)/64 -  (5*(self.e2)**3)/256
        A2 = 3/8 * (self.e2 + (self.e2)**2/4 + 15*(self.e2)**3/128)
        A4 = 15/256 * ( (self.e2)**2 + (3*((self.e2)**3))/4 )
        A6 = 35 * (self.e2)**3 / 3072
        sigma = self.a * ( A0 * fi - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi) )
           
        return(sigma)
    

        
        # XYZ ---> BLH - ALGORYTM HIRVONENA
        """
            Następujący algorytm przelicza współrzędne z układu ortokartezjańskiego na współrzędne geodezyjne.
        """
    def hirvonen(self, X, Y, Z):
        p = np.sqrt(X**2 + Y**2)
        fi = np.arctan(Z / (p * (1 - self.e2)))
        while True:
            N = self.Npu(fi)
            h = p / np.cos(fi) - N
            fip = fi     #fip - fi poprzednie, fi - fi nowe
            fi = np.arctan(Z / (p * (1 - N * self.e2 / (N + h))))
            if abs(fip - fi) < (0.000001/206265):
                break

            lam = np.arctan2(Y, X)
            return(fi, lam, h)



        # BLH ---> XYZ
        """
            Algorytm przelicza współrzędne geodezyjne (BLH) na współrzędne w układzie ortokartezjańskim (XYZ)
        """
    def filh2XYZ(self, fi, lam, h):
        while True:
            N = self.Npu(fi)
            X = (N + h) * np.cos(fi) * np.cos(lam)
            Xp = X
            Y = (N + h) * np.cos(fi) * np.sin(lam)
            Z = (N * (1 - self.e2) + h) * np.sin(fi)
            if abs(Xp - X) < (0.000001/206265):
                break
        return(X, Y, Z)



        # XYZ ---> NEU
        """
            Obliczenie macierzy Rneu
        """
    def Rneu(self, fi, lam):
        Rneu = np.array([[-np.sin(fi)*np.cos(lam), -np.sin(lam), np.cos(fi)*np.cos(lam)],
                         [-np.sin(fi)*np.sin(lam),  np.cos(lam), np.cos(fi)*np.sin(lam)],
                         [             np.cos(fi),            0,             np.sin(fi)]])
        return Rneu
    
    
        """
            Przeliczenie wsp XYZ na neu
        """
    def xyz2neup(self, X, Y, Z, X0, Y0, Z0):
        neu = []
        p = np.sqrt(X0**2 + Y0**2)
        fi = np.arctan(Z0 / (p*(1 - self.e2)))
        while True:
            N = self.Npu(fi)
            h = (p / np.cos(fi)) - N
            fi_poprzednia = fi
            fi = np.arctan((Z0 / p)/(1-((N * self.e2)/(N + h))))
            if abs(fi_poprzednia - fi) < (0.000001/206265):
                break 
        N = self.Npu(fi)
        h = p/np.cos(fi) - N
        lam = np.arctan(Y0 / X0)
        
        R_neu = self.Rneu(fi, lam)
        for X, Y, Z in zip(X, Y, Z):
            X_sr = [X - X0, Y - Y0, Z - Z0] 
            X_rneu = R_neu.T@X_sr
            neu.append(X_rneu.T)
            
        return neu



        # TRANSFORMACJA WSP BL ---> 1992
        """
            Algorytm przelicza współrzędne geodezyjne (BL) na współrzędne w układzie 1992 (XY)
        """
    def cale92(self, fi, lam):
        lam0 = (19*np.pi)/180
                
        m = 0.9993
                
        b2 = (self.a**2) * (1-self.e2)   #krotsza polos
                
        e2p = (self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
                
        dlam = lam - lam0
                
        t = np.tan(fi)
                
        ni = np.sqrt(e2p * (np.cos(fi))**2)
                
        N = self.Npu(fi)
                
        sigma = self.Sigma(fi)

                
                
        xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
        ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
                
        x92 = xgk*m - 5300000
        y92 = ygk*m + 500000
                
        return(x92, y92)
            
            
            
        # TRANSFORMACJA WSP BL ---> 2000
        """
            Następujący algorytm umożliwia przeliczenie współrzędnych geodezyjnych (BLH) na współrzędne w układzie 2000 (XY)
        """
    def cale00(self, fi, lam):
        m=0.999923
        lam0=0 
        strefa = 0
        if lam >np.deg2rad(13.5) and lam < np.deg2rad(16.5):
            strefa = 5
            lam0 = np.deg2rad(15)
        elif lam >np.deg2rad(16.5) and lam < np.deg2rad(19.5):
            strefa = 6
            lam0 = np.deg2rad(18)
        elif lam >np.deg2rad(19.5) and lam < np.deg2rad(22.5):
            strefa =7
            lam0 = np.deg2rad(21)
        elif lam >np.deg2rad(22.5) and lam < np.deg2rad(25.5):
            strefa = 8
            lam0 = np.deg2rad(24)
        else:
            print("Punkt poza strefami odwzorowawczymi układu PL-2000")        
             
        b2 = (self.a**2) * (1-self.e2)   #krotsza polos
             
        e2p = ( self.a**2 - b2 ) / b2   #drugi mimosrod elipsy
             
        dlam = lam - lam0
             
        t = np.tan(fi)
             
        ni = np.sqrt(e2p * (np.cos(fi))**2)
             
        N = self.Npu(fi)
             
        sigma = self.Sigma(fi)

             
             
        xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
        ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
             
        x00 = xgk * m
        y00 = ygk * m + strefa*1000000 + 500000
             
        return(x00, y00)  
    
    
    
    def pliczek(self, plik, funkcja: type = str):
    
            
        if funkcja == "XYZ_BLH":
            X = []
            Y = []
            Z = []

            data = open("plik", 'r')
            lines = data.read().splitlines()
            for el in lines:
                for i in el:
                    x = el.split( )[0]
                    y = el.split( )[1]
                    z = el.split( )[2]
                X.append(x)
                Y.append(y)
                Z.append(z)
            b, l, h = self.hirvonen(X, Y, Z)
            np.savetxt(f"WYNIK_{funkcja}.txt", np.rad2deg(b), np.rad2deg(l), "%0.3f" %h,  delimiter = ";")
        
        
        elif funkcja == "BLH_XYZ":
            fi = []
            lam = []
            h = []

            data = open("plik", 'r')
            lines = data.read().splitlines()
            for el in lines:
                for i in el:
                    f = np.deg2rad(el.split( )[0])
                    l = np.deg2rad(el.split( )[1])
                    hi = el.split( )[2]
                fi.append(f)
                lam.append(l)
                h.append(hi)
            X, Y, Z = self.filh2XYZ(fi, lam, h)
            np.savetxt(f"WYNIK_{funkcja}.txt", "%0.3f" %X, "%0.3f" %Y, "%0.3f" %Z,  delimiter = ";")
            
        elif funkcja == "XYZ_NEU":
            X = []
            Y = []
            Z = []
            X0 = []
            Y0 = []
            Z0 = []

            data = open("plik", 'r')
            lines = data.read().splitlines()
            for el in lines:
                for i in el:
                    x = el[1:].split( )[0]
                    y = el[1:].split( )[1]
                    z = el[1:].split( )[2]
                    x0 = el[0].split( )[0]
                    y0 = el[0].split( )[1]
                    z0 = el[0].split( )[2]
                X.append(x)
                Y.append(y)
                Z.append(z)
                X0.append(x0)
                Y0.append(y0)
                Z0.append(z0)
                    
            neu = self.xyz2neup(X, Y, Z, X0, Y0, Z0)
            np.savetxt(f"WYNIK_{funkcja}.txt", neu)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Podaj plik")
    parser.add_argument("-pliczek", type = str, help = "Podaj nazwę pliku, w którym znajdują się dane wejsciowe (ps. oprócz nazwy podaj rozszerzenie:)")
    parser.add_argument("-elip", type = str, help = "Wybierz elipsoidę, na której ma wykonać się transformacja, wpisz jedną: 'WGS84', 'GRS80', 'Elipsoida Krasowskiego' ")
    parser.add_argument("-trans", type = str, help = "Wybierz transformację jaką chcesz obliczyć: 'XYZ_BLH', 'BLH_XYZ', 'XYZ_NEU' ")
    args = parser.parse_args()
    elip = {'WGS84' : 'WGS84', 'GRS80' : 'GRS80', 'Elipsoida Krasowskiego' : 'Elipsoida Krasowskiego'}
    trans = {'XYZ_BLH' : 'hirvonen', 'BLH_XYZ' : 'filh2XYZ', 'XYZ_NEU' : 'xyz2neup'}
            
            
            
            
#if __name__ == "__main__":

#    parser = argparse.ArgumentParser(description='Proszę wpisać tutaj współrzędne w zależnosci od tego, co chcesz policzyć: XYZ->BLH, lub BLH->XYZ; FL->układ 1992, lub FL->układ 2000)')
#    parser.add_argument('-a', '--argument', help = "Ten argument jest opcjonalny", required = False)
#    parser.add_argument('X', type = float, help = "Współrzędna X")
#    parser.add_argument('Y', type = float, help = "Współrzędna Y")
#    parser.add_argument('Z', type = float, help = "Współrzędna Z")
#    args = parser.parse_args()
#    geo = Transformations(model = "WGS84")
#    FI, LAM, H = geo.hirvonen(X, Y, Z, output = 'dec_degree')
# quit() #nwm co to dodało mi się to jak zaimportowałam biblioteke argparse



