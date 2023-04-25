# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:09:13 2023

@author: 48531
"""
from math import *
import numpy as np
import argparse


#plik = 'belble.txt'
#p = open('belble.txt', "r")
#dane = p.readlines()

#for i in dane:
#    if "X" in i:
#        X = float(i[3:15])
#print('X = ', X)

#for i in dane:
#    if "Y" in i:
#        Y = float(i[20:32])
#print('Y = ', Y)

#for i in dane:
#    if "Z" in i:
#        Z = float(i[38:49])
#print('Z = ', Z)


class Transformations:
    def __init__(self, model: str = "WGS84"):
       # self.X = X
        #self.Y = Y
        #self.Z = Z
        """
        Parametry elipsoid:
        a - duża półos elipsoidy 
        e - mimosród elipsoidy
        + WGS84:
            + GRS80:
                + Elipsoida Krasowskiego:
        """
        if model == "WGS84":
            self.a = 6378137.000
            self.e2 = 0.00669438002290
        elif model == "GRS80":
            self.a = 6378137.000
            self.e2 = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
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
    
    def dms2deg(self, dms):
        '''
        This method convert deg, min, sek to decimal degree
        '''
        d = dms[0]; m = dms[1]; s = dms[2]
        decimal_degree = d+m/60+s/3600
        return (decimal_degree)
    
    def deg2dms(self, decimal_degree):
        '''
        This method convert decimal degree to deg, min, sek   
        '''
        decimal_degree = decimal_degree / 3600
        dms = np.array([])
        st = np.floor(decimal_degree)
        append(dms, st)
        m = np.floor((decimal_degree - st)*60)
        append(dms, m)
        sek = (decimal_degree - st - m/60)*3600
        append(dms, sek)
        print(dms)
        return (dms)
    
        
        
        
        # XYZ ---> BLH - ALGORYTM HIRVONENA
        """
            Następujący algorytm przelicza współrzędne z układu ortokartezjańskiego na współrzędne geodezyjne.
        """
    def hirvonen(self, X, Y, Z, output: 'dms'):
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
            Wyznaczenie azymutu i odleglosci miedzy punktami
        """
    def vincenty(self, fiA, lamA, fiB, lamB):
        b = self.a * np.sqrt(1 - self.e2)
        f = 1 - b / self.a
        UA = np.arctan((1 - f) * np.tan(fiA))
        UB = np.arctan((1 - f) * np.tan(fiB))
        dl = lamB - lamA
        
        L = dl
        while True:
            sinsigma = np.sqrt(((np.cos(UB) * np.sin(L))**2) + ((np.cos(UA) * np.sin(UB) - np.sin(UA) * np.cos(UB) * np.cos(L))**2))
            cossigma = np.sin(UA) * np.sin(UB) + np.cos(UA) * np.cos(UB) * np.cos(L)
            sigma = np.arctan(sinsigma / cossigma)
            sinAz = (np.cos(UA) * np.cos(UB) * np.sin(L)) / sinsigma
            cos2Az = 1 - (sinAz)**2
            cos2sigma = cossigma - ((2 * np.sin(UA) * np.sin(UB)) / cos2Az)
            C = (f / 16) * cos2Az * (4 + f * (4 - 3 * cos2Az))
            Lp = L
            L = dl + (1 - C) * f * sinAz * (sigma + C * sinsigma * (cos2sigma + C * cossigma * (-1 + 2 * (cos2sigma)**2)))
            if abs(L - Lp) < (0.000001/206265):
                break
                
            u2 = ((self.a**2 - b**2) / b**2) * cos2Az
            A = 1 + (u2/16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
            B = (u2/1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
            dsigma = B * sinsigma * (cos2sigma + 1/4 * B * (cossigma * (-1 + 2 * (cos2sigma)**2) - 1/6 * B * cos2sigma * (-3 + 4 * (sinsigma)**2) * (-3 + 4 * (cos2sigma)**2)))
            sAB = b * A * (sigma - dsigma)
            AzAB = np.arctan2((np.cos(UB) * np.sin(L)), ((np.cos(UA) * np.sin(UB)) - (np.sin(UA) * np.cos(UB) * np.cos(L))))
            AzBA = np.arctan2((np.cos(UA) * np.sin(L)), ((-np.sin(UA) * np.cos(UB)) + (np.cos(UA) * np.sin(UB) * np.cos(L)))) + np.pi
            if AzAB > 2 * np.pi:
                AzAB = AzAB - 2 * np.pi
            if AzAB < 0:  
                AzAB = AzAB + 2 * np.pi
            if AzBA > 2 * np.pi:
                AzBA = AzBA - 2 * np.pi
            if AzBA < 0:
                AzBA = AzBA + 2 * np.pi
                    
        return(sAB, AzAB, AzBA)

        """
            Wyznaczenie kata zenitalnego
        """
    def zenitalny(self, fiA, lamA, hA, fiB, lamB, hB):
        N_a = self.a/np.sqrt(1-self.e2*np.sin(fiA)**2)
        R_a = np.array([(N_a+hA)*np.cos(fiA)*np.cos(lamA),
                (N_a+hA)*np.cos(fiA)*np.sin(lamA), 
                (N_a*(1-self.e2)+hA)*np.sin(fiA)])
            
        N_b = self.a/np.sqrt(1-self.e2*np.sin(fiB)**2)
        R_b = np.array([(N_b+hB)*np.cos(fiB)*np.cos(lamB),
                       (N_b+hB)*np.cos(fiB)*np.sin(lamB),
                       (N_b*(1-self.e2)+hB)*np.sin(fiB)])
            
        R_delta = R_b - R_a
        R_dlugosc  = np.linalg.norm(R_delta, 2)
            
        R = R_delta / R_dlugosc
        
        n  = np.array([-1 * np.sin(fiA)*np.cos(lamA),
                       -1 * np.sin(fiA)*np.sin(lamA),
                       np.cos(fiA)])
            
        e  = np.array([-1 * np.sin(lamA),
                         np.cos(lamA),
                         0]) 
            
        u  = np.array([np.cos(fiA)*np.cos(lamA),
                  np.cos(fiA)*np.sin(lamA),
                  np.sin(fiA)])
            
        cos_z = float(np.dot(u, R))
        z = np.arccos(cos_z)
            
        return(z)

    """
            Wyznaczenie macierzy dneu
        """
    def saz2neu(sAB, AzAB, z):
        dneu = np.array([sAB * np.sin(z) * np.cos(AzAB),
                             sAB * np.sin(z) * np.sin(AzAB),
                             sAB * np.cos(z)])
        return(dneu)



        # TRANSFORMACJA WSP BL ---> 1992
        """
            Algorytm przelicza współrzędne geodezyjne (BL) na współrzędne w układzie 1992 (XY)
        """
    def cale92(self,fi, lam):
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

if __name__ == "__main__":
    # Tworzenie parsera argumentów
    
    parser = argparse.ArgumentParser(description="Przykładowy program z użyciem argparse")
    parser.add_argument("--model", type=str, default="WGS84", choices=["WGS84", "GRS80", "Krasowski"],
                        help="Model elipsoidy (domyślnie: WGS84)")
    parser.add_argument("--X", type=float, help="Wartość X")
    parser.add_argument("--Y", type=float, help="Wartość Y")
    parser.add_argument("--Z", type=float, help="Wartość Z")
    parser.add_argument("--output", type=str, default="dec_degree", choices=["dec_degree", "dms"],
                        help="Format wyników (domyślnie: dec_degree)")
    args = parser.parse_args()

    # Tworzenie instancji klasy Transformacje na podstawie podanych argumentów
    transformations = Transformations(args.model)

    # Wywołanie odpowiednich funkcji na podstawie przekazanych argumentów
    wynik = transformations.hirvonen(args.X, args.Y, args.Z, args.output)


