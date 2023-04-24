# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 11:09:13 2023

@author: 48531
"""
from math import *
import numpy as np

class Transformations:
    def __init__(self, model: str = "WGS84"):
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
    #Pomocnicze funkcje
        """
        Poniższe funkcje są funkcjami pomocniczymi dla obliczeń transformacji
        """
        def Npu(self, fi, output: 'dec_degree'):     #promien krzywizny w I wertykale
            N = self.a / np.sqrt(1 - self.e2 * np.sin(fi)**2)
        return(N)

        def Sigma(self, fi, output: 'dec_degree'):
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
        def hirvonen(self, X, Y, Z, output: 'dec_degree'):
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
        def filh2XYZ(self, fi, lam, h, output: 'dec_degree'):
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
            Przeliczenie wspolrzednych puktow A i B z ukladu kartezjanskiego (XYZ) na geodezyjny (fi lam h)
            """
        def hirvonen(X, Y, Z, a, e2):
            p = np.sqrt(X**2 + Y**2)
            fi = np.arctan(Z / (p * (1 - e2)))
            while True:
                N = Np(fi, a, e2)
                h = p / np.cos(fi) - N
                fip = fi     #fip - fi poprzednie, fi - fi nowe
                fi = np.arctan(Z / (p * (1 - N * e2 / (N + h))))
                if abs(fip - fi) < (0.000001/206265):
                    break
            l = np.arctan2(Y, X)
            return(fi, l, h)

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
        def cale00(self, fi, lam, output: 'dec_degree'):
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
         