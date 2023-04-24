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
            self_1 = 6378137.000
            self_2 = 0.00669438002290
        elif model == "GRS80":
            self_1 = 6378137.000
            self_2 = 6356752.31414036
        elif model == "Elipsoida Krasowskiego":
            self_1 = 6378245.000
            self_2 = 0.00669342162296
    #Pomocnicze funkcje
        """
        Poniższe funkcje są funkcjami pomocniczymi dla obliczeń transformacji
        """
        def Npu(fi, a ,e2):     #promien krzywizny w I wertykale
            N = a / np.sqrt(1 - e2 * np.sin(fi)**2)
        return(N)

        def Sigma(fi, a, e2):
            A0 = 1 - (e2/4) - (3*(e2)**2)/64 -  (5*(e2)**3)/256
            A2 = 3/8 * (e2 + (e2)**2/4 + 15*(e2)**3/128)
            A4 = 15/256 * ( (e2)**2 + (3*((e2)**3))/4 )
            A6 = 35 * (e2)**3 / 3072
            sigma = a * ( A0 * fi - A2 * np.sin(2*fi) + A4 * np.sin(4*fi) - A6 * np.sin(6*fi) )
           
            return(sigma)
        
        # XYZ ---> BLH - ALGORYTM HIRVONENA

        def hirvonen(X, Y, Z, a, e2):
            p = np.sqrt(X**2 + Y**2)
            fi = np.arctan(Z / (p * (1 - e2)))
            while True:
                N = Npu(fi, a, e2)
                h = p / np.cos(fi) - N
                fip = fi     #fip - fi poprzednie, fi - fi nowe
                fi = np.arctan(Z / (p * (1 - N * e2 / (N + h))))
                if abs(fip - fi) < (0.000001/206265):
                    break
            l = np.arctan2(Y, X)
            return(fi, l, h)




        # BLH ---> XYZ

        def filh2XYZ(fi, l, h, a, e2):
            while True:
                N = Npu(fi, a, e2)
                X = (N + h) * np.cos(fi) * np.cos(l)
                Xp = X
                Y = (N + h) * np.cos(fi) * np.sin(l)
                Z = (N * (1 - e2) + h) * np.sin(fi)
                if abs(Xp - X) < (0.000001/206265):
                    break
            return(X, Y, Z)




        # XYZ ---> NEU
        
        
        #Transformacja współrzędnych BL -> 2000
            """
            Następujący algorytm umożliwia przeliczenie współrzędnych geodezyjnych na współrzędne ortokartezjańskie w układzie 2000
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
             
             b2 = (a**2) * (1-e2)   #krotsza polos
             
             e2p = ( a**2 - b2 ) / b2   #drugi mimosrod elipsy
             
             dlam = lam - lam0
             
             t = np.tan(fi)
             
             ni = np.sqrt(e2p * (np.cos(fi))**2)
             
             N = Npu(fi, a, e2)
             
             sigma = Sigma(fi, a, e2)

             
             
             xgk = sigma + ((dlam**2)/2)*N*np.sin(fi)*np.cos(fi) * ( 1+ ((dlam**2)/12)*(np.cos(fi))**2 * ( 5 - (t**2)+9*(ni**2) + 4*(ni**4)     )  + ((dlam**4)/360)*(np.cos(fi)**4) * (61-58*(t**2)+(t**4) + 270*(ni**2) - 330*(ni**2)*(t**2))  )
             ygk = (dlam*N* np.cos(fi)) * (1+(((dlam)**2/6)*(np.cos(fi))**2) *(1-(t**2)+(ni**2))+((dlam**4)/120)*(np.cos(fi)**4)*(5-18*(t**2)+(t**4)+14*(ni**2)-58*(ni**2)*(t**2)) )
             
             x00 = xgk * m
             y00 = ygk * m + strefa*1000000 + 500000
             
             return(x00, y00)   
         