PROJEKT NR 1: TRANSFORMACJE


Dostępne transformacje:
```sh
   XYZ ---> BLH
   BLH ---> XYZ
   XYZ ---> NEU
   BL ---> PL1992
   BL ---> PL2000
```  
 
 Dostępne elipsoidy:
 ```sh
   GRS80
   WGS84
   Elipsoida Krasowskiego
 ```
 
 Wymagania:
 ```sh
   python 3.9 lub python 3.8.9
   biblioteka numpy
   biblioteka argparse
 ``` 
  
 Opis programu:
 
 Plik przyjmuje argumenty podane za pomocą następujących flag:
 ```sh
   -pliczek przyjmuje plik (koniecznie z rozszerzeniem), w którym znajdują się dane potrzebne do wykonania transformacji
   -elip przyjmuje nazwę modelu elipsoidy, na której chcemy dokonać transformacji
   -trans przyjmuje nazwę transformacji, którą chcemy wykonać
  ```
  
  Wybór elipsoidy możliwy jest poprzez wpisanie jednej z poniższych nazw:
  ```sh
   'WGS84'
   'GRS80'
   'Elipsoida Krasowskiego'
  ```
  
  Wybór transformacji możliwy jest poprzez wpisanie jednej z poniższych nazw:
  ```sh
   'XYZ_BLH'
   'BLH_XYZ'
   'XYZ_NEU'
   'BL_PL1992'
   'BL_PL2000'
  ```
  
  Po wyborze parametrów i załadowaniu pliku z danymi utworzy się plik tekstowy zawierający wyniki wykonanych obliczeń, a na konsoli pojawi się komunikat:
  ```sh
   Zapisano
  ```
  
  
   
 
