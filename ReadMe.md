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
   system operacyjny Windows 11 lub macOS 11.7.6
   python 3.9 lub python 3.8.9
   biblioteka numpy
   biblioteka argparse
 ``` 
  
 Opis programu:
 
 Plik przyjmuje argumenty podane za pomocą następujących flag:
 ```sh
   -plik przyjmuje plik (koniecznie z rozszerzeniem), w którym znajdują się dane potrzebne do wykonania transformacji
   -elip przyjmuje nazwę modelu elipsoidy, na której chcemy dokonać transformacji
   -funkcja przyjmuje nazwę transformacji, którą chcemy wykonać
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
  Pliku wynikowy zapisuje się pod nazwą:
  ```sh
  "WYNIK_{funkcja}.txt"
  (gdzie {funkcja} oznacza nazwę transformacji, którą chcemy wykonać)
  ```
  
  
  Przykładowe transformacje:
  
  XYZ ---> BLH
  dla danych z pliku (kolejno X[m], Y[m], Z[m])
  ```sh
  12345.789 65789.003 36674.123
  11111.222 22222.343 55555.332
  ```
  otrzymujemy wyniki (kolejno fi, lambda, h)
  ```sh
  4.498849871596222982e+01   7.937163647711760461e+01   -6.294189244106777012e+06
  7.564029626858634003e+01   6.343484465985566345e+01   -6.298088199815715663e+06
  ```
  
  BLH ---> XYZ
  dla danych z pliku (kolejno fi, lambda, h)
  ```sh
  52 21 319
  52 19 420
  ```
  otrzymujemy wyniki (kolejno X[m], Y[m], Z[m])
  ```sh
  3.673785422237344086e+06   1.410234096034315648e+06   5.003054720799141563e+06
  3.720822705660262145e+06   1.281182001713992795e+06   5.003134309885255992e+06
  ```
  
  XYZ,X0Y0Z0 ---> neu
  dla danych z pliku 'wsp_XYZ_NEU.txt' (kolejno fi, lambda h,)
  ```sh
   7.221181885654304642e+05
  -2.477294595237924659e+05 
  -1.045556128732857667e+07
  ```
  otrzymujemy wyniki (kolejno n, e, u)
  
  Ważne jest, aby współrzędne punktów podane zostały w odpowiedniej kolejności - jako pierwsze podać należy współrzędne początku układu NEU (x0, y0), a dopiero potem współrzędne, do których policzyć chcemy wektor. Stąd, żeby otrzymać jeden punkt wyjściowy należy wprowadzić dane aż dwóch punktów wejściowych.
  
  BL ---> XY PL1992
  dla danych z pliku 'bl-pl.txt' (kolejno fi, lambda)
  ```sh
  52  21
  52  17
  ```
  otrzymujemy wyniki (kolejno x92[m], y92[m])
  ```sh
  4.611972433942528442e+05   6.372531611049008789e+05
  4.611972433942528442e+05   3.627468388950995868e+05
  ```
  
  BL ---> XY PL2000
  dla danych z pliku 'bl-pl.txt' (kolejno fi, lambda)
  ```sh
  52  21
  52  17
  ```
  otrzymujemy wyniki (kolejno x00[m], y00[m])
  ```sh
  5.762899772909278050e+06   7.500000000000000000e+06
  5.763372029424777254e+06   6.431328112376346253e+06
  ```
   
   
   ZNANE BŁĘDY
   ```sh
      - zawarte w kodzie transformacje nie działają, gdy użytkownik próbuje wykonać je na elipsoidzie Krasowskiego
      - tranformacje PL1992 i PL2000 dają niepoprawne wyniki, gdy użytkownik próbuje wykonać je na elipspidzie Krasowskiego
   ```
 
