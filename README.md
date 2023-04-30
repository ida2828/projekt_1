INSTRUKCJA

1.Zasotsowanie programu:
-przeliczenie współrzędnych geocentrycznych XYZ na współrzędne geodezyjne BLH

-przeliczenie współrzędnych geodezyjnych BLH do współrzędnych geocentrycznych XYZ

-przeliczenie współrzędnych geocentrycznych XYZ do współrzędnych topocentrcznych NEU

-przeliczenie współrzędnych geodezyjnych BLHdo współrzędnych prostokątnych w układzie 2000

-przeliczenie współrzędnych geodezyjnych BLH do współrzędnych prostokątnych w układzie 1992



Program obsługuje parametry modeli elipsoid: GRS80, WGS84, Krassowskiego

2.Warunki działania programu:

-system operacyjny Windows

-zainstalowany na danym urządzeniu Python w wersji 3.9 oraz Spyder

-biblioteka Numpy oraz Argparse

3.Praca z programem:


Plik wejścowy może mieć rozszerzenie .txt i być zapisnay w tym folderze co program. 
Polecenia wykonujemy w cmd- wierszu poleceń. Program korzysta z biblioteki argparse - zatem przy wywołaniu podajemy argumenty. Przy wywołaniu argumentu "-m" wybieramy model elipsoid (wgs84, grs80, Krasowskiego) na której dokonujemy obliczeń. Argument "-t" oznacza rodzaj transformacji z której chcemy skorzystać (XYZ2BLH , BLH2XYZ, XYZ2NEU, u2000, u1992) oraz "-f" ścieżkę do pliku z danymi do transformacji
Dane zapisujemy w określony sposób:

a)XYZ2BLH
Współrzędne geocentryczne koljeno: X , Y , Z oddzielone od siebie przecinakmi, a części dziesiętne oddzielone kropką przykładowo: 3782580.000,1084640.000,5002880.000

b)BLH2XYZ
Współrzędne geodezyjne kolejno fi, lambda, ha (w radianach) oddzielone od siebie przecinakmi, a części dziesiętne oddzielone kropką przykładowo:0.7100,1.2915,46.15

c)XYZ2neu
Współrzędne geocentryczne XYZ początka odcinka i współrzędnych geocentrycznych XYZ oraz  współrzędne geodezyjne fi, lambda, ha. Zapisane w następującej kolejności: Xpoczątkowe, Ypoczątkowe, Zpoczątkowe, Xkońcowe, Ykońcowe, Zkońcowe, fi, lambda, ha. Przy czym współrzędne musza być one oddzielone od siebie przecinkami, a części dziesiętne oddzielone od jedności kropką. Przykładowo:3782520.000,1084820.000,5003000.000,3782647.600,1084331.060,5003550.000,66,19,150

d)u2000
Współrzędne geodezyjne kolejno fi, lambda, ha (w radianach) oddzielone od siebie przecinakmi, a części dziesiętne oddzielone kropką przykładowo: 0.3320,2.0001,12.9900

e)u1992
Współrzędne geodezyjne kolejno fi, lambda (w radianach) oddzielone od siebie przecinakmi, a części dziesiętne oddzielone kropką przykładowo: 0.9200,0.2643
