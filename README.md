INSTRUKCJA

1.Zasotsowanie programu:
-przeliczenie współrzędnych geocentrycznych XYZ na współrzędne geodezyjne BLH
-przeliczenie współrzędnych geodezyjnych BLHdo współrzędnych prostokątnych w układzie 2000
-przeliczenie współrzędnych geodezyjnych BLH do współrzędnych prostokątnych w układzie 1992
-przeliczenie współrzędnych geodezyjnych BLH do współrzędnych geocentrycznych XYZ
-przeliczenie współrzędnych geocentrycznych XYZ do współrzędnych topocentrcznych NEU
-przeliczenie współrzędnych geodezyjnych BLH do współrzędnych geocentrycznych XYZ

Program obsługuje parametry modeli elipsoid: GRS80, WGS84, Krassowskiego

2.Warunki działania programu:
-system operacyjny Windows
-zainstalowany na danym urządzeniu Python w wersji 3.9 oraz Spyder
-biblioteka Numpy

3.Praca z programem:


Plik wejścowy muzi mieć rozszerzenie .txt i być zapisnay w tym folderze co program. 
Polecenia wykonujemy w cmd- wierszu poleceń. Program korzysta z biblioteki argoarse - zatem przy wywołaniu podajemy argumenty. Przy wywołaniu argumentu "-m" wybieramy model elipsoid (wgs84, grs80, Krasowskiego) na której dokonujemy obliczeń. Argument "-t" oznacza rodzaj transformacji z której chcemy skorzystać (XYZ2BLH , BLH2XYZ, XYZ2NEU, u2000, u1992) oraz "-f" ścieżkę do pliku z danymi do transformacji
