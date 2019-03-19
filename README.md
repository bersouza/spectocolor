# spectocolor
a program to convert absorbance or emission spectrum to its color following CIE stadards

To use, just run ./spectocolor filename.txt

The .txt file should be composed of two columns, x and y. The first line will be not read. 

As for the x axis, the default scale is nanometers. If one has a spectrum in cm-1 or eV, one should use the options:

./spectocolor filename.txt -cm-1 or ./spectocolor filename.txt -eV

Further references are on the .cpp file.
