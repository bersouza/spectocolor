/////////////////////////////////////////////////
///  This is a program to convert a spectrum  ///
///  file (intensity by wavelength) into its  ///
///  color following the CIE standards.       ///
///                                           ///
///  Programmed by Bernardo de Souza, 2017,   ///
///  Brazil.                                  ///
/////////////////////////////////////////////////

#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{
    //Introductory phrase
    cout << "\n\t\t**************************\t\t" << endl;
    cout << "\t\t* Welcome to spectocolor *\t\t" << endl;
    cout << "\t\t**************************\t\t\n" << endl;

    //Check if a filename was given
    if (argc < 2)
    {
        cout << "We need a spectrum file as input!" << endl;
        return 1;
    }

    //Open spectrum file
    ifstream f (argv[1]);

cout << "Could not open the file " << argv[1] << "!" << endl;

    //If file is not open, close
    if (f.is_open())
    {
        cout << "Could not open the file " << argv[1] << "!" << endl;
        return 1;
    }

    return 0;
}
