/////////////////////////////////////////////////
///  This is a program to convert a spectrum  ///
///  file (intensity by wavelength) into its  ///
///  color, following the CIE color matching  ///
///  functions from 1931, using their         ///
///  analytic form as given by Wiman, C et    ///
///  all on Journal of Computer Graphics      ///
///  Tecniches, v.2, n.2, 2013.               ///
///                                           ///
///  Programmed by Bernardo de Souza, 2017,   ///
///  Brazil.                                  ///
/////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#define eV_cm 8065.540106923572

using namespace std;

//
// Function to count the number of lines (c) in file f
//
int CountLines(ifstream &f,auto &c)
{
    string text;

    while (!f.eof())
    {
        f >> text >> text;
        c++;
    }

    //After counting
    f.clear();
    f.seekg(0);
    getline(f,text);

    return 0;
}

//
// This function converts the X axis to nm
//
int ConvertXaxis(char *argv, double *xaxis, auto &c)
{
    char inp[512];

    //Convert cm-1 to nm
    sprintf(inp,"-cm-1");
    if (argv == inp)
    for (auto i=0;i<c;++i)
        xaxis[i] = 1e7/xaxis[i];

    //Convert eV to nm
    sprintf(inp,"-eV");
    if (argv == inp)
        for (auto i=0;i<c;++i)
            xaxis[i] = 1e7/(xaxis[i]*eV_cm);

    return 0;
}

//
// This function normalizes the Y axis
int NormalizeYaxis(double *yaxis, auto &c)
{
    double max=0;
    //Find the largest element
    for (auto i=0;i<c;++i)
        max = (yaxis[i]>max)?yaxis[i]:max;

    //Nomalize
    for (auto i=0;i<c;++i)
        yaxis[i] /= max;

    return 0;
}

//
// Gets the color matching function value for a given wavenumber in nm
//
double CMFX(double x)
{
    double t1 = (x-442.0)*((x<442.0)?0.0624:0.0374);
    double t2 = (x-599.8)*((x<599.8)?0.0264:0.0323);
    double t3 = (x-501.1)*((x<501.1)?0.0490:0.0382);
    return 0.362*exp(-0.5*t1*t1) + 1.056*exp(-0.5*t2*t2)
            - 0.065*exp(-0.5*t3*t3);
}

//
// Gets the color matching function value for a given wavenumber in nm
//
double CMFY(double x)
{
    double t1 = (x-568.8)*((x<568.8)?0.0213:0.0247);
    double t2 = (x-530.9)*((x<530.9)?0.0613:0.0322);
    return 0.821*exp(-0.5*t1*t1) + 0.286*exp(-0.5*t2*t2);
}

//
// Gets the color matching function value for a given wavenumber in nm
//
double CMFZ(double x)
{
    double t1 = (x-437.0)*((x<437.0)?0.0845:0.0278);
    double t2 = (x-459.0)*((x<459.0)?0.0385:0.0725);
    return 1.217*exp(-0.5*t1*t1) + 0.681*exp(-0.5*t2*t2);
}

//
// Generates the color matching functions
// that correspond to xaxis
//
int GenerateCMFs(double *xaxis,double *cmf_x,double *cmf_y,double *cmf_z, auto &c)
{
    for (auto i=0;i<c;++i)
    {
        cmf_x[i] = CMFX(xaxis[i]);
        cmf_y[i] = CMFY(xaxis[i]);
        cmf_z[i] = CMFZ(xaxis[i]);
    }

    return 0;
}

//
// This calculates the products yaxis*(color matching functions)
//
int CalcProducts(double *yaxis,double *cmf_x, double *cmf_y, double *cmf_z,
                 double *pX,double *pY,double *pZ,auto &c)
{
    for (auto i=0;i<c;++i)
    {
        pX[i] = yaxis[i]*cmf_x[i];
        pY[i] = yaxis[i]*cmf_y[i];
        pZ[i] = yaxis[i]*cmf_z[i];
    }

    return 0;
}

//
// This integrates p over the X axis
//
double Integrate(double *xaxis,double *p, auto &c)
{
    auto sum=0.0;
    for (auto i=0;i<(c-1);++i)
        sum += ((p[i]+p[i+1])/2.0)*fabs((xaxis[i+1]-xaxis[i]));

    return sum;
}

int ConverttoRBG(double X, double Y, double Z,
                 double &R, double &G, double &B)
{
    X /= Y;
    Z /= Y;
    Y /= Y;

    R = X*3.2404542 + Y*-1.5371385 + Z*-0.4985312;
    G = X*-0.9692660 + Y*1.8760108 + Z*0.0415560;
    B = X*0.0556434 + Y*-0.2040259 + Z*1.0572252;

    if (R > 0.0031308)
        R = 1.055*(pow(R,(1.0/2.4)) - 0.055);
    else
        R = 12.92*R;

    if (G > 0.0031308)
        G = 1.055*(pow(G,(1.0/2.4)) - 0.055);
    else
        G = 12.92*G;

    if (B > 0.0031308)
        B = 1.055*(pow(B,(1.0/2.4)) - 0.055);
    else
        B = 12.92*B;

    R = R*255;
    G = G*255;
    B = B*255;

    return 0;
}

int ConverttoHSL(double R,double G,double B,
                 double &H, double &S, double &L)
{
    R /= 255;
    G /= 255;
    B /= 255;

    auto Min = min(min(R,G),B);                    //Min. value of RGB
    auto Max = max(max(R,G),B);                    //Max. value of RGB

    auto del_Max = Max - Min;    //Delta RGB value

    L = (Max + Min)/ 2;

    if (del_Max == 0)                     //This is a gray, no chroma...
    {
        H = 0;
        S = 0;
    }
    else                                    //Chromatic data...
    {
       if (L < 0.5)
           S = del_Max/(Max + Min);
       else
           S = del_Max/(2 - Max - Min);

       auto del_R = (((Max - R)/6.0) + (del_Max/2.0))/del_Max;
       auto del_G = (((Max - G)/6.0) + (del_Max/2.0))/del_Max;
       auto del_B = (((Max - B)/6.0) + (del_Max/2.0))/del_Max;

       if (R == Max)
           H = del_B - del_G;
       else

           if (G == Max )
               H = (1.0/3.0) + del_R - del_B;
           else
               if (B == Max )
                   H = (2.0/3.0) + del_G - del_R;

        if (H < 0)
            H += 1;
        if (H > 1)
            H -= 1;
    }

    return 0;
}

int ConverttoCMYK(double R,double G,double B,
                  double &C,double &M,double &Y,double &K)
{
    C = 1 - (R/255);
    M = 1 - (G/255);
    Y = 1 - (B/255);

    K = 1;

    if (C < K)
        K = C;
    if (M < K)
        K = M;
    if (Y < K)
        K = Y;

    if (K == 1)
    {
        C = 0;          //Black only
        M = 0;
        Y = 0;
    }
    else
    {
        C = (C - K)/(1 - K);
        M = (M - K)/(1 - K);
        Y = (Y - K)/(1 - K);
    }

    C *= 100;
    M *= 100;
    Y *= 100;
    K *= 100;

    return 0;
}

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

    //If file is not open, close
    if (!f.is_open())
    {
        cout << "Could not open the file " << argv[1] << "!" << endl;
        return 1;
    }

    //If everyting is OK, skip the first line and read all the rest
    string text;
    getline(f,text);

    //Count the number of lines
    auto c=-1;
    CountLines(f,c);

    //Read the input spectrum
    double xaxis[c],yaxis[c];
    for (auto i=0;i<c;++i)
        f >> xaxis[i] >> yaxis[i];

    //If cm-1 or eV, convert the X axis to nm
    ConvertXaxis(argv[2],xaxis,c);

    //Normalize the Y axis
    NormalizeYaxis(yaxis,c);

    //Create the corresponding color matching functions
    double cmf_x[c],cmf_y[c],cmf_z[c];
    GenerateCMFs(xaxis,cmf_x,cmf_y,cmf_z,c);

    //Compute the products pX,pY and pZ
    double pX[c],pY[c],pZ[c];
    CalcProducts(yaxis,cmf_x,cmf_y,cmf_z,pX,pY,pZ,c);

    //Compute integrals X, Y and Z
    double X,Y,Z;
    X = Integrate(xaxis,pX,c);
    Y = Integrate(xaxis,pY,c);
    Z = Integrate(xaxis,pZ,c);

    //Convert to Yxy
    double x,y;
    x = X/(X+Y+Z);
    y = Y/(X+Y+Z);

    //Convert to RBG
    double R,G,B;
    ConverttoRBG(X,Y,Z,R,G,B);

    //Convert to HSL
    double H,S,L;
    ConverttoHSL(R,G,B,H,S,L);

    //Convert to CMYK
    double C,M,K;
    ConverttoCMYK(R,G,B,C,M,Y,K);

    /*for (auto i=0;i<c;++i)
        cout << xaxis[i] << "\t" << yaxis[i] << "\t" << pX[i] << "\t" << pY[i] << "\t" << pZ[i]
                << "\t" << endl;

    cout << X << " " << Y << " " << Z << endl;*/

    //Print results
    cout << "Calculating the color from the spectra...\n" << endl;

    cout << "XYZ color space:\n" << endl;
    cout << "\tX: " << X << endl;
    cout << "\tY: " << Y << endl;
    cout << "\tZ: " << Z << endl;

    cout << "\nYxy color space:\n" << endl;
    cout << "\tY: " << Y << endl;
    cout << "\tx: " << x << endl;
    cout << "\ty: " << y << endl;

    cout << "\nsRGB color space [0-255]:\n" << endl;
    cout << "\tR: " << R << endl;
    cout << "\tG: " << G << endl;
    cout << "\tB: " << B << endl;

    cout << "\nHSL color space [0-1]:\n" << endl;
    cout << "\tH: " << H << endl;
    cout << "\tS: " << S << endl;
    cout << "\tL: " << L << endl;

    cout << "\nCMYK color space [0-1]:\n" << endl;
    cout << "\tC: " << C << endl;
    cout << "\tM: " << M << endl;
    cout << "\tY: " << Y << endl;
    cout << "\tK: " << K << endl;

    return 0;
}
