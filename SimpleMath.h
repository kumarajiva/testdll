//---------------------------------------------------------------------------

#if defined(BUILD_DLL)
# define DLL_EXP __declspec(dllexport)
#else
# if defined(BUILD_APP)
# define DLL_EXP __declspec(dllimport)
# else
# define DLL_EXP
# endif
#endif

//---------------------------------------------------------------------------


//--------------Approximating Functions--------------------------------------
void DLL_EXP LineApp(int N, double** (&Line), double& a, double& b);
void DLL_EXP SqrApp(int N, double** (&Func), double& a, double& b, double& c);
//---------------------------------------------------------------------------

//--------------Sorting Functions--------------------------------------------
void DLL_EXP QuickSort(double** (&items), int left, int right);
//---------------------------------------------------------------------------

//--------------Point Approximating Functions--------------------------------
double DLL_EXP Trap(double a, double b, double c);
double DLL_EXP Lin(double x1, double x2, double y1, double y2, double x);
//---------------------------------------------------------------------------

//--------------Matrix Operating Functions-----------------------------------
double DLL_EXP det3(double mtx[3][3]);
//---------------------------------------------------------------------------

//--------------Fourier Transform Functions----------------------------------
void DLL_EXP FourFrec(int Left, int Right, double** (&inputLSF), double f, double delta, double& re,double& im);
//---------------------------------------------------------------------------

