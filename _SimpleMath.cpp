//---------------------------------------------------------------------------

#define BUILD_DLL
#include <windows.h>
#include <math.h>
#pragma hdrstp
#include "SimpleMath.h"

//---------------------------------------------------------------------------

double Trap(double a, double b, double c)
{
  double s=(a+b)*c/2;
  return s;
}

//---------------------------------------------------------------------------

double Lin(double x1, double x2, double y1, double y2, double x)
{
  double y=(y2-y1)*(x-x1)/(x2-x1)+y1;
  return y;
}

//---------------------------------------------------------------------------

double det3(double mtx[3][3])
{
  double d=0;

  for (int i=0; i<3; i++)
  {d+=mtx[i][0]*(mtx[(i+1)%3][1]*mtx[(i+2)%3][2]-mtx[(i+2)%3][1]*mtx[(i+1)%3][2]);}

  return d;
}

//---------------------------------------------------------------------------

void SqrApp(int N, double** (&Func), double& a, double& b, double& c)
{
   // String to be changed in order to test git

   double xs[5];
   double ys[3];
   double mtx[3][3];

   double da=0;
   double db=0;
   double dc=0;
   double dd=0;

   for (int i=0; i<5; i++)
     {xs[i]=0;}
   for (int i=0; i<3; i++)
     {ys[i]=0;}

   for (int i=0; i<N; i++)
   {
       xs[0]++;
       xs[1]+=Func[0][i];
       xs[2]+=Func[0][i]*Func[0][i];
       xs[3]+=Func[0][i]*Func[0][i]*Func[0][i];
       xs[4]+=Func[0][i]*Func[0][i]*Func[0][i]*Func[0][i];

       ys[0]+=Func[1][i];
       ys[1]+=Func[1][i]*Func[0][i];
       ys[2]+=Func[1][i]*Func[0][i]*Func[0][i];
   }

   mtx[0][0]=xs[4]; mtx[1][0]=xs[3]; mtx[2][0]=xs[2];
   mtx[0][1]=xs[3]; mtx[1][1]=xs[2]; mtx[2][1]=xs[1];
   mtx[0][2]=xs[2]; mtx[1][2]=xs[1]; mtx[2][2]=xs[0];
   dd=det3(mtx);
   mtx[0][0]=ys[2]; mtx[1][0]=xs[3]; mtx[2][0]=xs[2];
   mtx[0][1]=ys[1]; mtx[1][1]=xs[2]; mtx[2][1]=xs[1];
   mtx[0][2]=ys[0]; mtx[1][2]=xs[1]; mtx[2][2]=xs[0];
   da=det3(mtx);
   mtx[0][0]=xs[4]; mtx[1][0]=ys[2]; mtx[2][0]=xs[2];
   mtx[0][1]=xs[3]; mtx[1][1]=ys[1]; mtx[2][1]=xs[1];
   mtx[0][2]=xs[2]; mtx[1][2]=ys[0]; mtx[2][2]=xs[0];
   db=det3(mtx);
   mtx[0][0]=xs[4]; mtx[1][0]=xs[3]; mtx[2][0]=ys[2];
   mtx[0][1]=xs[3]; mtx[1][1]=xs[2]; mtx[2][1]=ys[1];
   mtx[0][2]=xs[2]; mtx[1][2]=xs[1]; mtx[2][2]=ys[0];
   dc=det3(mtx);

   if (dd!=0)
   {
     a=da/dd;
     b=db/dd;
     c=dc/dd;
   }
}

//---------------------------------------------------------------------------

void LineApp(int N, double** (&Line), double& a, double& b)
{
   // Linear approximation

   double sumX=0;
   double sumY=0;
   double sumX2=0;
   double sumXY=0;

   for (int i=0; i<N; i++)
   {
       sumX+=Line[0][i];
       sumY+=Line[1][i];
       sumX2+=Line[0][i]*Line[0][i];
       sumXY+=Line[0][i]*Line[1][i];
   }

   double dd=sumX2*N-sumX*sumX;
   double da=sumXY*N-sumX*sumY;
   double db=sumX*sumXY-sumY*sumX2;

   if (dd!=0)
   {
     a = da/dd;
     b = da/dd;
   }
}

//---------------------------------------------------------------------------

void QuickSort(double** (&items), int left, int right)
{
  register int i, j;
  double x, y;

  i = left;
  j = right;

  double l = items[0][(left+right)/2];

  do {
    while((items[0][i] < l) && (i < right)) i++;
    while((l < items[0][j]) && (j > left)) j--;

    if(i <= j) {
      x = items[0][i];
      y = items[1][i];
      items[0][i] = items[0][j];
      items[1][i] = items[1][j];
      items[0][j] = x;
      items[1][j] = y;
      i++; j--;
    }
  } while(i <= j);

  if(left < j) QuickSort(items, left, j);
  if(i < right) QuickSort(items, i, right);
}

//---------------------------------------------------------------------------

void FourFrec(int Left, int Right, double** (&inputLSF), double f, double delta, double& re,double& im)
{
  double ret=0;
  double imt=0;

  for (int j=Left; j<Right; j++)
  {
     ret+=inputLSF[1][j]*cos(2.0*M_PI*f*inputLSF[0][j])*delta;
     imt+=-inputLSF[1][j]*sin(2.0*M_PI*f*inputLSF[0][j])*delta;
  }

  re=ret;
  im=imt;
}

//------------------------------------------------------------------------------

#pragma argsused
BOOL WINAPI DllMain(HINSTANCE hinstDLL, DWORD fwdreason, LPVOID lpvReserved)
{
        return 1;
}

//---------------------------------------------------------------------------
