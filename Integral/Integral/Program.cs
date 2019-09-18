using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Integral
{
    class Program
    {
        static void LUm(double[,] A, ref double[,] L, ref double[,] U, ref int[,] P, int n)
        {
            double max = Math.Abs(A[0, 0]);
            int q = 0;

            for (int i = 0; i < n; i++)
            {
                P[i, i] = 1;
            }

            for (int j = 1; j < n; j++)
            {
                if (Math.Abs(A[j, 0]) > max)
                {
                    max = Math.Abs(A[j, 0]);
                    q = j;
                }
            }
            if (q != 0)
            {
                for (int i = 0; i < n; i++)
                {
                    double c;
                    int d;
                    d = P[0, i];
                    P[0, i] = P[q, i];
                    P[q, i] = d;
                    c = A[0, i];
                    A[0, i] = A[q, i];
                    A[q, i] = c;
                }
            }

            U[0, 0] = A[0, 0];
            L[0, 0] = 1;
            for (int j = 1; j < n; j++)
            {
                U[0, j] = A[0, j];
                L[j, 0] = A[j, 0] / U[0, 0];
            }
            for (int i = 1; i < n; i++)
            {
                max = Math.Abs(A[i, i]);
                q = i;
                for (int j = i + 1; j < n; j++)
                {
                    if (Math.Abs(A[j, i]) > max)
                    {
                        max = Math.Abs(A[j, i]);
                        q = j;
                    }
                }
                if (q != i)
                {
                    for (int j = 0; j < n; j++)
                    {
                        int d;
                        double c;
                        d = P[i, j];
                        P[i, j] = P[q, j];
                        P[q, j] = d;
                        c = A[i, j];
                        A[i, j] = A[q, j];
                        A[q, j] = c;
                        c = L[q, j];
                        L[q, j] = L[i, j];
                        L[i, j] = c;
                    }
                }

                U[i, i] = A[i, i];
                L[i, i] = 1;
                for (int k = 0; k <= i - 1; k++)
                    U[i, i] -= L[i, k] * U[k, i];
                for (int j = i + 1; j < n; j++)
                {
                    U[i, j] = A[i, j];
                    L[j, i] = A[j, i];
                    for (int k = 0; k <= i - 1; k++)
                    {
                        U[i, j] -= L[i, k] * U[k, j];
                        L[j, i] -= L[j, k] * U[k, i];
                    }
                    L[j, i] /= U[i, i];
                }
            }
        }


        static void LU(ref double[] x, double[,] A, double[] b, int n)
        {
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            double[] Pb = new double[n];
            double[] y = new double[n];
            int[,] P = new int[n, n];

            double[,] A1 = new double[n, n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    A1[i, j] = A[i, j];

            LUm(A1, ref L, ref U, ref P, n);


            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++)
                    Pb[i] += P[i, k] * b[k];


            for (int i = 0; i < n; i++)
            {
                double sum = Pb[i];
                for (int j = 0; j <= i - 1; j++)
                    sum -= L[i, j] * y[j];
                y[i] = sum;
            }

            for (int i = n - 1; i >= 0; i--)
            {
                double sum = y[i];
                for (int j = n - 1; j >= i; j--)
                    sum -= U[i, j] * x[j];
                sum /= U[i, i];
                x[i] = sum;
            }

        }

        static double[] Func(double[] x, int n)
        {
            double[] F = new double[n];
            for (int i = 0; i < n; i++)
            {
                F[i] = 2 * Math.Cos(2.5 * x[i]) * Math.Exp(x[i] / 3) + 4 * Math.Sin(3.5 * x[i]) * Math.Exp(-3 * x[i]) + x[i];
            }
            return F;
        }



        static double ICF(double[] x,double a1,double a,double b,double alfa,double beta,int n)
        {
            
            double[] mu = new double[n];
            mu[0] = 1 / (Math.Pow(b - a1, alfa - 1) * (1 - alfa)) - 1 / (Math.Pow(a - a1, alfa - 1) * (1 - alfa));
            for (int i = 1; i < n; i++)
            {
                double bb = Math.Pow(b - a1, 1 - alfa) * Math.Pow(b, Convert.ToDouble(i)) / (1 - alfa);
                double aa = Math.Pow(a - a1, 1 - alfa) * Math.Pow(a, Convert.ToDouble(i)) / (1 - alfa);
                mu[i] = (bb - aa + Convert.ToDouble(i) * a1 / (1 - alfa) * mu[i - 1]) / (1 + Convert.ToDouble(i) / (1 - alfa));
            }
            double[] A = new double[n];
            double[,] B = new double[n, n];
            for (int i=0;i<n;i++)
            {
                for (int j=0;j<n;j++)
                {
                    if (i == 0)
                        B[i, j] = 1;
                    else
                        B[i, j] = Math.Pow(x[j], i);
                }
            }
            LU(ref A, B, mu, n);
            double[] F = Func(x, n);
            double I = 0;
            for (int i=0;i<n;i++)
            {
                I += A[i] * F[i];
            }
            return I;
        }


        static double G(double[] x, double a1, double a, double b, double alfa, double beta, int n)
        {
            double[] mu = new double[2 * n];
            mu[0] = 1 / (Math.Pow(b - a1, alfa - 1) * (1 - alfa)) - 1 / (Math.Pow(a - a1, alfa - 1) * (1 - alfa));
            for (int i = 1; i <= 2 * n - 1; i++)
            {
                double bb = Math.Pow(b - a1, 1 - alfa) * Math.Pow(b, Convert.ToDouble(i)) / (1 - alfa);
                double aa = Math.Pow(a - a1, 1 - alfa) * Math.Pow(a, Convert.ToDouble(i)) / (1 - alfa);
                mu[i] = (bb - aa + Convert.ToDouble(i) * a1 / (1 - alfa) * mu[i - 1]) / (1 + Convert.ToDouble(i) / (1 - alfa));
            }
            double[,] A = new double[n, n];
            double[] B = new double[n];
            double[] ar = new double[n];
            for (int i = 0; i < n; i++)
            {
                B[i] = -mu[n + i];
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = mu[i + j];
                }
            }
            LU(ref ar, A, B, n);
            /*for (int i = n-1; i >=0; i--)
            {
                Console.WriteLine(ar[i]);
            }*/
            double[] x1 = new double[n];
            {
                double p = -ar[2] * ar[2] / 3 + ar[1];
                double q = 2 * Math.Pow(ar[2] / 3.0, 3.0) - ar[2] * ar[1] / 3.0 + ar[0];
                //double Q = Math.Pow(p / 3, 3) + Math.Pow(q/2, 2);
                double fi = Math.Acos(-(q / 2.0) * Math.Pow(3 / (-p), 1.5));

                x1[0] = 2 * Math.Pow(-p / 3.0, 0.5) * Math.Cos(fi / 3.0) - ar[2] / 3;
                x1[1] = 2 * Math.Pow(-p / 3, 0.5) * Math.Cos(fi / 3.0 + (2.0 / 3.0) * Math.PI) - ar[2] / 3;
                x1[2] = 2 * Math.Pow(-p / 3, 0.5) * Math.Cos(fi / 3.0 - (2.0 / 3.0) * Math.PI) - ar[2] / 3;
            }
            double[,] X1 = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                B[i] = mu[i];
                for (int j = 0; j < n; j++)
                {
                    if (i == 0)
                        X1[i, j] = 1;
                    else
                    {
                        double sum = 0;
                        sum = Math.Pow(x1[j], i);
                        X1[i, j] = Math.Pow(x1[j], i);
                    }
                    //Console.Write("{0} ", X1[i, j]);
                }
                //Console.WriteLine(B[i]);
            }
            double[] AA = new double[n];
            LU(ref AA, X1, B, n);
            /*for (int i=0;i<n;i++)
            {
                Console.WriteLine(AA[i]);
            }*/
            double[] F = Func(x, n);
            double I = 0;
            for (int i = 0; i < n; i++)
            {
                I += AA[i] * F[i];
            }
            return I;
        }




        static double NCI(double a, double b, double alfa, double beta, int n, int k)
        {
            double delta = (b - a) / k;
            double I = 0;
            for (int i = 1; i <= k; i++)
            {
                double[] x1 = new double[] { a + (i - 1) * delta, (2 * a + (2 * i - 1) * delta) / 2, a + i * delta };
                I += ICF(x1, a, a + (i - 1) * delta, a + i * delta, alfa, beta, n);
            }
            return I;
        }


        static double NCG(double a, double b, double alfa, double beta, int n, int k)
        {
            double delta = (b - a) / k;
            double I = 0;
            for (int i = 1; i <= k; i++)
            {
                double[] x1 = new double[] { a + (i - 1) * delta, (2 * a + (2 * i - 1) * delta) / 2, a + i * delta };
                I += G(x1, a, a + (i - 1) * delta, a + i * delta, alfa, beta, n);
            }
            return I;
        }



        

        static double NCIh(double a, double b, double alfa, double beta, int n, double h)
        {
            double I = 0;
            for (int i = 0; a+i*h <= b; i++)
            {
                double[] x1 = new double[n];
                if (a + (i + 1) * h < b)
                {
                    x1 = new double[] { a + i * h, (2 * a + (2 * i + 1) * h) / 2, a + (i + 1) * h };
                    I += ICF(x1, a, a + i * h, a + (i + 1) * h, alfa, beta, n);
                }
                else
                {
                    x1 = new double[] { a + i * h, ((a + i * h) + b) / 2, b };
                    I += ICF(x1, a, a + i * h, b, alfa, beta, n);
                }
            }
            return I;
        }

        static double NCGh(double a, double b, double alfa, double beta, int n, double h)
        {
            double I = 0;
            for (int i = 0; a + i * h <= b; i++)
            {
                double[] x1 = new double[n];
                if (a + (i + 1) * h < b)
                {
                    x1 = new double[] { a + i * h, (2 * a + (2 * i + 1) * h) / 2, a + (i + 1) * h };
                    I += G(x1, a, a + i * h, a + (i + 1) * h, alfa, beta, n);
                }
                else
                {
                    x1 = new double[] { a + i * h, ((a + i * h) + b) / 2, b };
                    I += G(x1, a, a + i * h, b, alfa, beta, n);
                }
            }
            return I;
        }

        static double ANCI(double a, double b, double alfa, double beta,double Eps, int n,bool opt)
        {
            double L = 2.0;
            double h1 = 1.0;
            double h2 = h1 / L;
            double h3 = h2 / L;
            double R = 2 * Eps;
            double S1=0, S2=0, S3=0;
            if (opt)
            {
                S1 = NCIh(a, b, alfa, beta, n, h1);
                S2 = NCIh(a, b, alfa, beta, n, h2);
                S3 = NCIh(a, b, alfa, beta, n, h3);
                double m = Atkin(S1, S2, S3, L);
                h3 = h1 * Math.Pow(Eps * (1 - Math.Pow(L, -m)) / Math.Abs(S2 - S1), 1.0 / m);
                h2 = h3 * L;
                h1 = h3 * (L * L);
                Console.WriteLine("hopt {0}", h3);
            }
            int count = 0;
            while (R>Eps)
            {
                S1 = NCIh(a, b, alfa, beta, n, h1);
                S2 = NCIh(a, b, alfa, beta, n, h2);
                S3 = NCIh(a, b, alfa, beta, n, h3);
                double m = Atkin(S1, S2, S3, L);
                R = Math.Abs(Richard(S2, S3, L, m));
                h1 = h2;
                h2 = h3;
                h3 = h3 / L;
                count++;
            }
            Console.WriteLine(count);
            Console.WriteLine(h3*L);
            Console.WriteLine(R);
            return S3;
        }

        /*static double ANCG(double a, double b, double alfa, double beta, double Eps, int n, bool opt)
        {
            double L = 2.0;
            double h1 = 1.0;
            double h2 = h1 / L;
            double h3 = h2 / L;
            double R = 2 * Eps;
            double S1 = 0, S2 = 0, S3 = 0;
            if (opt)
            {
                S1 = NCGh(a, b, alfa, beta, n, h1);
                S2 = NCGh(a, b, alfa, beta, n, h2);
                S3 = NCGh(a, b, alfa, beta, n, h3);
                double m = Atkin(S1, S2, S3, L);
                h3 = h1 * Math.Pow(Eps * (1 - Math.Pow(L, -m)) / Math.Abs(S2 - S1), 1.0 / m);
                h2 = h3 * L;
                h1 = h3 * (L * L);
                Console.WriteLine("hopt {0}", h3);
            }
            while (R > Eps)
            {
                S1 = NCGh(a, b, alfa, beta, n, h1);
                S2 = NCGh(a, b, alfa, beta, n, h2);
                S3 = NCGh(a, b, alfa, beta, n, h3);
                double m = Atkin(S1, S2, S3, L);

                R = Math.Abs(Richard(S2, S3, L, m));
                h1 = h2;
                h2 = h3;
                h3 = h3 / L;
            }
            Console.WriteLine(h3 * L);
            Console.WriteLine(R);
            return S3;
        }*/


        static double[] Kordano(double[] arr,int n )
        {
            double[] x = new double[n];
            double b = arr[2];
            double c = arr[1];
            double d =arr[0];
            double q = (Math.Pow(b, 2) - 3 * c) / 9.0;
            double r = (2 * Math.Pow(b, 3) - 9 * b * c + 27 * d) / 54.0;
            double s = Math.Pow(q, 3) - Math.Pow(r, 2);
            double phi = Math.Acos(r / (Math.Pow(q, (3.0 / 2.0)))) / 3.0;
            x[0] = -2 * (Math.Pow(q, (1.0 / 2.0)) * Math.Cos(phi)) - b / 3.0;
            x[1] = -2 * (Math.Pow(q, (1.0 / 2.0)) * Math.Cos(phi + 2 * Math.PI / 3.0)) - b / 3.0;
            x[2] = -2 * (Math.Pow(q, (1.0 / 2.0)) * Math.Cos(phi - 2 * Math.PI / 3.0)) - b / 3.0;
            return x;
        }

        static double Atkin(double S1,double S2,double S3,double L)
        {
            double m = -Math.Log((S3 - S2) / (S2 - S1)) / Math.Log(L);
            return m;
        }

        static double Richard(double S1,double S2,double L,double m)
        {
            double R = (S2 - S1) / (-1 + Math.Pow(L, m));
            return R;
        }

        static void Main(string[] args)
        {
            int n = 3;
            double a = 1.5;
            double a1 = a;
            double b = 3.3;
            double alfa = 1.0 / 3.0;
            double beta = 0;
            double J = 7.07703;
            double[] x = new double[] { a, (a + b) / 2, b };
            //Console.WriteLine("ИКФ {0}",ICF(x, a,a, b, alfa, beta,n));
            //Console.WriteLine("СКФ НК {0}",NCI(a, b, alfa, beta, n, 10));
            Console.WriteLine("СКФ НК нахождение с заданной точности");
            double S=ANCI(a, b, alfa, beta, 0.000001, n,false);
            Console.WriteLine(S);
            Console.WriteLine("СКФ НК hopt");
            S = ANCI(a, b, alfa, beta, 0.000001, n, true);
            Console.WriteLine(S);
            //Console.WriteLine(ANCG(a, b, alfa, beta, 0.000001, n, false));
            //Console.WriteLine(NCG(a, b, alfa, beta, n,1));
            Console.ReadKey();
        }
    }
}
