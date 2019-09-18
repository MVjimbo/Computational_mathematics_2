using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace MetN
{
    class Program
    {

        static void LUm(double[,] A, ref double[,] L, ref double[,] U, ref int[,] P, int n)
        {
            double max = Math.Abs(A[0, 0]);
            int q = 0;

            for (int i=0;i<n;i++)
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


        static void LU (ref double[] x,double[,] A,double[] b,int n)
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

        static void LU2(double[,] L,double[,] U, int[,] P,ref double[] x, double[,] A, double[] b, int n)
        {
            double[] Pb = new double[n];
            double[] y = new double[n];

            double[,] A1 = new double[n, n];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    A1[i, j] = A[i, j];


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


        static double[,] Yakobi(double[] x,int n)
        {
            double[,] cg = new double[n, n];
            cg[0, 0] = -Math.Sin(x[0] * x[1]) * x[1];
            cg[0, 1] = -Math.Sin(x[0] * x[1]) * x[0];
            cg[0, 2] = 3.0 * Math.Exp(-3.0 * x[2]);
            cg[0, 3] = x[4] * x[4];
            cg[0, 4] = 2.0 * x[3] * x[4];
            cg[0, 5] = -1.0;
            cg[0, 6] = 0.0;
            cg[0, 7] = -2.0 * Math.Cosh((2.0 * x[7])) * x[8];
            cg[0, 8] = -Math.Sinh((2.0 * x[7]));
            cg[0, 9] = 2.0;
            cg[1, 0] = Math.Cos(x[0] * x[1]) * x[1];
            cg[1, 1] = Math.Cos(x[0] * x[1]) * x[0];
            cg[1, 2] = x[8] * x[6];
            cg[1, 3] = 0.0;
            cg[1, 4] = 6.0 * x[4];
            cg[1, 5] = (-1.0) * Math.Exp(-x[9] + x[5]) - x[7] - 1;
            cg[1, 6] = x[2] * x[8];
            cg[1, 7] = (-1.0) * x[5];
            cg[1, 8] = x[2] * x[6];
            cg[1, 9] = Math.Exp(-x[9] + x[5]);
            cg[2, 0] = 1.0;
            cg[2, 1] = -1.0;
            cg[2, 2] = 1.0;
            cg[2, 3] = -1.0;
            cg[2, 4] = 1.0;
            cg[2, 5] = -1.0;
            cg[2, 6] = 1.0;
            cg[2, 7] = -1.0;
            cg[2, 8] = 1.0;
            cg[2, 9] = -1.0;
            cg[3, 0] = -x[4] * 1.0 / ((x[2] + x[0]) * (x[2] + x[0]));
            cg[3, 1] = -2.0 * Math.Cos(x[1] * x[1]) * x[1];
            cg[3, 2] = -x[4] * 1.0 / ((x[2] + x[0]) * (x[2] + x[0]));
            cg[3, 3] = -2.0 * Math.Sin(-x[8] + x[3]);
            cg[3, 4] = 1.0 / (x[2] + x[0]);
            cg[3, 5] = 0.0;
            cg[3, 6] = -2.0 * Math.Cos(x[6] * x[9]) * Math.Sin(x[6] * x[9]) * x[9];
            cg[3, 7] = -1.0;
            cg[3, 8] = 2.0 * Math.Sin(-x[8] + x[3]);
            cg[3, 9] = -2.0 * Math.Cos(x[6] * x[9]) * Math.Sin(x[6] * x[9]) * x[6];
            cg[4, 0] = 2.0 * x[7];
            cg[4, 1] = -2.0 * Math.Sin(x[1]);
            cg[4, 2] = 2.0 * x[7];
            cg[4, 3] = 1.0 / ((-x[8] + x[3]) * (-x[8] + x[3]));
            cg[4, 4] = Math.Cos(x[4]);
            cg[4, 5] = x[6] * Math.Exp(-x[6] * (-x[9] + x[5]));
            cg[4, 6] = -(x[9] - x[5]) * Math.Exp(-x[6] * (-x[9] + x[5]));
            cg[4, 7] = (2.0 * x[2]) + 2.0 * x[0];
            cg[4, 8] = -1.0 / ((-x[8] + x[3]) * (-x[8] + x[3]));
            cg[4, 9] = -x[6] * Math.Exp(-x[6] * (-x[9] + x[5]));
            cg[5, 0] = Math.Exp(x[0] - x[3] - x[8]);
            cg[5, 1] = -1.5 * Math.Sin(3.0 * x[9] * x[1]) * x[9];
            cg[5, 2] = -x[5];
            cg[5, 3] = -Math.Exp(x[0] - x[3] - x[8]);
            cg[5, 4] = 2.0 * x[4] / x[7];
            cg[5, 5] = -x[2];
            cg[5, 6] = 0.0;
            cg[5, 7] = -x[4] * x[4] / x[7] / x[7];
            cg[5, 8] = -Math.Exp(x[0] - x[3] - x[8]);
            cg[5, 9] = -1.5 * Math.Sin(3.0 * x[9] * x[1]) * x[1];
            cg[6, 0] = Math.Cos(x[3]);
            cg[6, 1] = 3.0 * x[1] * x[1] * x[6];
            cg[6, 2] = 1.0;
            cg[6, 3] = -(x[0] - x[5]) * Math.Sin(x[3]);
            cg[6, 4] = Math.Cos(x[9] / x[4] + x[7]) * x[9] / x[4] / x[4];
            cg[6, 5] = -Math.Cos(x[3]);
            cg[6, 6] = Math.Pow(x[1], 3.0);
            cg[6, 7] = -Math.Cos(x[9] / x[4] + x[7]);
            cg[6, 8] = 0.0;
            cg[6, 9] = -Math.Cos(x[9] / x[4] + x[7]) / x[4];
            cg[7, 0] = 2.0 * x[4] * (x[0] - 2.0 * x[5]);
            cg[7, 1] = -x[6] * Math.Exp(x[1] * x[6] + x[9]);
            cg[7, 2] = -2.0 * Math.Cos(-x[8] + x[2]);
            cg[7, 3] = 1.5;
            cg[7, 4] = Math.Pow(x[0] - 2.0 * x[5], 2.0);
            cg[7, 5] = -4.0 * x[4] * (x[0] - 2.0 * x[5]);
            cg[7, 6] = -x[1] * Math.Exp(x[1] * x[6] + x[9]);
            cg[7, 7] = 0.0;
            cg[7, 8] = 2.0 * Math.Cos(-x[8] + x[2]);
            cg[7, 9] = -Math.Exp(x[1] * x[6] + x[9]);
            cg[8, 0] = -3.0;
            cg[8, 1] = -2.0 * x[7] * x[9] * x[6];
            cg[8, 2] = 0.0;
            cg[8, 3] = Math.Exp((x[4] + x[3]));
            cg[8, 4] = Math.Exp((x[4] + x[3]));
            cg[8, 5] = -7.0 / x[5] / x[5];
            cg[8, 6] = -2.0 * x[1] * x[7] * x[9];
            cg[8, 7] = -2.0 * x[1] * x[9] * x[6];
            cg[8, 8] = 3.0;
            cg[8, 9] = -2.0 * x[1] * x[7] * x[6];
            cg[9, 0] = x[9];
            cg[9, 1] = x[8];
            cg[9, 2] = -x[7];
            cg[9, 3] = Math.Cos(x[3] + x[4] + x[5]) * x[6];
            cg[9, 4] = Math.Cos(x[3] + x[4] + x[5]) * x[6];
            cg[9, 5] = Math.Cos(x[3] + x[4] + x[5]) * x[6];
            cg[9, 6] = Math.Sin(x[3] + x[4] + x[5]);
            cg[9, 7] = -x[2];
            cg[9, 8] = x[1];
            cg[9, 9] = x[0];
            return cg;
        }

        static double[] Func(double[] x,int n)
        {
            double[] F = new double[n];
            F[0] = -(Math.Cos(x[0] * x[1]) - Math.Exp(-3 * x[2]) + x[3] * Math.Pow(x[4], 2) - x[5] - Math.Sinh(2 * x[7]) * x[8] + 2 * x[9] + 2.0004339741653854440);
            F[1] = -(Math.Sin(x[0] * x[1]) + x[2] * x[8] * x[6] - Math.Exp(-x[9] + x[5]) + 3 * Math.Pow(x[4], 2) - x[5] * (x[7] + 1) + 10.886272036407019994);
            F[2] = -(x[0] - x[1] + x[2] - x[3] + x[4] - x[5] + x[6] - x[7] + x[8] - x[9] - 3.1361904761904761904);
            F[3] = -(2 * Math.Cos(-x[8] + x[3]) + x[4] / (x[2] + x[0]) - Math.Sin(x[1] * x[1]) + Math.Pow(Math.Cos(x[6] * x[9]), 2) - x[7] - 0.1707472705022304757);
            F[4] = -(Math.Sin(x[4]) + 2 * x[7] * (x[2] + x[0]) - Math.Exp(-x[6] * (-x[9] + x[5])) + 2 * Math.Cos(x[1]) - 1 / (x[3] - x[8]) - 0.368589627310127786);
            F[5] = -(Math.Exp(x[0] - x[3] - x[8]) + x[4] * x[4] / x[7] + 0.5 * Math.Cos(3 * x[9] * x[1]) - x[5] * x[2] + 2.049108601677187511);
            F[6] = -x[1] * x[1] * x[1] * x[6];
            F[6] -= -Math.Sin(x[9] / x[4] + x[7]) + (x[0] - x[5]) * Math.Cos(x[3]) + x[2] - 0.738043007620279801;
            F[7] = -(x[4] * Math.Pow(x[0] - 2 * x[5], 2) - 2 * Math.Sin(-x[8] + x[2]) + 1.5 * x[3] - Math.Exp(x[1] * x[6] + x[9]) + 3.566832198969380904);
            F[8] = -(7 / x[5] + Math.Exp(x[4] + x[3]) - 2 * x[1] * x[7] * x[9] * x[6] + 3 * x[8] - 3 * x[0] - 8.439473450838325749);
            F[9] = -(x[9] * x[0] + x[8] * x[1] - x[7] * x[2] + Math.Sin(x[3] + x[4] + x[5]) * x[6] - 0.7823809523809523809);
            return F;
        }


        static void MetN(ref double[] x,int n,double Eps,ref int iterrations)
        {
            double norm;

            do
            {
                double[] deltax = new double[n];
                double[,] cg = Yakobi(x,n);
                double[] F = Func(x,n);
                


                for (int i = 0; i < n; i++)
                    deltax[i] = 0;


                LU(ref deltax, cg, F, n);


                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    sum = x[i] + deltax[i];
                    x[i] = sum;
                }


                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);

                iterrations++;

            } while (norm>Eps);
        }




        static void MetNM(ref double[] x, int n, double Eps, ref int iterrations)
        {
            double[,] cg = Yakobi(x, n);
            double norm;
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            int[,] P = new int[n, n];
            for (int i = 0; i < n; i++)
                P[i, i] = 1;
            LUm(cg, ref L, ref U, ref P, n);
            do
            {
                double[] deltax = new double[n];
                double[] F = Func(x, n);




                LU2(L,U,P,ref deltax, cg, F, n);
                for (int i = 0; i < n; i++)
                    x[i] += deltax[i];


                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);
                iterrations++;

            } while (norm > Eps);
        }


        static void MetNMC(ref double[] x, int n, double Eps,int count, ref int iterrations)
        {
            double[,] cg = new double[n, n];
            double norm;
            int itt = 0;
            do
            {
                itt++;
                double[] deltax = new double[n];
                double[] F = Func(x, n);

                cg = Yakobi(x, n);

                LU(ref deltax, cg, F, n);
                for (int i = 0; i < n; i++)
                    x[i] += deltax[i];



                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);
                iterrations++;

            } while ((norm > Eps) & (itt<count));

            cg = Yakobi(x, n);

            do
            {
                double[] deltax = new double[n];
                double[] F = Func(x,n);


                LU(ref deltax, cg, F, n);
                for (int i = 0; i < n; i++)
                    x[i] += deltax[i];


                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);
                iterrations++;

            } while (norm > Eps);
        }




        static void MetNComb(ref double[] x, int n, double Eps,int count, ref int iterrations)
        {
            double norm;
            double[,] cg = new double[n, n];

            do
            {
                if ((iterrations/count)%2==0)
                {
                    double[] deltax = new double[n];
                    double[] F = Func(x, n);
                    cg = Yakobi(x, n);

                    for (int i = 0; i < n; i++)
                        deltax[i] = 0;


                    LU(ref deltax, cg, F, n);


                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        sum = x[i] + deltax[i];
                        x[i] = sum;
                    }


                    norm = 0;
                    for (int i = 0; i < n; i++)
                        if (norm < Math.Abs(deltax[i]))
                            norm = Math.Abs(deltax[i]);

                    iterrations++;
                }
                else if (iterrations%count==0)
                {
                    double[] deltax = new double[n];
                    double[] F = Func(x, n);
                    cg = Yakobi(x, n);


                    for (int i = 0; i < n; i++)
                        deltax[i] = 0;


                    LU(ref deltax, cg, F, n);


                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        sum = x[i] + deltax[i];
                        x[i] = sum;
                    }


                    norm = 0;
                    for (int i = 0; i < n; i++)
                        if (norm < Math.Abs(deltax[i]))
                            norm = Math.Abs(deltax[i]);

                    iterrations++;
                }
                else
                {
                    double[] deltax = new double[n];
                    double[] F = Func(x, n);


                    for (int i = 0; i < n; i++)
                        deltax[i] = 0;


                    LU(ref deltax, cg, F, n);


                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        sum = x[i] + deltax[i];
                        x[i] = sum;
                    }


                    norm = 0;
                    for (int i = 0; i < n; i++)
                        if (norm < Math.Abs(deltax[i]))
                            norm = Math.Abs(deltax[i]);

                    iterrations++;
                }
            } while (norm > Eps);
        }


        static void MetNMC2(ref double[] x, int n, double Eps, int count, ref int iterrations, int iterations2)
        {
            double[,] cg = new double[n, n];
            double norm;
            int itt = 0;
            do
            {
                itt++;
                double[] deltax = new double[n];
                double[] F = Func(x, n);

                cg = Yakobi(x, n);

                LU(ref deltax, cg, F, n);
                for (int i = 0; i < n; i++)
                    x[i] += deltax[i];



                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);
                iterrations++;

            } while ((norm > Eps) & (itt < count));

            cg = Yakobi(x, n);

            do
            {
                double[] deltax = new double[n];
                double[] F = Func(x, n);


                LU(ref deltax, cg, F, n);
                for (int i = 0; i < n; i++)
                    x[i] += deltax[i];


                norm = 0;
                for (int i = 0; i < n; i++)
                    if (norm < Math.Abs(deltax[i]))
                        norm = Math.Abs(deltax[i]);
                iterrations++;

            } while ((norm > Eps) & (iterrations<iterations2));
        }


        static void Main(string[] args)
        {
            int iterations = 0;
            int n = 10;
            double[] x=new double[] { 0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 };
           


            Stopwatch t1 = new Stopwatch();
            t1.Start();
            MetN(ref x, n, 0.000001,ref iterations);
            t1.Stop();

            Console.WriteLine("Обычный");
            Console.WriteLine("x");
            for (int i = 0; i < n; i++)
                Console.WriteLine(x[i]);

            double[] F = Func(x,n);

            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(F[i]) > 0.000001)
                    Console.WriteLine("Bad {0}",i);
            }

            Console.WriteLine("Time");
            Console.WriteLine(t1.ElapsedMilliseconds);
            Console.WriteLine("Iterations");
            Console.WriteLine(iterations);



            int iterations1 = 0;
            double[] x1 = new double[] { 0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 };

            t1.Restart();
            MetNM(ref x1, n, 0.000001, ref iterations1);
            t1.Stop();

            Console.WriteLine("Модифицированный");
            Console.WriteLine("x");
            for (int i = 0; i < n; i++)
                Console.WriteLine(x1[i]);

            F = Func(x1, n);

            for (int i = 0; i < n; i++)
            {
                if (Math.Abs(F[i]) > 0.00001)
                    Console.WriteLine("Bad {0}", i);
            }


            Console.WriteLine("Time");
            Console.WriteLine(t1.ElapsedMilliseconds);
            Console.WriteLine("Iterations");
            Console.WriteLine(iterations1);
            Console.WriteLine();

            //int iterations2 = 0;



            Console.WriteLine();
            Console . WriteLine("Пункт с");
            Console.WriteLine();
            for (int k=1;k<iterations;k++)
            {
                Console.Write("K ");
                Console.WriteLine(k);
                double[] x2 = new double[] { 0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 };
                int iterations2 = 0;
                t1.Restart();
                MetNMC(ref x2, n, 0.000001, k, ref iterations2);
                t1.Stop();

                Console.WriteLine("x");
                for (int i = 0; i < n; i++)
                    Console.WriteLine(x2[i]);

                F = Func(x2, n);

                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(F[i]) > 0.000001)
                        Console.WriteLine("Bad {0}", i);
                }

                Console.WriteLine("Time");
                Console.WriteLine(t1.ElapsedMilliseconds);
                Console.WriteLine("Iterations");
                Console.WriteLine(iterations2);
                Console.WriteLine();
            }


            Console.WriteLine();
            Console.WriteLine("Пункт д");
            Console.WriteLine();
            for (int k = 1; k < iterations; k++)
            {
                Console.Write("K ");
                Console.WriteLine(k);
                double[] x2 = new double[] { 0.5, 0.5, 1.5, -1.0, -0.5, 1.5, 0.5, -0.5, 1.5, -1.5 };
                int iterations2 = 0;
                t1.Restart();
                MetNComb(ref x2, n, 0.000001, k, ref iterations2);
                t1.Stop();

                Console.WriteLine("x");
                for (int i = 0; i < n; i++)
                    Console.WriteLine(x2[i]);

                F = Func(x2, n);

                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(F[i]) > 0.000001)
                        Console.WriteLine("Bad {0}", i);
                }

                Console.WriteLine("Time");
                Console.WriteLine(t1.ElapsedMilliseconds);
                Console.WriteLine("Iterations");
                Console.WriteLine(iterations2);
                Console.WriteLine();
            }



            x = new[]{ 0.5, 0.5, 1.5, -1.0, -0.2, 1.5, 0.5, -0.5, 1.5, -1.5 };
            int iterations3 = 0;

            Console.WriteLine();
            Console.WriteLine("Пункт е");
            Console.WriteLine();

            MetN(ref x, n, 0.000001, ref iterations3);
            Console.WriteLine("x");
            for (int i = 0; i < n; i++)
                Console.WriteLine(x[i]);
            Console.WriteLine("Iterations");
            Console.WriteLine(iterations3);

            F = Func(x, n);

            for (int i=0;i<n;i++)
            {
                if (Math.Abs(F[i]) > 0.000001)
                    Console.WriteLine("Bad {0}",i);
            }

            Console.WriteLine();
            x = new[] { 0.5, 0.5, 1.5, -1.0, -0.2, 1.5, 0.5, -0.5, 1.5, -1.5 };
            int iterations4 = 0;

            Console.WriteLine("Модифицированный метод");
            MetNM(ref x, n, 0.000001,ref iterations4);
            Console.WriteLine("x");
            for (int i = 0; i < n; i++)
                Console.WriteLine(x[i]);

            Console.WriteLine();
            Console.WriteLine("Модифицированный+обычный метод");
            for (int k = 1; k < iterations3; k++)
            {
                Console.Write("K ");
                Console.WriteLine(k);
                double[] x2 = new double[] { 0.5, 0.5, 1.5, -1.0, -0.2, 1.5, 0.5, -0.5, 1.5, -1.5 };
                int iterations2 = 0;
                MetNMC2(ref x2, n, 0.000001, k, ref iterations2,iterations4);

                Console.WriteLine("x");
                for (int i = 0; i < n; i++)
                    Console.WriteLine(x2[i]);

                F = Func(x2, n);

                for (int i = 0; i < n; i++)
                {
                    if (Math.Abs(F[i]) > 0.000001)
                        Console.WriteLine("Bad {0}", i);
                }

                Console.WriteLine("Iterations");
                Console.WriteLine(iterations2);
                Console.WriteLine();
            }

            Console.ReadKey();
        }
    }
}
