using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace LU
{
    class Program
    {
        static void LUm(double[,] A, ref double[,] L,ref double [,] U, ref int[,] P, int n,ref int count)
        {
            double max = Math.Abs(A[0, 0]);
            int q = 0;


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
                count++;
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
                    count++;
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

        

        static void LUP (ref double [,]  A,ref double[,] LU,ref int [,] P,int n,ref int count)
        {
            double max = Math.Abs(A[0, 0]);
            int p = 0, q=0 ;


            for (int j=1;j<n;j++)
            {
                if (Math.Abs(A[j, 0]) > max)
                {
                    max = Math.Abs(A[j, 0]);
                    q = j;
                }
            }
            if (q != 0)
            {
                count++;
                for (int i = 0; i < n; i++)
                {
                    double c;
                    int d;
                    d = P[p, i];
                    P[p, i] = P[q, i];
                    P[q, i] = d;
                    c = A[p, i];
                    A[p, i] = A[q, i];
                    A[q, i] = c;
                    
                }
            }

            LU[0, 0] = A[0, 0];
            for (int j=1;j<n;j++)
            {
                LU[0, j] = A[0, j];
                LU[j, 0] = A[j, 0] / LU[0, 0];
            }
            for (int i=1;i<n;i++)
            {
                max = Math.Abs(A[i, i]);
                p = i;q = p;
                for (int j = i+1; j < n; j++)
                {
                    if (Math.Abs(A[j, i]) > max)
                    {
                        max = Math.Abs(A[j, i]);
                        q = j;
                    }
                }
                if (q != p)
                {
                    count++;
                    for (int j = 0; j < n; j++)
                    {
                        double c;
                        int d;
                        double v;
                        d = P[p, j];
                        P[p, j] = P[q, j];
                        P[q, j] = d;
                        c = A[p, j];
                        A[p, j] = A[q, j];
                        A[q, j] = c;
                        v = LU[p, j];
                        LU[p, j] = LU[q, j];
                        LU[q, j] = v;
                    }
                }

                LU[i, i] = A[i, i];
                for (int k = 0; k <= i - 1; k++)
                    LU[i, i] -= LU[i, k] * LU[k, i];
                for (int j=i+1;j<n;j++)
                {
                    LU[i, j] = A[i, j];
                    LU[j, i] = A[j, i];
                    for (int k = 0; k <= i - 1; k++)
                    {
                        LU[i, j] -= LU[i, k] * LU[k, j];
                        LU[j, i] -= LU[j, k] * LU[k, i];
                    }
                    LU[j, i] /= LU[i, i];
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
            {
                P[i, i] = 1;
                for (int j = 0; j < n; j++)
                    A1[i, j] = A[i, j];
            }

            int count=0;

            LUm(A1, ref L, ref U, ref P, n,ref count);


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




        static void Main(string[] args)
        {
            Random rand = new Random();
            int n = Convert.ToInt32(Console.ReadLine());
            int [] diap=new int [2];
            for (int i=0;i<2;i++)
                diap[i]= Convert.ToInt32(Console.ReadLine());
            double[,] A = new double[n, n];
            double[,] A2 = new double[n, n];
            double[,] PA = new double[n, n];
            double[,] A4 = new double[n, n];
            double[,] A1 = new double[n, n];
            double[,] Aobr = new double[n, n];

            double[] b = new double [n];
            double[] Pb = new double[n];

            double[] x = new double[n];
            double[] y = new double[ n];
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            double[,] L1 = new double[n, n];
            double[,] U1 = new double[n, n];
            double[,] Ucopy = new double[n, n];
            int[,] P = new int[n, n];
            int[,] P1 = new int[n, n];

            for (int i = 0; i < n; i++)
            {
                P[i, i] = 1;
                P1[i, i] = 1;
            }
            for (int i=0;i<n;i++)
            {
                for (int j = 0; j < n; j++)
                {
                    A[i, j] = rand.Next(diap[0], diap[1]);
                    A2[i, j] = A[i, j];
                    A4[i, j] = A[i, j];
                    while (A[i, j] == 0)
                    {
                        A[i, j] = rand.Next(diap[0], diap[1]);
                        A2[i, j] = A[i, j];
                        A4[i, j] = A[i, j];
                    }
                
                    }
                b[i] = rand.Next(diap[0], diap[1]);
                while (b[i]==0)
                {
                    b[i]= rand.Next(diap[0], diap[1]);
                }
            }

            /*Console.WriteLine("A");
            for (int i=0;i<n;i++)
            {
                for (int j=0;j<n;j++)
                {
                    A[i, j] = Convert.ToDouble(Console.ReadLine());
                    A2[i, j] = A[i, j];
                    A4[i, j] = A[i, j];
                }
            }


            Console.WriteLine("b");
            for (int i=0;i<n;i++)
            {
                b[i] = Convert.ToDouble(Console.ReadLine());
            }*/


            /*for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(A[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }

            Console.WriteLine();

            int count = 0;
            int count1 = 0;

            //Поиск вывод LU L U
            LUP(ref A, ref LU, ref P, n,ref count);
            LUm(ref A4, ref L, ref U, ref P1, n,ref count1);
            if (count == count1)
                Console.WriteLine(true);

            Console.WriteLine("P");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(P[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            

            for (int i=0;i<n;i++)
            {
                L[i, i] = 1;
            }

            Console.WriteLine("LU");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(LU[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(L[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(U[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();


            Console.WriteLine("PA");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {

                        PA[i, j] += P[i, k] * A2[k, j];
                    }
                    Console.Write(PA[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            Console.WriteLine("LU");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i>j)
                    {
                        for (int k = 0; k <= j; k++)
                            A1[i, j] += LU[i, k] * LU[k, j];
                    }
                    else
                    {
                        for (int k = 0; k < i; k++)
                            A1[i, j] += LU[i, k] * LU[k, j];
                        A1[i, j] += LU[i, j];

                    }
                }
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(A1[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            

            for (int i = 0; i < n; i++)
                for (int k = 0; k < n; k++)
                    Pb[i] += P[i, k] * b[k];

            Console.WriteLine("b");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(b[i]);
            }
            Console.WriteLine();

            Console.WriteLine("Pb");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(Pb[i]);
            }
            Console.WriteLine();

            for (int i=0;i<n;i++)
            {
                double sum = Pb[i];
                for (int j = 0; j <= i - 1; j++)
                    sum -= L[i, j] * y[j];
                y[i] = sum;
            }

            Console.WriteLine("y");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(y[i]);
            }
            Console.WriteLine();

            for (int i=n-1;i>=0;i--)
            {
                double sum = y[i];
                for (int j = n - 1; j >= i; j--)
                    sum -= U[i, j] * x[j];
                sum /= U[i, i];
                x[i] = sum;
            }

            Console.WriteLine("x");
            for (int i=0;i<n;i++)
            {
                Console.WriteLine(x[i]);
            }
            Console.WriteLine();


            //det A
            double detA = 1;
            for (int i = 0; i < n; i++)
                detA *= U[i, i];

            detA = detA*Math.Pow(-1,count);



            Console.Write("detA= ");
            Console.WriteLine(detA);
            Console.WriteLine();



            //ПРоверка PA=LU
            int q = 1;

            for (int i=0;i<n;i++)
            {
                for (int j = 0; j < n; j++)
                    if (Math.Abs(A1[i, j] - PA[i, j]) > 0.000001)
                    {
                        q = 0;
                        break;
                    }
                if (q == 0)
                    break;
            }
            if (q == 1)
                Console.WriteLine("PA=LU");
            else
                Console.WriteLine("PA!=LU");
            Console.WriteLine();*/


            //Проверка Ax-b=0
            Stopwatch t = new Stopwatch();
            t.Start();
            LU(ref x, A, b, n);
            t.Stop();
            Console.WriteLine(t.ElapsedMilliseconds);

            int q = 1;
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int k = 0; k < n; k++)
                    sum += A2[i, k] * x[k];
                if (Math.Abs(sum-b[i]) > 0.000001)
                {
                    q = 0;
                    break;
                }
            }

            if (q == 1)
                Console.WriteLine("Ax-b=0");
            else
                Console.WriteLine("Ax-b!=0");
            Console.WriteLine();


            /*
            //Поиск L^-1
            for (int i = 0; i < n; i++)
            {
                L1[i, i] = 1;
            }

            for (int k=0;k<n;k++)
                for (int i = k + 1; i < n; i++)
                    for (int j = 0; j <= k; j++)
                        L1[i,j] -= L[i, k] * L1[k, j];

            double[,] E = new double[n, n];

            double[,] E1 = new double[n, n];

            for (int i = 0; i < n; i++)
                E1[i, i] = 1;

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                        E[i, j] += L[i, k] * L1[k, j];
                }
            }


            q = 1;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(E[i, j] - E1[i, j]) > 0.0000001)
                    {
                        q = 0;
                        break;
                    }
                }
                if (q == 0)
                    break;
            }
            if (q == 1)
                Console.WriteLine("NiceL");
            else
                Console.WriteLine("HoustonWehaveProblemL");
            Console.WriteLine();



            //Поиск U^-1
            for (int i = 0; i < n; i++)
                for (int j = i; j < n; j++)
                    Ucopy[i, j] = U[i, j] / U[i, i];

            for (int i = 0; i < n; i++)
                U1[i, i] = 1/U[i,i];

            for (int k = n-1; k >= 0; k--)
                for (int i = k - 1; i >=0 ; i--)
                    for (int j = k; j <n; j++)
                        U1[i, j] -= Ucopy[i, k] * U1[k, j];


            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                        sum += U[i, k] * U1[k, j];
                    E[i, j] = sum;
                }

            q = 1;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(E[i, j] - E1[i, j]) > 0.0000001)
                    {
                        q = 0;
                        break;
                    }
                }
                if (q == 0)
                    break;
            }
            if (q == 1)
                Console.WriteLine("NiceU");
            else
                Console.WriteLine("HoustonWehaveProblemU");
            Console.WriteLine();



            //Поиск A^-1
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n ; j++)
                    for (int k = n-1; k >=i && k>=j ; k--)
                        Aobr[i, j] += U1[i, k] * L1[k, j];

            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                    for (int k = 0; k < n; k++)
                        Aobr2[i, j] += Aobr[i, k] * P[k, j];

            Console.WriteLine("A^-1");
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(Aobr2[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();

            for (int i=0;i<n;i++)
                for (int j=0;j<n;j++)
                {
                    E[i, j] = 0;
                    double sum = 0;
                    for (int k = 0; k < n; k++)
                         sum+= A2[i, k] * Aobr2[k, j];
                    E[i, j] = sum;
                }


            q = 1;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (Math.Abs(E[i, j] - E1[i, j]) > 0.0000001)
                    {
                        q = 0;
                        break;
                    }
                }
                if (q == 0)
                    break;
            }
            if (q == 1)
                Console.WriteLine("A*A^-1=E");
            else
                Console.WriteLine("HoustonWehaveProblemAobr");
            Console.WriteLine();

            double maxA = -1;
            double maxAobr = -1;

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                double sumobr = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += Math.Abs(A[i, j]);
                    sumobr += Math.Abs(Aobr[i, j]);
                }
                if (maxA < sum)
                    maxA = sum;
                if (maxAobr < sumobr)
                    maxAobr = sumobr;
            }

            maxA *= maxAobr;
            Console.Write("cond =");
            Console.WriteLine(maxA);*/
            Console.ReadKey();
        }
    }
}
