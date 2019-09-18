using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace itecmethod
{
    class Program
    {

        static void LUm(double[,] A, ref double[,] L, ref double[,] U, ref int[,] P, int n)
        {
            double max = Math.Abs(A[0, 0]);
            int q = 0;


            for (int i = 0; i < n; i++)
                P[i, i] = 1;

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

        static void Yacoby(double[,] A, double[] b,ref double[] x,int n,ref int itYacob,int diap,double Eps)
        {



        //Find norm B
            double normB= 0;
            for (int i=0;i<n;i++)
            {
                double sum = 0;
                for (int j=0;j<n;j++)
                {
                    if (i!=j)
                        sum += Math.Abs(A[i, j]/A[i,i]);
                }
                if (normB<sum)
                    normB = sum;
            }



            //Find Epsilon
            double Epsilon = (1 - normB) / normB * Eps;

            Random rand = new Random();
            //Find x^(1)
            for (int i = 0; i < n; i++)
            {
                x[i] = b[i] / A[i, i];
            }



            double[] x1 = new double[n];

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (i!=j)
                        sum -= A[i, j] * x[j];
                }

                sum += b[i];
                sum /= A[i, i];
                x1[i] = sum;
            }


            //Find ||x-x^(1)||
            double normX = 0;
            for (int i=0;i<n;i++)
            {
                if (normX < Math.Abs(x1[i] - x[i]))
                    normX = Math.Abs(x1[i] - x[i]);
            }

            itYacob++;

            while (normX>Epsilon)
            {
                //Find x^(1)
                for (int i = 0; i < n; i++)
                    x[i] = x1[i];

                for (int i = 0; i < n; i++)
                {
                    x1[i] = 0;
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                            sum -= A[i, j] * x[j];
                    }

                    sum += b[i];
                    sum /= A[i, i];
                    x1[i] = sum;
                }


                //Find ||x-x^(q)||
                normX = 0;
                for (int i = 0; i < n; i++)
                {
                    if (normX < Math.Abs(x1[i] - x[i]))
                        normX = Math.Abs(x1[i] - x[i]);
                }

                itYacob++;
            }


            double normAxb = 0;

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x1[j];
                }
                sum = Math.Abs(sum - b[i]);
                if (sum > normAxb)
                   normAxb=sum;   
            }

            while (normAxb>Eps)
            {
                //Find x^(1)
                for (int i = 0; i < n; i++)
                    x[i] = x1[i];

                for (int i = 0; i < n; i++)
                {
                    x1[i] = 0;
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        if (i != j)
                            sum -= A[i, j] * x[j];
                    }

                    sum += b[i];
                    sum /= A[i, i];
                    x1[i] = sum;
                }

                normAxb = 0;

                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += A[i, j] * x1[j];
                    }
                    sum = Math.Abs(sum - b[i]);
                    if (sum > normAxb)
                        normAxb = sum;
                }

                itYacob++;
            }

            for (int i = 0; i < n; i++)
                x[i] = x1[i];

        }


        static void Seidel(double[,] A,double [] b, ref double[] x,int n,ref int count,double Eps)
        { 
            //Find norm B
            double normB = 0;
           for (int i = 0; i < n; i++)
           {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                        sum += Math.Abs(A[i, j] / A[i, i]);
                }
                if (normB < sum)
                    normB = sum;
           }



                //Find Epsilon
                double Epsilon = (1 - normB) / normB * Eps;


                //Find x^(1)
                for (int i = 0; i < n; i++)
                {
                    x[i] = b[i] / A[i, i];
                }



                double[] x1 = new double[n];
            for (int i = 0; i < n; i++)
                x1[i] = x[i];

                for (int i = 0; i < n; i++)
                {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum -= A[i, j] * x1[j];
                for (int j = i + 1; j < n; j++)
                    sum -= A[i, j] * x[j];
                sum += b[i];
                sum /= A[i, i];
                x1[i] =sum;
                }

                //Find ||x-x^(1)||
                double normX = 0;
                for (int i = 0; i < n; i++)
                {
                    if (normX < Math.Abs(x1[i] - x[i]))
                        normX = Math.Abs(x1[i] - x[i]);
                }

                count++;

                while (normX > Epsilon)
                {
                //Find x^(q)
                for (int i = 0; i < n; i++)
                    x[i] = x1[i];
                    
                    for (int i = 0; i < n; i++)
                    {
                        double sum = 0;
                        for (int j = 0; j < i; j++)
                            sum -= A[i, j] * x1[j];
                        for (int j = i + 1; j < n; j++)
                            sum -= A[i, j] * x[j];
                        sum += b[i];
                        sum /= A[i, i];
                        x1[i] = sum ;
                    }


                    //Find ||x-x^(q)||
                    normX = 0;
                    for (int i = 0; i < n; i++)
                    {
                        if (normX < Math.Abs(x1[i] - x[i]))
                            normX = Math.Abs(x1[i] - x[i]);
                    }

                    count++;
                }

            double normAxb = 0;

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x1[j];
                }
                sum = Math.Abs(sum - b[i]);
                if (sum > normAxb)
                    normAxb = sum;
            }

            while (normAxb > Eps)
            {
                //Find x^(1)
                for (int i = 0; i < n; i++)
                    x[i] = x1[i];

                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < i; j++)
                        sum -= A[i, j] * x1[j];
                    for (int j = i + 1; j < n; j++)
                        sum -= A[i, j] * x[j];
                    sum += b[i];
                    sum /= A[i, i];
                    x1[i] = sum;
                }

                normAxb = 0;

                for (int i = 0; i < n; i++)
                {
                    double sum = 0;
                    for (int j = 0; j < n; j++)
                    {
                        sum += A[i, j] * x1[j];
                    }
                    sum = Math.Abs(sum - b[i]);
                    if (sum > normAxb)
                        normAxb = sum;
                }

                count++;
            }

            for (int i = 0; i < n; i++)
                x[i] = x1[i];

        }



        static double findnormC(double[,] A,int n)
        {
            double norm = 0;
            for (int i=0;i<n;i++)
            {
                double sum = 0;
                for (int j=0;j<n;j++)
                    if (i != j)
                        sum += Math.Abs(A[i, j] / A[i, i]);
                if (norm < sum)
                    norm = sum;
            }
            return norm;
        }

        static double findnormD(double[,] A,double[]b,int n)
        {
            double norm = 0;
            for (int i=0;i<n;i++)
            {
                if (norm < Math.Abs(b[i] / A[i, i]))
                    norm = Math.Abs(b[i] / A[i, i]);
            }
            return norm;
        }









        static void Main(string[] args)
        {
            Random rand = new Random();
            int n = Convert.ToInt32(Console.ReadLine());
            int diap = Convert.ToInt32(Console.ReadLine());
            int q = Convert.ToInt32(Console.ReadLine());
            double[,] A = new double[n, n];
            double[,] B = new double[n, n];
            double[,] L = new double[n, n];
            double[,] U = new double[n, n];
            int[,] P = new int[n, n];


            double[] b = new double[n];
            double[] d = new double[n];
            double[] x = new double[n];
            double[] x1 = new double[n];
            double[] x2 = new double[n];
            double[] Pb = new double[n];
            double[] y = new double[n];

            double Eps = 0.000001;

            for (int i=0;i<n;i++)
            {
                for (int j=0;j<n;j++)
                {
                    if (i!=j)
                    {
                        A[i, j] = rand.Next(-diap/n,diap/n);
                    }
                }
            }

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                    if (i != j)
                        sum += Math.Abs(A[i, j]);
                A[i, i] = rand.Next(Convert.ToInt32(sum)*q , diap*q);
                b[i] = rand.Next(-diap*q, diap*q);
                /*int r = rand.Next(2);
                if (r == 1)
                    A[i, i] = -A[i, i];*/
            }

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(A[i, j]);
                    Console.Write(" ");
                }
                Console.WriteLine();
            }
            Console.WriteLine();


            Console.WriteLine("b");
            for (int i = 0; i < n; i++)
                Console.WriteLine(b[i]);

            Console.WriteLine();

            LUm(A,ref L,ref U,ref P, n);

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
                    sum -= U[i, j] * x2[j];
                sum /= U[i, i];
                x2[i] = sum;
            }

            Console.WriteLine("LU result");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(x2[i]);
            }
            Console.WriteLine();




            int itYacob=0;
            Yacoby(A, b,ref x, n,ref itYacob,diap,Eps);

            Console.WriteLine("Yacoby");
            Console.WriteLine("X");

            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(x[i]);
            }
            Console.Write("Number of iterations ");
            Console.WriteLine(itYacob);
            Console.WriteLine();
            Console.WriteLine("Result");
            for (int i=0;i<n;i++)
            {
                double sum = 0;
                for (int j=0;j<n;j++)
                {
                    sum += A[i, j] * x[j];
                }
                Console.WriteLine(sum);
            }
            Console.WriteLine();
            Console.WriteLine("|Ax-b|");
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x[j];
                }
                sum = Math.Abs(sum- b[i]);
                Console.WriteLine(sum);
            }
            Console.WriteLine();
            Console.WriteLine("LU-Yacoby");
            for (int i=0;i<n;i++)
            {
                Console.WriteLine(Math.Abs(x2[i] - x[i]));
            }
            Console.WriteLine();

            Console.WriteLine();
            Console.WriteLine("Seidel");
            int itSeid = 0;
            Seidel(A, b, ref x1, n,ref itSeid,Eps);
            Console.WriteLine("X");

            for (int i=0;i<n;i++)
            {
                Console.WriteLine(x1[i]);
            }
            Console.WriteLine();
            Console.Write("Number of iterations ");
            Console.WriteLine(itSeid);
            Console.WriteLine();
            Console.WriteLine("Result");
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x1[j];
                }
                Console.WriteLine(sum);
            }
            Console.WriteLine();
            Console.WriteLine("|Ax-b|");
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    sum += A[i, j] * x1[j];
                }
                sum = Math.Abs(sum - b[i]);
                Console.WriteLine(sum);
            }
            Console.WriteLine();
            Console.WriteLine("LU-Seidel");
            for (int i = 0; i < n; i++)
            {
                Console.WriteLine(Math.Abs(x2[i] - x1[i]));
            }
            Console.WriteLine();

            double normC = findnormC(A, n);
            double normd = findnormD(A, b, n);


            double sas = (1 - normC)*Eps;

            double apprior = Math.Log( sas / normd);
             apprior/= Math.Log(normC);


            //Console.WriteLine(apprior);

            Console.ReadKey();
        }
    }
}
