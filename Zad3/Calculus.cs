using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Zad3
{
    class Calculus
    {
        const double m0 = 9.1e-31;
        const double hb = 1.05459e-34;
        protected int _n;
        protected double _k;

        public Calculus(int n, double k)
        {
            _n = n;
            _k = k;
        }
        public void Normalisation(double[] vector)
        {
            double nor = 0;
            for (int i = 0; i <_n; i++)
            {
                nor += vector[i] * vector[i];
            }
            nor = Math.Sqrt(nor);
            for (int i = 0; i < _n; i++)
            {
                vector[i] = vector[i] / nor;
            }
        }
        public double Factorial(int ent)
        {
            if (ent == 0) return 1.0;
            if (ent == 1) return 1.0;
            return Convert.ToDouble(ent) * Factorial(ent - 1);
        }
        public double HermitePolynominal(int n, double x)
        {
            if (n == 0) return 1.0;
            if (n == 1) return 2.0 * x;
            return 2.0 * x * HermitePolynominal(n - 1, x) - 2 * (n - 1) * HermitePolynominal(n - 2, x);
        }
        public double AnaliticalVector(int n, double x, double K)
        {
            double m = 9.109E-31;
            double hb = 1.055E-34;
            x = 2.0 * x;

            double lOmega = Math.Sqrt(K / m);
            double tA = 1.0 / Math.Sqrt(Math.Pow(2.0, n) * Factorial(n));
            double tB = Math.Pow((m * lOmega / Math.PI / hb), 0.25);
            double tC = Math.Exp(-(m * lOmega * x * x) / (2.0 * hb));
            double tPsi = tA * tB * tC * HermitePolynominal(n, Math.Sqrt(m * lOmega / hb) * x);
            return tPsi;
        }
        public double[] GenerateAnalitycalEigenVector(int n, double Width, int nz, double K)
        {
            double[] eigvec = new double[nz];
            double dZ = Width / Convert.ToDouble(nz - 1);
            for (int i = 0; i < nz; i++)
            {
                eigvec[i] = AnaliticalVector(n, Convert.ToDouble(i) * dZ - Width,K);
            }
            Normalisation(eigvec);
            return eigvec;
        }
        public double AnalitycalEnergies(double ent)
        {
            double energy;
            double omega = Math.Sqrt(_k / m0);
            energy = hb * omega * (ent + 0.5);
            return energy;
        }
        public int PositionOfMinimum(double[]vector)
        {
            int pos = 0;
            double minVec = vector[0];
            for (int i = 0; i <_n; i++)
            {
                if (vector[i] < minVec)
                {
                    minVec = vector[i];
                    pos = i;
                }
            }
            return pos;
        }
        public int NegativeU(double z, double[]e, double[]d)
        {
            int negativeUCount = 0;
            double u = d[0] - z;
            if (u < 0.0) negativeUCount++;
            for (int i = 1; i < d.Length; i++)
            {
                u = d[i] - z - Math.Pow(e[i - 1], 2) / u;
                if (u < 0.0) negativeUCount++;
            }
            return negativeUCount;
        }
        public double CalculateEigenValue(double[] e, double[] d, int M, double prec)
        {
            double minZ = -1.0;
            while (NegativeU(minZ, e, d) > M)
            {
                minZ *= 2.0;
            }
            double maxZ = 1.0;
            while (NegativeU(maxZ,e,d)<=M)
            {
                maxZ *= 2.0;
            }
            while ((maxZ - minZ) > prec)
            {
                double midZ = (maxZ + minZ) / 2.0;
                if (NegativeU(midZ, e, d) <= M) minZ = midZ;
                else maxZ = midZ;
            }
            return (maxZ + minZ) / 2.0;
        }
        public double[] CalculatedEigenVector(double[] e, double[] d,double nEigenValue,int pivpos)
        {
            int nz = d.Length;
            double[] eigvec = new double[nz];
            for (int i = 0; i < nz; i++)
            {
                eigvec[i] = 0.0;
            } //zeroing eigvec
            double[] l = new double[nz];
            double[] r = new double[nz];

            l[0] = 1.0 / (d[0] - nEigenValue);
            for (int j = 1; j < nz; j++)
            {
                l[j] = 1.0 / (d[j] - nEigenValue - Math.Pow(e[j], 2) * l[j - 1]);
            }

            r[nz - 1] = 1.0 / (d[nz - 1] - nEigenValue);
            for (int i = nz-2; i >=0 ; i--)
            {
                r[i] = 1.0 / (d[i]-nEigenValue-Math.Pow(e[i+1],2)*r[i+1]);
            }

            eigvec[pivpos] = 10.0;
            for (int j = pivpos+1; j < nz; j++)
            {
                eigvec[j] = -e[j] * r[j] * eigvec[j - 1];
            }
            for (int j = pivpos-1; j >= 0; j--)
            {
                eigvec[j] = -e[j +1] * l[j] * eigvec[j + 1];
            }
            Normalisation(eigvec);
            return eigvec;
        }

    }
}
