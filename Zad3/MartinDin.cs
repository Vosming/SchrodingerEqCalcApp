using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Zad3
{
    public class MartinDin
    {
        const double q = 1.0e-19;
        const double m0 = 9.1e-31;
        const double hb = 1.05459e-34;
        const int nz = 1000;
        const int number = 10;

        protected double _k;
        protected double[,] envelope_functions = new double[nz, number];// 2 wymiarowa tablica przechowujaca funkcje wlasne
        protected double[] energies = new double[nz];// tablica przechowujaca obliczone energie

        public double Energies (int n)
        {
            double energy = energies[n];
            return energy;
        }

        double[] m = new double[nz];
        double[] v = new double[nz];
        // for loop for creating m and v

        public MartinDin(double k)
        {
            _k = k;
            for (int i = 0; i < nz; i++)
            {
                m[i] = 1.0 * m0;// masa w kg
                x = -5.0 + 10.0 * (float)i / nz;
                v[i] = _k * x * x * q;
                 

            }
            double omega = Math.Sqrt(_k / m0);

        }

        int i, km, n, nu, i0, k, ii;
        double[,] d = new double[nz, 2];
        double[] enn = new double[number];
        double x0, x1, xk, dx, xroot, xw, y1, y2, dok, w1, w2, delta, xmx, xmy, g0, x;
        double[] del_m = new double[nz];
        double[] del_p = new double[nz];
        double[] psi_dwsz = new double[nz];
        double[] pr = new double[nz];
        double energia, calka;
        double potmin, potmax, max_psi;
        double delta_z = (10.0 / (double)nz)*1e-9;

        public void Normalisation(int ent,double[]vector)
        {
            double nor = 0;
            for (int i = 0; i < ent; i++)
            {
                nor += vector[i] * vector[i];
            }
            nor = Math.Sqrt(nor);
            for (int i = 0; i < ent; i++)
            {
                v[i] = v[i] / nor;
            }
        }
        public double AnalitycalEigenValue(double ix, double ent)
        {
            double omega = Math.Sqrt(_k / m0);
            double eigenValue;
            double factorial=1;
            for (int i = 2; i <= ent; i++)
            {
                factorial = factorial * n;
            }
            eigenValue = (1 / (Math.Sqrt(Math.Pow(2, ent) * factorial))) * Math.Pow((m0 * omega) / (Math.PI * hb), 0.25)
                *Math.Exp(-(m0*omega*Math.Pow(ix,2))/(2*hb))
                *alglib.hermitecalculate(Convert.ToInt32(ent),ix)
                *(Math.Sqrt((m0*omega)/hb)*ix);
            return eigenValue;
        }
        public double []GenerateAnalitycalEiggenVector(int n, double Width, int nz, double K)
        {
            double[] eigvec = new double[nz];
            double dZ = Width / Convert.ToDouble(nz - 1);
            for (int i = 0; i < nz; i++)
            {
                eigvec[i] = AnalitycalEigenValue(n, Convert.ToDouble(i) * dZ - Width);
            }
            Normalisation(n, eigvec);
            return eigvec;
        }
        public double AnalitycalEnergies(double ent)
        {
            double energy;
            double omega = Math.Sqrt(_k / m0);
            energy = hb * omega*(ent + 0.5);
            return energy;
        }
            





        double y=0;
        public double CalcY(double x)
        {
                
            double i = -5.0 + (10.0 * x / nz);//dysktetyzacja
            y = _k * (i * i);//obliczenie potencjau w eV

            return y;
        }

        public double CalcN(int x, int n)
        {
            Calculate();
            double y = envelope_functions[x, n]*10;
            return y;
        }

        public double Energy(int i)
        {
            double energy = energies[i] / energies[0];
            return energy;
        }
        
        public void Calculate()
        {
            g0 = 2.0 * m0 / hb / hb * (1.0e-3 * q) * 1.0e-20;



            max_psi = 100.0;
            dok = 1.0e-24;
            delta = delta_z / 1.0e-10;


            potmin = v[0];
            potmax = v[0];
            k = nz / 2;

            for (i = 0; i < nz; i++)
            {
                if (v[i] < potmin)
                {
                    potmin = v[i];
                    k = i;
                }
                if (v[i] > potmax)
                {
                    potmax = v[i];
                }
            }

            potmin = (potmin) * delta * delta * g0 / q * 1.0e3;
            potmax = (potmax) * delta * delta * g0 / q * 1.0e3;

            xmx = potmin;
            xmy = potmax;
            x0 = xmx;
            x1 = xmy;
            km = nz - 1;
            x = 1.0;


            for (i = 0; i < nz - 1; i++)
            {

                w1 = 1.0 / m[i] * m0;
                w2 = 1.0 / m[i + 1] * m0;
                d[i, 1] = w1 + w2;
                d[i, 1] = d[i, 1] + delta * delta * g0 * v[i] / q * 1.0e3;
                d[i + 1, 0] = -w2;

            }

            d[nz - 1, 1] = d[nz - 2, 1];
            d[nz - 1, 0] = d[nz - 2, 0];




            for (i = 1; i < number + 1; i++)
            {

                nu = 1;
            l7401:
                xk = 0.0;
                nu = nu + 1;
                if (nu > 100) { goto l7513; }
                xw = (double)i;
                dx = x1 - x0;
                if (dx < dok) { goto l7513; }
                x = (x0 + x1) / 2.0;
                y1 = d[0, 1] - x;
                if (y1 > 0.0) { goto l7133; }
                xk = xk + 1.0;
                if (xk >= xw) { goto l7403; }
            l7133:

                for (n = 1; n < km; n++)
                {
                    y2 = d[n, 1] - x - d[n, 0] * d[n, 0] / y1;
                    if (y2 > 0.0) { goto l7200; }
                    xk = xk + 1.0;
                    if (xk >= xw) { goto l7403; }
                l7200:
                    y1 = y2;
                }


                y2 = d[nz - 1, 1] - x - d[nz - 1, 0] * d[nz - 1, 0] / y1;
                if (y2 > 0.0) { goto l7202; }
                xk = xk + 1.0;
                if (xk > xw) { goto l7403; }
            l7202:
                dx = x1 - x0;
                if (dx < dok) { goto l7513; }
                x0 = x;
                goto l7401;
            l7403:
                dx = x1 - x0;
                if (dx < dok) { goto l7513; }
                x1 = x;
                goto l7401;
            l7513:

                xroot = (x0 + x) / 2.0;
                pr[i - 1] = xroot;
                enn[i - 1] = xroot / delta / delta / g0;


                energies[i - 1] = enn[i - 1] * q * 1.0e-3;


                x1 = xmy;
                x0 = xroot;
            }



            for (ii = 0; ii < number; ii++)
            {

                energia = enn[ii] * delta * delta * g0;
                max_psi = 1.0;

                for (i = 0; i < nz; i++)
                {
                    psi_dwsz[i] = 0.0;
                }


            l1000:
                i0 = k;
                psi_dwsz[i0] = 10.0;
                del_m[0] = 1.0 / (d[0, 1] - energia);
                del_p[nz - 1] = 1.0 / (d[nz - 1, 1] - energia);
                for (i = 1; i < nz; i++)
                {
                    del_m[i] = 1.0 / (d[i, 1] - energia - d[i, 0] * d[i, 0] * del_m[i - 1]);
                }

                for (i = 1; i < nz - 1; i++)
                {

                    del_p[nz - i - 1] = 1.0 / (d[nz - i - 1, 1] - energia - d[nz - i, 0] * d[nz - i, 0] * del_p[nz - i]);

                }
                for (i = i0; i < nz - 1; i++)
                {
                    psi_dwsz[i + 1] = -d[i + 1, 0] * del_p[i + 1] * psi_dwsz[i];

                    if (psi_dwsz[i + 1] > max_psi)
                    {
                        max_psi = psi_dwsz[i - 1];
                        k = i + 1;
                    }

                }


                for (i = i0; i >= 1; i--)
                {
                    psi_dwsz[i - 1] = -d[i, 0] * del_m[i - 1] * psi_dwsz[i];

                }

                if (i0 != k) { goto l1000; }


                calka = 0.0e0;
                for (i = 0; i < nz; i++)
                {
                    calka = calka + psi_dwsz[i] * psi_dwsz[i];

                }

                for (i = 0; i < nz; i++)
                {
                    psi_dwsz[i] = psi_dwsz[i] / System.Math.Sqrt(calka);
                    envelope_functions[i, ii] = psi_dwsz[i];
                }
            }
            
        }


    }
}
