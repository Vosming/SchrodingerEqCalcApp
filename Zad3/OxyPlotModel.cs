using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Drawing;
using System.Runtime.CompilerServices;
using System.Windows.Forms;
using System.IO;
namespace Zad3

{
    using System;
    using OxyPlot;
    using OxyPlot.Series;
    using Prism.Mvvm;
    using Prism.Commands;
    using System.Text;


    public class OxyPlotModel : BindableBase
    {
        private OxyPlot.PlotModel plotModel;
        private double k;
        private int n;
        private double mdEnergy;
        private double anEnergy;

        public OxyPlot.PlotModel PlotModel
        {
            get
            {
                return plotModel;
            }
            set
            {
                SetProperty(ref plotModel, value);
            }
        }

        public double  Coefficient_K
        {
            get
            {
                return k;
            }

            set
            {
                SetProperty(ref k, value);
            }
        }

        public int Coefficient_N
        {
            get
            {
                return n;
            }

            set
            {
                SetProperty(ref n, value);
            }
        }
        public double MDEnergy
        {
            get
            {
                return mdEnergy;
            }
            set
            {
                SetProperty(ref mdEnergy, value);
            }
        }
        public double AnEnergy
        {
            get
            {
                return anEnergy;
            }
            set
            {
                SetProperty(ref anEnergy, value);
            }
        }



        public DelegateCommand GenerateChartCommand { get; private set; }
        public DelegateCommand GenerateFileCommand { get; private set; }
        public DelegateCommand GenerateChartNewCommand { get; private set; }

        public OxyPlotModel()
        {
            PlotModel = new PlotModel();
            GenerateChartCommand = new DelegateCommand(GenerateChart);
            GenerateFileCommand = new DelegateCommand(GenerateFile);
            GenerateChartNewCommand = new DelegateCommand(GenerateChartNew);


        }

        private void GenerateChart()
        {
            var md = new MartinDin(Coefficient_K);
            PlotModel.Series.Clear();
            var points = new LineSeries {
                Title="Martin Deen"
            };
            for (int i=-500; i<500;i++)
            {
                points.Points.Add(new DataPoint(i, md.CalcN(i+500, Coefficient_N)+md.Energies(Coefficient_N)+md.AnalitycalEnergies(Convert.ToDouble(Coefficient_N) )));;
            }
            var pointsAn = new LineSeries { 
                Title="Analitycal"
            };
            for (int i = -500; i < 500; i++)
            {
                pointsAn.Points.Add(new DataPoint(i,md.AnalitycalEigenValue(i,Coefficient_N)+md.AnalitycalEnergies(Convert.ToDouble(Coefficient_N))));
            }
            PlotModel.Series.Add(new FunctionSeries(x => md.CalcY(x+500), -500, 500, 0.5));
            PlotModel.Series.Add(points);
            PlotModel.Series.Add(pointsAn);
            PlotModel.InvalidatePlot(true);


        }

        private void GenerateFile()
        {
            var md = new MartinDin(Coefficient_K);
            md.Calculate();
            string path = @"C:\Users\Vosming\source\repos\Zad3\Zad3.txt";
            string[] lines = new string[12+Coefficient_N];
            lines[0] = "Eigen Energies";
            for(int i=1;i<11;i++)
            {
                lines[i] = $"{md.Energy(i-1)}";
            }
            lines[11] = "EigenFunctions";
            for(int i = 0; i<Coefficient_N;i++)
            {
                lines[12 + i] = $" {i}  {md.CalcY(i)}";
            }
            System.IO.File.WriteAllLines(path, lines);




        }
        private void GenerateChartNew()
        {
            PlotModel.Series.Clear();
            double Width = 3.0E-9; // szerokosc badanego przedziału [m]
            int nz = 201; // liczba punktow siatki
            var md = new Calculus(Coefficient_N,Coefficient_K);
            double M = 9.109E-31; //electron mass [kg]
            double Hbar = 1.055E-34; //Dirac constant [J*s]
            double dZ = Width / Convert.ToDouble(nz-1); // net range [m]
            double Vc = 1.0; //unit exchange [J]
            double Lc = 1.0; //unit exchange [m]
            double precision = 1E-14; //enegry calculation precision
            double beta = 2.0 * M * Vc * Math.Pow(Lc * dZ / Hbar, 2); //from out of scale to SI

            double[] v = new double[nz];
            double[] e = new double[nz];
            double[] d = new double[nz];

            for (int i = 0; i < nz; i++)
            {
                double Z = (-0.5 * Width + Convert.ToDouble(i) * dZ);
                v[i]=0.5 * Coefficient_K * Z * Z;
                d[i] = 2 + beta * v[i];
                e[i] = -1.0;
            }
            double[] mdVec = new double[nz];
            mdVec = md.CalculatedEigenVector(e,d,md.CalculateEigenValue(e,d,Coefficient_N,precision),md.PositionOfMinimum(v)+Coefficient_N);
            //for (int i = 0; i < nz; i++)
            //{
            //    mdVec[i] = (mdVec[i]/Math.Pow(10,22))+(MDEnergy*Math.Pow(10,22));
            //}

            double[] anVec = new double[nz];
            anVec = md.GenerateAnalitycalEigenVector(Coefficient_N, Width, nz, Coefficient_K);

            MDEnergy = md.CalculateEigenValue(e, d, Coefficient_N, precision)/beta;
            AnEnergy = Hbar * Math.Sqrt(Coefficient_K / M) * (Convert.ToDouble(Coefficient_N) + 0.5);

            var pointsmd = new LineSeries { 
                Title="Martin Deen"
            };
            for (int i = 0; i < nz; i++)
            {
                pointsmd.Points.Add(new DataPoint(i,(mdVec[i]/500)+MDEnergy));
            }
            
            var pointsan = new LineSeries
            {
                Title = "Analitical"
            };
            for (int i = 0; i < nz; i++)
            {
                pointsan.Points.Add(new DataPoint(i, (anVec[i]/3500)+AnEnergy));
            }
            var potentialpoints = new LineSeries
            {
                Title = "Parabollic Potential"
            };
            for (int i =0; i <nz; i++)
            {
                //potentialpoints.Points.Add(new DataPoint(i,v[i]));
                potentialpoints.Points.Add(new DataPoint(i, v[i]));

            }
            PlotModel.Series.Add(pointsmd);
            //PlotModel.Series.Add(pointsan);
            PlotModel.Series.Add(potentialpoints);
            PlotModel.InvalidatePlot(true);

            string path = @"C:\Users\Vosming\source\repos\Zad3\Values.txt";
            string[] lines = new string[nz+1];
            lines[0] = "MD_Energy;MD_EigenVector;Analitycal_Energy;Analitycal_EigenVector";
            for (int i = 1; i < nz+1; i++)
            {
                if (i == 1) lines[i] = $"{MDEnergy};{mdVec[i-1]};{AnEnergy};{anVec[i-1]}";
                else lines[i] = $" ;{mdVec[i - 1]}; ;{anVec[i - 1]}";
            }
          
            System.IO.File.WriteAllLines(path, lines);


        }
    }
}
