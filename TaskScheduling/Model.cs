using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using LinqToExcel;

namespace TaskScheduling
{

    public class Model
    {
        public int NumberOfFamilies { get; set; } //F
        public int NumberOfProducts { get; set; } //N

        public int Kmin { get; set; }
        public int Kmax { get; set; }

        public int[] DelayImportanceFactor { get; set; } //WTj
        public double[] DelayOfJobs { get; set; } //Tj
        public double[] DueTimeOfJobs { get; set; } //dj
        public int[,] JobBelongsFamily { get; set; } // h[j,f]==1 when job j belongs to family f
        public double[] ProcessTimeOfJobsInStep1 { get; set; } // pjs process time of job j in step 1
        public double[] ProcessTimeOfJobsInStep2 { get; set; } // pjs process time of job j in step 2
        public int[] SizeOfJobs { get; set; } //Sj
        public int MaxNumberOfBatches { get; set; } // B = Ceiling(N/Kmin)
        public int Ww { get; set; } //Ww = a constant cost
        public int[] IgnoreImportanceFactor { get; set; } //sigma j 
        public int[] IgnoranceBinary { get; set; } //Pi j (equals 1 if Rj>0 else equals 0)
        public int NumberOfNonEmptyBatches { get; set; } // B = Ceiling(Kmin/N)
        public int NumberOfMachinesInStep1 { get; set; }
        public int NumberOfMachinesInStep2 { get; set; }

        public Model()
        {

        }
        public Model(Model model)
        {
            this.NumberOfProducts = model.NumberOfProducts;
            this.NumberOfFamilies = model.NumberOfFamilies;
            this.Ww = model.Ww;
            this.Kmin = model.Kmin;
            this.Kmax = model.Kmax;
            this.MaxNumberOfBatches = model.MaxNumberOfBatches;
            this.NumberOfMachinesInStep1 = model.NumberOfMachinesInStep1;
            this.NumberOfMachinesInStep2 = model.NumberOfMachinesInStep2;
            this.DelayImportanceFactor = model.DelayImportanceFactor;
            this.ProcessTimeOfJobsInStep1 = model.ProcessTimeOfJobsInStep1;
            this.ProcessTimeOfJobsInStep2 = model.ProcessTimeOfJobsInStep2;
            this.SizeOfJobs = model.SizeOfJobs;
            this.DueTimeOfJobs = model.DueTimeOfJobs;
            this.JobBelongsFamily = model.JobBelongsFamily;
            this.IgnoreImportanceFactor = model.IgnoreImportanceFactor;
            this.IgnoranceBinary = model.IgnoranceBinary;
            this.DelayOfJobs = model.DelayOfJobs;
        }

        /// <summary>
        /// Model from File
        /// </summary>
        /// <param name="pathToExcelFile"> the input file directory
        ///  </param>
        public void CreateModel(string pathToExcelFile)
        {

            var excelFile = new ExcelQueryFactory(pathToExcelFile);
            var sheetRows = excelFile.Worksheet("init");

            int ww = 0;
            int kMin = 0;
            int kMax = 0;
            int numberOfProducts = 0;
            int numberOfFamilies = 0;
            int numberOfMachines1 = 0;
            int numberOfMachines2 = 0;

            foreach (var row in sheetRows)
            {
                ww = int.Parse(row["ww"]);
                kMin = int.Parse(row["kmin"]);
                kMax = int.Parse(row["kmax"]);
                numberOfProducts = int.Parse(row["NoP"]);
                numberOfFamilies = int.Parse(row["NoF"]);
                numberOfMachines1 = int.Parse(row["NoM1"]);
                numberOfMachines2 = int.Parse(row["NoM2"]);
            }

            int B = numberOfProducts / kMin + 1;
            int[] wTj = new int[numberOfProducts];
            double[] p1j = new double[numberOfProducts];
            double[] p2j = new double[numberOfProducts];
            int[] sj = new int[numberOfProducts];
            double[] dj = new double[numberOfProducts];
            double[] Tj = new double[numberOfProducts];
            int[,] hjf = new int[numberOfFamilies, numberOfProducts];
            int[] sigmaj = new int[numberOfProducts];
            int[] PiJ = new int[numberOfProducts];

            #region Values Assignments

            sheetRows = excelFile.Worksheet("arrays");

            List<int> wT = new List<int>();
            List<double> p1 = new List<double>();
            List<double> p2 = new List<double>();
            List<int> ss = new List<int>();
            List<double> dd = new List<double>();
            List<int> sigj = new List<int>();


            foreach (var row in sheetRows)
            {
                wT.Add(int.Parse(row["wTj"]));
                p1.Add(double.Parse(row["p1j"]));
                p2.Add(double.Parse(row["p2j"]));
                ss.Add(int.Parse(row["Sj"]));
                dd.Add(Convert.ToDouble(row["dj"]));
                sigj.Add(int.Parse(row["SigmaJ"]));

            }

            wTj = wT.ToArray();
            p1j = p1.ToArray();
            p2j = p2.ToArray();
            sj = ss.ToArray();
            dj = dd.ToArray();
            sigmaj = sigj.ToArray();

            sheetRows = excelFile.Worksheet("hj");

            int i = 0;
            foreach (var row in sheetRows)
            {
                for (int j = 0; j < row.Count; j++)
                {
                    hjf[j, i] = int.Parse(row[j]);
                }
                i++;
            }

            #endregion

            this.NumberOfProducts = numberOfProducts;
            this.NumberOfFamilies = numberOfFamilies;
            this.Ww = ww;
            this.Kmin = kMin;
            this.Kmax = kMax;
            this.MaxNumberOfBatches = B;
            this.NumberOfMachinesInStep1 = numberOfMachines1;
            this.NumberOfMachinesInStep2 = numberOfMachines2;
            this.DelayImportanceFactor = wTj;
            this.ProcessTimeOfJobsInStep1 = p1j;
            this.ProcessTimeOfJobsInStep2 = p2j;
            this.SizeOfJobs = sj;
            this.DueTimeOfJobs = dj;
            this.JobBelongsFamily = hjf;
            this.IgnoreImportanceFactor = sigmaj;
            this.IgnoranceBinary = PiJ;
            this.DelayOfJobs = Tj;
        }

        public double CostFunction()
        {
            double z = 0;

            int[] wT = this.DelayImportanceFactor;
            int[] sigma = this.IgnoreImportanceFactor;
            int[] Pi = this.IgnoranceBinary;
            int Vb = this.NumberOfNonEmptyBatches;
            int Ww = this.Ww;
            double[] T = this.DelayOfJobs;

            for (int j = 0; j < wT.Length; j++)
                z += ((wT[j] * T[j]) + (sigma[j] * Pi[j]));

            z += Ww * Vb;

            return z;

        }

    }

}
