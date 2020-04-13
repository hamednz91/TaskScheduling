using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.Dynamic;
using System.IO;
using System.Linq;
using System.Net;
using System.Reflection.Emit;
using System.Runtime.Remoting.Messaging;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using LinqToExcel;
using Remotion.Data.Linq.Clauses;

namespace TaskScheduling
{
    class Program
    {
        struct Sol
        {

            public List<Batch> BatchesAllocatedToMachines;

            public double[] TimeofMachinesStep1;

            public double[] TimeofMachinesStep2;

            public double[] Tj;
        }

        struct Model
        {
            public int NumberOfFamilies; //F
            public int NumberOfProducts; //N

            public int Kmin;
            public int Kmax;

            public int[] DelayImportanceFactor; //WTj
            public double[] DelayOfJobs; //Tj
            public double[] DueTimeOfJobs; //dj
            public int[,] JobBelongsFamily; // h[j,f]==1 when job j belongs to family f
            public double[] ProcessTimeOfJobsInStep1; // pjs process time of job j in step 1
            public double[] ProcessTimeOfJobsInStep2; // pjs process time of job j in step 2
            public int[] SizeOfJobs; //Sj

            public int MaxNumberOfBatches; // B = Ceiling(N/Kmin)

            public int Ww; //Ww = a constant cost

            public int[] IgnoreImportanceFactor; //sigma j 

            public int[] IgnoranceBinary; //Pi j (equals 1 if Rj>0 else equals 0)

            public int NumberOfNonEmptyBatches; // B = Ceiling(Kmin/N)

            public int NumberOfMachinesInStep1;

            public int NumberOfMachinesInStep2;
        }

        struct Batch
        {
            public List<int> BatchCandidateList;

            public List<double> Eta1jCandidates;

            public List<double> Eta2jCandidates;

            public List<double> Eta3jCandidates;

            public List<double> Eta4jCandidates;

            public List<double> TauJBCandidates;

            public List<double> TauJCandidates;

            public List<double> JobCandidateSelectionProbability;

            public List<int> JobsIndice;

            public List<int> SizeOfJobs;

            public List<double> UrgentMetric;

            //public int SumOfJobSizes;

            public int[] machineNumber; //machineNumber[in_step1,in_step2]

            public int[,] Pjbs; //job j proccessing time of batch b in step s (1 and 2 [step1,step2] )

            public double[] Pbs; //batch proccessing time in step 1 and 2 [step1,step2]

            public int Family;

            public int batchIndex;

            public double AverageDueTimeofJobToDelayImportanceFactor; // dj/wj


        }

        struct Ant
        {
            public List<Batch> Tour;

            public double Cost;

            public Sol sol;

            public bool[] SelectedJobs;

            public bool[] SelectedJobsForEmptyBatches;
        }

        /// <summary>
        /// Model from File
        /// </summary>
        /// <param name="a"> an arbitrary parameter for using this method in the body of code. 
        ///     for using 5 jobs model an arbitrary number is required to invoke this method. 
        /// ex. CreateModel(1)
        ///  </param>
        static Model CreateModel(string pathToExcelFile)
        {


            var excelFile = new ExcelQueryFactory(pathToExcelFile);
            var sheetRows = excelFile.Worksheet("init").Select(a => a);


            Model model = new Model();

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

            sheetRows = excelFile.Worksheet("arrays").Select(a => a);

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

            sheetRows = excelFile.Worksheet("hj").Select(a => a);

            int f = 0;
            foreach (var row in sheetRows)
            {
                for (int j = 0; j < row.Count; j++)
                    hjf[f, j] = int.Parse(row[j]);
                f++;
            }

            #endregion

            model.NumberOfProducts = numberOfProducts;

            model.NumberOfFamilies = numberOfFamilies;

            model.Ww = ww;

            model.Kmin = kMin;

            model.Kmax = kMax;

            model.MaxNumberOfBatches = B;

            model.NumberOfMachinesInStep1 = numberOfMachines1;

            model.NumberOfMachinesInStep2 = numberOfMachines2;

            model.DelayImportanceFactor = wTj;

            model.ProcessTimeOfJobsInStep1 = p1j;

            model.ProcessTimeOfJobsInStep2 = p2j;

            model.SizeOfJobs = sj;

            model.DueTimeOfJobs = dj;

            model.JobBelongsFamily = hjf;

            model.IgnoreImportanceFactor = sigmaj;

            model.IgnoranceBinary = PiJ;

            model.DelayOfJobs = Tj;


            return model;
        }

        /// <summary>
        /// Model for 50 jobs
        /// </summary>
        static bool CheckExistInBatch(Batch batch, int jobIndex)
        {
            return batch.JobsIndice.Any(t => t == jobIndex);
        }

        static bool CheckInSameFamily(int[,] belong2family, int job1, int job2, int family)
        {
            int[,] hjf = belong2family;

            if (hjf[family, job1] == 1 && hjf[family, job2] == 1)
            {
                return true;
            }
            return false;
        }

        static Sol Algorithm1(int option, List<Batch> noneEmptybatches, double[] t1, double[] t2, double[] Tj, double[] d,
              double t_Now)
        {
            Sol sol = new Sol();

            switch (option)
            {
                #region Case 1

                case 1:
                    List<Batch> myBatches = noneEmptybatches;

                    for (int i = 0; i < myBatches.Count; i++)
                    {

                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (myBatches[i].Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += myBatches[i].Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + myBatches[i].Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += myBatches[i].Pbs[1];
                        }

                        myBatches[i].machineNumber[0] = minIndexT1;

                        myBatches[i].machineNumber[1] = minIndexT2;

                        foreach (int j in myBatches[i].JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);
                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    sol.BatchesAllocatedToMachines = noneEmptybatches;


                    break;


                #endregion

                #region Case 2

                case 2:
                    List<Batch> nonEmptyBatchesSortedByP1j = noneEmptybatches;

                    nonEmptyBatchesSortedByP1j = nonEmptyBatchesSortedByP1j.OrderBy(a => a.Pbs[0]).ToList();

                    foreach (var batch in nonEmptyBatchesSortedByP1j)
                    {
                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (batch.Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += batch.Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + batch.Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += batch.Pbs[1];
                        }

                        batch.machineNumber[0] = minIndexT1;

                        batch.machineNumber[1] = minIndexT2;

                        foreach (int j in batch.JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);
                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    sol.BatchesAllocatedToMachines = nonEmptyBatchesSortedByP1j;

                    break;


                #endregion

                #region Case 3

                case 3:
                    List<Batch> nonEmptyBatchesSortedByP1jPlusP2j = noneEmptybatches;

                    nonEmptyBatchesSortedByP1jPlusP2j =
                        nonEmptyBatchesSortedByP1jPlusP2j.OrderBy(a => a.Pbs[0] + a.Pbs[1]).ToList();

                    foreach (var batch in nonEmptyBatchesSortedByP1jPlusP2j)
                    {
                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (batch.Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += batch.Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + batch.Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += batch.Pbs[1];
                        }

                        batch.machineNumber[0] = minIndexT1;

                        batch.machineNumber[1] = minIndexT2;

                        foreach (int j in batch.JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);
                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    sol.BatchesAllocatedToMachines = nonEmptyBatchesSortedByP1jPlusP2j;
                    break;


                #endregion

                #region Case 4

                case 4:

                    List<Batch> nonEmptyBatchesSortedByMeanUrgent = noneEmptybatches;

                    for (int b = 0; b < nonEmptyBatchesSortedByMeanUrgent.Count; b++)
                    {
                        for (int j = 0; j < nonEmptyBatchesSortedByMeanUrgent[b].UrgentMetric.Count; j++)
                        {
                            nonEmptyBatchesSortedByMeanUrgent[b].UrgentMetric[j] =
                                (double)(d[nonEmptyBatchesSortedByMeanUrgent[b].JobsIndice[j]] - t_Now) /
                                (double)(nonEmptyBatchesSortedByMeanUrgent[b].Pbs[0] +
                                          nonEmptyBatchesSortedByMeanUrgent[b].Pbs[1]);
                        }
                    }

                    nonEmptyBatchesSortedByMeanUrgent =
                        nonEmptyBatchesSortedByMeanUrgent.OrderBy(a => a.UrgentMetric.Average()).ToList();

                    foreach (var batch in nonEmptyBatchesSortedByMeanUrgent)
                    {
                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (batch.Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += batch.Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + batch.Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += batch.Pbs[1];
                        }

                        t_Now = t2[minIndexT2];

                        batch.machineNumber[0] = minIndexT1;

                        batch.machineNumber[1] = minIndexT2;


                        foreach (int j in batch.JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);
                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    sol.BatchesAllocatedToMachines = nonEmptyBatchesSortedByMeanUrgent;

                    break;

                #endregion

                #region Case 5

                case 5:

                    List<Batch> allocatedBatchesSortedByBatchIndex =
                        noneEmptybatches.OrderBy(a => a.batchIndex).ToList();

                    foreach (var batch in allocatedBatchesSortedByBatchIndex)
                    {
                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (batch.Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += batch.Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + batch.Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += batch.Pbs[1];
                        }

                        batch.machineNumber[0] = minIndexT1;

                        batch.machineNumber[1] = minIndexT2;

                        foreach (int j in batch.JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);

                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    sol.BatchesAllocatedToMachines = allocatedBatchesSortedByBatchIndex;

                    break;


                    #endregion

            }

            sol.Tj = Tj;

            return sol;
        }

        static double CostFunction(Model model)
        {
            double z = 0;

            int[] wT = model.DelayImportanceFactor;

            int[] sigma = model.IgnoreImportanceFactor;

            int[] Pi = model.IgnoranceBinary;

            int Vb = model.NumberOfNonEmptyBatches;

            int Ww = model.Ww;

            double[] T = model.DelayOfJobs;

            for (int j = 0; j < wT.Length; j++)
                z += ((wT[j] * T[j]) + (sigma[j] * Pi[j]));

            z += Ww * Vb;

            return z;

        }

        static int RouletteWheelSelection(double[] pr)
        {
            double rand = new Random().NextDouble();

            double sum = 0;

            var output = pr.Select(w => sum += w).ToArray();

            for (int i = 0; i < output.Length; i++)
            {
                if (rand < output[i])
                {
                    return i;
                }
            }

            Console.WriteLine("Error In Roulette Wheel Selection !!!");

            return 0;
        }

        static void Run_ACO(string inputFilePath)
        {
            System.Diagnostics.Stopwatch sw = new System.Diagnostics.Stopwatch();

            sw.Start();

            #region Output File Preparation 

            try
            {
                File.Delete("Output.csv");
            }
            catch (FileNotFoundException)
            {
                throw;
            }

            StreamWriter streamWriter = new StreamWriter("Output.csv", true);

            #endregion

            #region Problem Definition

            Model model = CreateModel(inputFilePath);

            int[,] hjf = model.JobBelongsFamily;

            int kMax = model.Kmax;

            int kMin = model.Kmin;

            int[] Sj = model.SizeOfJobs;

            int N = model.NumberOfProducts;

            int B = model.MaxNumberOfBatches;

            int[] wTj = model.DelayImportanceFactor;

            int[] sigmaJ = model.IgnoreImportanceFactor;

            int[] PiJ = model.IgnoranceBinary;

            double[] d = model.DueTimeOfJobs;

            double[] Tj = model.DelayOfJobs;

            double[] p1 = model.ProcessTimeOfJobsInStep1;

            double[] p2 = model.ProcessTimeOfJobsInStep2;

            Random r = new Random();

            int rand = r.Next(N);

            bool[] selectedJobs = new bool[N];

            bool[] selectedJobsForEmptyBatches = new bool[N];

            double[] t1 = new double[model.NumberOfMachinesInStep1];

            double[] t2 = new double[model.NumberOfMachinesInStep2];

            int A = N;

            #endregion

            #region ACO Prarameters 

            int maxIteration = 1000;

            int numberOfAnts = 70;

            int Q = 10;

            double alpha = 1;

            double beta1 = 1;
            double beta2 = 1;
            double beta3 = 1;
            double beta4 = 1;

            double rho = 0.05;

            #endregion

            #region ACO Initialization

            double[] eta1J = new double[N];

            double[] eta2J = new double[N];

            double[] eta3J = new double[N];

            double[] eta4J = new double[N];

            double[,] phiJX = new double[N, N];

            int[,] mJX = new int[N, N];

            double[] R = new double[N];

            double[] tauJB = new double[N];

            double[] tauJ = new double[N];

            for (int i = 0; i < N; i++)
            {
                eta1J[i] = wTj[i];

                eta2J[i] = (double)1 / (d[i] - (double)(p1[i] + p2[i]));

                eta3J[i] = 1;

                eta4J[i] = sigmaJ[i];

                tauJB[i] = 1;

                tauJ[i] = 1;

                for (int j = 0; j < N; j++)
                {
                    phiJX[i, j] = 0.1;
                }
            }

            Ant[] ant = new Ant[numberOfAnts];

            Ant bestAnt = new Ant();

            bestAnt.Tour = new List<Batch>();

            bestAnt.Cost = double.MaxValue;

            bestAnt.SelectedJobs = new bool[N];

            bestAnt.SelectedJobsForEmptyBatches = new bool[N];

            Sol sol = new Sol();

            sol.BatchesAllocatedToMachines = new List<Batch>();

            sol.TimeofMachinesStep1 = new double[t1.Length];

            sol.TimeofMachinesStep2 = new double[t2.Length];

            #endregion

            #region ACO Main Loop

            for (int it = 0; it < maxIteration; it++)
            {
                int t_now = 0;

                R = new double[N];

                // Ants Movement
                for (int k = 0; k < numberOfAnts; k++)
                {
                    sol = new Sol();

                    ant[k] = new Ant();

                    ant[k].Tour = new List<Batch>();

                    ant[k].sol = new Sol();

                    ant[k].SelectedJobs = new bool[N];

                    ant[k].SelectedJobsForEmptyBatches = new bool[N];

                    t1 = new double[model.NumberOfMachinesInStep1];

                    t2 = new double[model.NumberOfMachinesInStep2];

                    selectedJobs = new bool[N];

                    selectedJobsForEmptyBatches = new bool[N];

                    #region Batches Initialization

                    Batch[] batches = new Batch[B];

                    Batch[] virtualBatches = new Batch[model.NumberOfFamilies];

                    for (int i = 0; i < batches.Length; i++)
                    {
                        batches[i] = new Batch();

                        batches[i].BatchCandidateList = new List<int>();

                        batches[i].JobCandidateSelectionProbability = new List<double>();

                        batches[i].Eta1jCandidates = new List<double>();

                        batches[i].Eta2jCandidates = new List<double>();

                        batches[i].Eta3jCandidates = new List<double>();

                        batches[i].Eta4jCandidates = new List<double>();

                        batches[i].TauJBCandidates = new List<double>();

                        batches[i].TauJCandidates = new List<double>();

                        batches[i].JobsIndice = new List<int>();

                        batches[i].SizeOfJobs = new List<int>();

                        batches[i].UrgentMetric = new List<double>();

                        batches[i].Pjbs = new int[N, 2]; //jobs of batch processing time in 2 steps

                        batches[i].Pbs = new double[2]; //batch processing time in step 1 & 2

                        //batches[i].SumOfJobSizes = 0;

                        batches[i].Family = -1;

                        batches[i].machineNumber = new int[] { -1, -1 };

                        batches[i].batchIndex = -1;

                        batches[i].AverageDueTimeofJobToDelayImportanceFactor = -1;
                    }

                    for (int i = 0; i < virtualBatches.Length; i++)
                    {
                        virtualBatches[i] = new Batch();

                        virtualBatches[i].BatchCandidateList = new List<int>();

                        virtualBatches[i].JobCandidateSelectionProbability = new List<double>();

                        virtualBatches[i].Eta1jCandidates = new List<double>();

                        virtualBatches[i].Eta2jCandidates = new List<double>();

                        virtualBatches[i].Eta3jCandidates = new List<double>();

                        virtualBatches[i].Eta4jCandidates = new List<double>();

                        virtualBatches[i].TauJBCandidates = new List<double>();

                        virtualBatches[i].TauJCandidates = new List<double>();

                        virtualBatches[i].JobsIndice = new List<int>();

                        virtualBatches[i].SizeOfJobs = new List<int>();

                        virtualBatches[i].UrgentMetric = new List<double>();

                        virtualBatches[i].Pjbs = new int[N, 2]; //jobs of batch processing time in 2 steps

                        virtualBatches[i].Pbs = new double[2]; //batch processing time in step 1 & 2

                        //virtualBatches[i].SumOfJobSizes = 0;

                        virtualBatches[i].Family = -1;

                        virtualBatches[i].machineNumber = new int[] { -1, -1 };

                        virtualBatches[i].batchIndex = -1;

                        virtualBatches[i].AverageDueTimeofJobToDelayImportanceFactor = -1;

                    }

                    List<Batch> nonEmptyBatches = new List<Batch>();



                    #endregion

                    #region Construct Batches From Batch Candidate Lists

                    for (int b = 0; b < batches.Length; b++)
                    {
                        #region Random Job Selection

                        int numberOfSelectedJobs = selectedJobs.Count(item => item);

                        int numberOfSelectedJobsForEmptyBatches = selectedJobsForEmptyBatches.Count(item => item);

                        if (numberOfSelectedJobs + numberOfSelectedJobsForEmptyBatches >= N)
                            break;

                        do
                        {
                            rand = r.Next(N);
                        } while (selectedJobsForEmptyBatches[rand] || selectedJobs[rand]);


                        batches[b].JobsIndice.Add(rand);

                        batches[b].SizeOfJobs.Add(Sj[rand]);

                        // batches[b].SumOfJobSizes = Sj[rand];

                        selectedJobs[rand] = true;

                        double maxP1j = p1[rand];

                        double maxP2j = p2[rand];

                        for (int f = 0; f < hjf.GetLength(0); f++)
                            if (hjf[f, rand] == 1)
                            {
                                batches[b].Family = f;
                                break;
                            }


                        #endregion

                        #region Create Jobs Candidate List for Random Job

                        for (int i = 0; i < N; i++)
                        {
                            if (!CheckExistInBatch(batches[b], i) &&
                                CheckInSameFamily(hjf, rand, i,
                                    batches[b].Family) && !selectedJobs[i])
                            {

                                batches[b].BatchCandidateList.Add(i);

                                batches[b].Eta1jCandidates.Add(eta1J[i]);

                                batches[b].Eta2jCandidates.Add(eta2J[i]);

                                batches[b].Eta3jCandidates.Add(eta3J[i]);

                                batches[b].Eta4jCandidates.Add(eta4J[i]);

                                batches[b].TauJBCandidates.Add(tauJB[i]);

                                batches[b].TauJCandidates.Add(tauJ[i]);

                                batches[b].JobCandidateSelectionProbability.Add(1);

                            }

                        }

                        for (int j = 0; j < batches[b].BatchCandidateList.Count; j++)
                        {
                            eta3J[batches[b].BatchCandidateList[j]] = (double)1 /
                                                                      (double)
                                                                      ((kMax - batches[b].SizeOfJobs.Sum()) -
                                                                       Sj[batches[b].BatchCandidateList[j]] + 1);

                            double sumOfPhiJXs = 0;

                            for (int l = 0; l < batches[b].JobsIndice.Count; l++)
                            {
                                // if (l != j)
                                //{
                                sumOfPhiJXs += phiJX[batches[b].BatchCandidateList[j], batches[b].JobsIndice[l]];
                                //  }
                            }

                            tauJB[batches[b].BatchCandidateList[j]] = sumOfPhiJXs / batches[b].JobsIndice.Count;

                        }


                        #endregion

                        #region Pick From Candidate List

                        int counter = 0;

                        int batchCandidateListCount = batches[b].BatchCandidateList.Count;

                        //تا زمانی که دسته جا داره از لیست کاندید بردار و به دسته اضافه کن
                        while (batches[b].SizeOfJobs.Sum() <= kMax && counter < batchCandidateListCount)
                        {
                            double sumOfEta1 = 0;

                            double sumOfEta2 = 0;

                            double sumOfEta3 = 0;

                            double sumOfEta4 = 0;

                            double sumOfTauJB = 0;

                            double sumOfTauJ = 0;

                            double sumOfSum = 0;

                            counter++;

                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                            {

                                batches[b].Eta1jCandidates[l] = eta1J[batches[b].BatchCandidateList[l]];

                                batches[b].Eta2jCandidates[l] = eta2J[batches[b].BatchCandidateList[l]];

                                batches[b].Eta3jCandidates[l] = eta3J[batches[b].BatchCandidateList[l]];

                                batches[b].Eta4jCandidates[l] = eta4J[batches[b].BatchCandidateList[l]];

                                batches[b].TauJBCandidates[l] = tauJB[batches[b].BatchCandidateList[l]];

                                batches[b].TauJCandidates[l] = tauJ[batches[b].BatchCandidateList[l]];

                                sumOfEta1 += batches[b].Eta1jCandidates[l];

                                sumOfEta2 += batches[b].Eta2jCandidates[l];

                                sumOfEta3 += batches[b].Eta3jCandidates[l];

                                sumOfEta4 += batches[b].Eta4jCandidates[l];

                                sumOfTauJB += batches[b].TauJBCandidates[l];

                                sumOfTauJ += batches[b].TauJCandidates[l];
                            }

                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                            {
                                batches[b].Eta1jCandidates[l] = batches[b].Eta1jCandidates[l] / sumOfEta1;

                                batches[b].Eta2jCandidates[l] = batches[b].Eta2jCandidates[l] / sumOfEta2;

                                batches[b].Eta3jCandidates[l] = batches[b].Eta3jCandidates[l] / sumOfEta3;

                                batches[b].Eta4jCandidates[l] = batches[b].Eta4jCandidates[l] / sumOfEta4;

                                batches[b].TauJBCandidates[l] = batches[b].TauJBCandidates[l] / sumOfTauJB;

                                batches[b].TauJCandidates[l] = batches[b].TauJCandidates[l] / sumOfTauJ;
                            }
                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                            {
                                batches[b].JobCandidateSelectionProbability[l] =
                                    (Math.Pow((batches[b].TauJBCandidates[l] + batches[b].TauJCandidates[l]), alpha) *
                                     Math.Pow(batches[b].Eta1jCandidates[l], beta1) *
                                     Math.Pow(batches[b].Eta2jCandidates[l], beta2) *
                                     Math.Pow(batches[b].Eta3jCandidates[l], beta3) *
                                     Math.Pow(batches[b].Eta4jCandidates[l], beta4)
                                    ) /
                                    (Math.Pow((sumOfTauJB + sumOfTauJB), alpha) *
                                     Math.Pow(sumOfEta1, beta1) *
                                     Math.Pow(sumOfEta2, beta2) *
                                     Math.Pow(sumOfEta3, beta3) *
                                     Math.Pow(sumOfEta4, beta4));
                                sumOfSum += batches[b].JobCandidateSelectionProbability[l];
                            }
                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                                batches[b].JobCandidateSelectionProbability[l] =
                                    batches[b].JobCandidateSelectionProbability[l] /
                                    sumOfSum;


                            int ind = RouletteWheelSelection(batches[b].JobCandidateSelectionProbability.ToArray());

                            if (Sj[batches[b].BatchCandidateList[ind]] + batches[b].SizeOfJobs.Sum() > kMax)
                                continue;

                            if (!selectedJobs[batches[b].BatchCandidateList[ind]])
                            {
                                batches[b].JobsIndice.Add(batches[b].BatchCandidateList[ind]);

                                batches[b].SizeOfJobs.Add(Sj[batches[b].BatchCandidateList[ind]]);

                                //batches[b].SizeOfJobs.Sum() += Sj[batches[b].BatchCandidateList[ind]];

                                //for (int l = 0; l < batches[b].Eta3jCandidates.Count; l++)
                                //    eta3J[batches[b].BatchCandidateList[l]] = (double)1 /
                                //                                              (double)
                                //                                              ((kMax - batches[b].SizeOfJobs.Sum()) + 1);

                                selectedJobs[batches[b].BatchCandidateList[ind]] = true;

                                if (p1[batches[b].BatchCandidateList[ind]] > maxP1j)
                                {
                                    maxP1j = p1[batches[b].BatchCandidateList[ind]];
                                }
                                if (p2[batches[b].BatchCandidateList[ind]] > maxP2j)
                                {
                                    maxP2j = p2[batches[b].BatchCandidateList[ind]];
                                }
                                batches[b].BatchCandidateList.RemoveAt(ind);

                                batches[b].JobCandidateSelectionProbability.RemoveAt(ind);

                                batches[b].Eta1jCandidates.RemoveAt(ind);

                                batches[b].Eta2jCandidates.RemoveAt(ind);

                                batches[b].Eta3jCandidates.RemoveAt(ind);

                                batches[b].Eta4jCandidates.RemoveAt(ind);

                                batches[b].TauJBCandidates.RemoveAt(ind);

                                batches[b].TauJCandidates.RemoveAt(ind);

                            }
                        }

                        #endregion

                        #region Remove Batches with Members Lower than Kmin

                        if (batches[b].JobsIndice.Count < kMin)
                        {
                            selectedJobsForEmptyBatches[batches[b].JobsIndice[0]] = true;

                            foreach (int t in batches[b].JobsIndice)
                                selectedJobs[t] = false;

                            batches[b].BatchCandidateList = new List<int>();

                            batches[b].JobCandidateSelectionProbability = new List<double>();

                            batches[b].Eta1jCandidates = new List<double>();

                            batches[b].Eta2jCandidates = new List<double>();

                            batches[b].Eta3jCandidates = new List<double>();

                            batches[b].Eta4jCandidates = new List<double>();

                            batches[b].TauJBCandidates = new List<double>();

                            batches[b].TauJCandidates = new List<double>();

                            batches[b].JobsIndice = new List<int>();

                            batches[b].SizeOfJobs = new List<int>();

                            batches[b].UrgentMetric = new List<double>();

                            batches[b].Pjbs = new int[N, 2];

                            batches[b].Pbs = new double[2];

                            // batches[b].SumOfJobSizes = 0;

                            batches[b].Family = -1;

                            batches[b].batchIndex = -1;

                            batches[b].machineNumber = new int[] { -1, -1 };

                        }

                        #endregion

                        #region Add Pbs

                        batches[b].Pbs[0] = maxP1j;

                        batches[b].Pbs[1] = maxP2j;

                        #endregion


                    }

                    #endregion

                    //for (int j = 0; j < tauJ.Length; j++)
                    //    tauJ[j] = (double)1 / (double)(R[j] + 1);

                    #region Fill Non-EmptyBatches

                    int count = 0;



                    for (int b = 0; b < batches.Length; b++)
                        if (batches[b].JobsIndice.Count > 0)
                        {
                            double average = 0;

                            for (int j = 0; j < batches[b].JobsIndice.Count; j++)
                            {
                                batches[b].UrgentMetric.Add(1);

                                average += d[batches[b].JobsIndice[j]] / wTj[batches[b].JobsIndice[j]];
                            }
                            average = average / (double)batches[b].JobsIndice.Count;

                            batches[b].AverageDueTimeofJobToDelayImportanceFactor = average;

                            batches[b].batchIndex = count++;


                            nonEmptyBatches.Add(batches[b]);
                        }


                    #endregion

                    for (int j = 0; j < selectedJobs.Length; j++)
                        PiJ[j] = !selectedJobs[j] ? 1 : 0;

                    model.IgnoranceBinary = PiJ;

                    model.NumberOfNonEmptyBatches = nonEmptyBatches.Count;

                    sol = Algorithm1(4, nonEmptyBatches, t1, t2, Tj, d, t_now);

                    model.DelayOfJobs = sol.Tj;

                    ant[k].Cost = CostFunction(model);

                    ant[k].sol = sol;

                    ant[k].SelectedJobs = selectedJobs;

                    ant[k].SelectedJobsForEmptyBatches = selectedJobsForEmptyBatches;

                    if (ant[k].Cost < bestAnt.Cost)
                    {
                        bestAnt = ant[k];
                    }

                    nonEmptyBatches = sol.BatchesAllocatedToMachines;

                    #region myNonEmptyBatch Init

                    Batch[] myNonEmptyBatches = new Batch[nonEmptyBatches.Count];

                    for (int x = 0; x < myNonEmptyBatches.Length; x++)
                    {
                        myNonEmptyBatches[x] = nonEmptyBatches[x];

                    }

                    #endregion

                    #region Virtual Batch Init

                    for (int j = 0; j < virtualBatches.Length; j++)
                    {
                        virtualBatches[j] = new Batch();
                        virtualBatches[j].BatchCandidateList = new List<int>();
                        virtualBatches[j].JobCandidateSelectionProbability = new List<double>();
                        virtualBatches[j].Eta1jCandidates = new List<double>();
                        virtualBatches[j].Eta2jCandidates = new List<double>();
                        virtualBatches[j].Eta3jCandidates = new List<double>();
                        virtualBatches[j].Eta4jCandidates = new List<double>();
                        virtualBatches[j].TauJBCandidates = new List<double>();
                        virtualBatches[j].TauJCandidates = new List<double>();
                        virtualBatches[j].JobsIndice = new List<int>();
                        virtualBatches[j].SizeOfJobs = new List<int>();
                        virtualBatches[j].UrgentMetric = new List<double>();
                        virtualBatches[j].Pjbs = new int[N, 2]; //jobs of batch processing time in 2 steps
                        virtualBatches[j].Pbs = new double[2]; //batch processing time in step 1 & 2
                                                               // virtualBatches[j].SumOfJobSizes = 0;
                        virtualBatches[j].Family = -1;
                        virtualBatches[j].machineNumber = new int[] { -1, -1 };
                        virtualBatches[j].batchIndex = -1;
                    }


                    #endregion

                    for (int j = 0; j < selectedJobs.Length; j++)
                    {
                        if (!selectedJobs[j])
                        {
                            //R[j]++;
                            for (int i = 0; i < model.NumberOfFamilies; i++)
                            {
                                if (hjf[i, j] == 1)
                                {
                                    virtualBatches[i].JobsIndice.Add(j);

                                    virtualBatches[i].SizeOfJobs.Add(Sj[j]);

                                    // virtualBatches[i].SumOfJobSizes += Sj[j];

                                    virtualBatches[i].Family = i;

                                    virtualBatches[i].UrgentMetric.Add(1);

                                    if (p1[j] > virtualBatches[i].Pbs[0])
                                    {
                                        virtualBatches[i].Pbs[0] = p1[j];
                                    }
                                    if (p2[j] > virtualBatches[i].Pbs[1])
                                    {
                                        virtualBatches[i].Pbs[1] = p2[j];
                                    }
                                }

                            }
                        }
                    }



                    for (int i = 0; i < A; i++)
                    {
                        int opSelector = RouletteWheelSelection(new[] { .3, .2, .2, .1, .2 });

                        //opSelector = 2;

                        switch (opSelector)
                        {
                            case 0:

                                #region OP0 Drop

                                List<Batch> BatchesGreaterThanKmin = new List<Batch>();

                                foreach (var item in nonEmptyBatches)
                                    if (item.JobsIndice.Count > kMin)
                                        BatchesGreaterThanKmin.Add(item);

                                bool stopOP1flag = BatchesGreaterThanKmin.Count == 0;

                                if (stopOP1flag) break;

                                int selectedBatchIndex = r.Next(BatchesGreaterThanKmin.Count);

                                //selectedBatchIndex = BatchesGreaterThanKmin[selectedBatchIndex].batchIndex;
                                // real index in nonemptyBatches using batch index field

                                int selectedBatchFamily = BatchesGreaterThanKmin[selectedBatchIndex].Family;

                                int selectedBatchLength = BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice.Count;

                                int a = r.Next(1, Math.Max(1, selectedBatchLength - kMin));

                                for (int j = 0; j < a; j++)
                                {

                                    int jobIndex = r.Next(selectedBatchLength);

                                    if (selectedBatchLength > 0)
                                    {
                                        #region Add Selected Jobs to Virtual Batch

                                        virtualBatches[selectedBatchFamily].JobsIndice.Add(
                                            BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]);

                                        virtualBatches[selectedBatchFamily].SizeOfJobs.Add(
                                            Sj[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]]);

                                        //virtualBatches[selectedBatchFamily].SizeOfJobs.Sum() +=
                                        //    Sj[nonEmptyBatches[selectedBatchIndex].JobsIndice[jobIndex]];

                                        virtualBatches[selectedBatchFamily].Family = selectedBatchFamily;

                                        virtualBatches[selectedBatchFamily].UrgentMetric.Add(1);

                                        #endregion

                                        #region Remove Selected Jobs from NonEmptyBatches and Update Batch (size of jobs, Urgent metric, pbs if needed)

                                        bool flag1 = p1[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] >=
                                                     BatchesGreaterThanKmin[selectedBatchIndex].Pbs[0];

                                        bool flag2 = p2[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] >=
                                                     BatchesGreaterThanKmin[selectedBatchIndex].Pbs[1];

                                        selectedJobs[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] = false;

                                        BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice.RemoveAt(jobIndex);

                                        BatchesGreaterThanKmin[selectedBatchIndex].SizeOfJobs.RemoveAt(jobIndex);

                                        BatchesGreaterThanKmin[selectedBatchIndex].UrgentMetric.RemoveAt(jobIndex);
                                        
                                        //-------------also remove from nonEmptyBatches------------------------
                                        //nonEmptyBatches[selectedBatchIndex].JobsIndice.RemoveAt(jobIndex);

                                        //nonEmptyBatches[selectedBatchIndex].SizeOfJobs.RemoveAt(jobIndex);

                                        //nonEmptyBatches[selectedBatchIndex].UrgentMetric.RemoveAt(jobIndex);

                                        selectedBatchLength = nonEmptyBatches[selectedBatchIndex].JobsIndice.Count;

                                        if (selectedBatchLength <= 0) continue;

                                        if (flag1)
                                        {
                                            double maxP1 = p1[nonEmptyBatches[selectedBatchIndex].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatches[selectedBatchIndex].JobsIndice)
                                                if (p1[t] > maxP1)
                                                    maxP1 = p1[t];

                                            nonEmptyBatches[selectedBatchIndex].Pbs[0] = maxP1;
                                        }
                                        if (flag2)
                                        {
                                            double maxP2 = p2[nonEmptyBatches[selectedBatchIndex].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatches[selectedBatchIndex].JobsIndice)
                                                if (p2[t] > maxP2)
                                                    maxP2 = p2[t];

                                            nonEmptyBatches[selectedBatchIndex].Pbs[1] = maxP2;
                                        }

                                        #endregion
                                    }
                                }

                                #endregion

                                break;
                            case 1:

                                #region OP1 Exchange between 2 nonempty

                                bool stopOP2flag = nonEmptyBatches.Count < 2;

                                if (stopOP2flag) break;

                                bool[] selectbatches2 = new bool[nonEmptyBatches.Count];

                                int selectedBatchIndex1 = r.Next(nonEmptyBatches.Count);

                                selectbatches2[selectedBatchIndex1] = true;

                                int selectedFamily = nonEmptyBatches[selectedBatchIndex1].Family;

                                int numberOfSameFamilyBatches2 =
                                    nonEmptyBatches.Count(item => item.Family == selectedFamily);

                                int selectedBatchIndex2 = -1;


                                do
                                {
                                    selectedBatchIndex2 = r.Next(nonEmptyBatches.Count);

                                    if (!selectbatches2[selectedBatchIndex2])
                                        selectbatches2[selectedBatchIndex2] = true;

                                } while ((selectedBatchIndex1 == selectedBatchIndex2 ||
                                          selectedFamily != nonEmptyBatches[selectedBatchIndex2].Family) &&
                                         selectbatches2.Count(item => item) < nonEmptyBatches.Count);

                                if (selectedFamily != nonEmptyBatches[selectedBatchIndex2].Family) break;

                                int selectedBatchLength1 = nonEmptyBatches[selectedBatchIndex1].JobsIndice.Count;
                                int selectedBatchLength2 = nonEmptyBatches[selectedBatchIndex2].JobsIndice.Count;

                                //
                                //can number 2 be chosen in the random selection below?
                                //
                                int n = r.Next(1, Math.Min(selectedBatchLength1, selectedBatchLength2));

                                for (int j = 0; j < n; j++)
                                {
                                    int batchSize1, batchSize2 = 0;

                                    int jobIndex1 = r.Next(selectedBatchLength1);

                                    int jobIndex2 = -1;

                                    int numberOfSelectedBatches2 = 0;

                                    bool[] selectJobIndex2 = new bool[selectedBatchLength2];

                                    do
                                    {
                                        int numberofSelectedJobfromBatch2 = selectJobIndex2.Count(item => item);
                                        if (numberofSelectedJobfromBatch2 >= selectedBatchLength2)
                                        {
                                            do
                                            {

                                                selectedBatchIndex2 = r.Next(nonEmptyBatches.Count);

                                                // in batch entekhab nashode bashe
                                            } while (selectbatches2[selectedBatchIndex2] ||
                                                     selectedBatchIndex1 == selectedBatchIndex2 ||
                                                     selectedFamily != nonEmptyBatches[selectedBatchIndex2].Family);

                                            if (!selectbatches2[selectedBatchIndex2])
                                                selectbatches2[selectedBatchIndex2] = true;

                                            selectedBatchLength2 =
                                                nonEmptyBatches[selectedBatchIndex2].JobsIndice.Count;

                                            selectJobIndex2 = new bool[selectedBatchLength2];
                                        }

                                        do jobIndex2 = r.Next(selectedBatchLength2); while (selectJobIndex2[jobIndex2]);

                                        int a1 = nonEmptyBatches[selectedBatchIndex1].SizeOfJobs.Sum() -
                                                 nonEmptyBatches[selectedBatchIndex1].SizeOfJobs[jobIndex1];

                                        batchSize1 = nonEmptyBatches[selectedBatchIndex2].SizeOfJobs[jobIndex2] +
                                                     a1;

                                        int a2 = nonEmptyBatches[selectedBatchIndex2].SizeOfJobs.Sum() -
                                                 nonEmptyBatches[selectedBatchIndex2].SizeOfJobs[jobIndex2];

                                        batchSize2 = nonEmptyBatches[selectedBatchIndex1].SizeOfJobs[jobIndex1] +
                                                     a2;

                                        if (!selectJobIndex2[jobIndex2])
                                            selectJobIndex2[jobIndex2] = true;

                                        numberOfSelectedBatches2 = selectbatches2.Count(item => item);

                                    } while (numberOfSelectedBatches2 < numberOfSameFamilyBatches2 - 1
                                             && (batchSize1 > kMax || batchSize2 > kMax));

                                    if (numberOfSelectedBatches2 < numberOfSameFamilyBatches2 - 1)
                                    {

                                        int job1 = nonEmptyBatches[selectedBatchIndex1].JobsIndice[jobIndex1];
                                        int job2 = nonEmptyBatches[selectedBatchIndex2].JobsIndice[jobIndex2];

                                        int jobSize1 = nonEmptyBatches[selectedBatchIndex1].SizeOfJobs[jobIndex1];
                                        int jobSize2 = nonEmptyBatches[selectedBatchIndex2].SizeOfJobs[jobIndex2];


                                        nonEmptyBatches[selectedBatchIndex1].JobsIndice.Add(job2);
                                        nonEmptyBatches[selectedBatchIndex2].JobsIndice.Add(job1);

                                        nonEmptyBatches[selectedBatchIndex1].SizeOfJobs.Add(jobSize2);
                                        nonEmptyBatches[selectedBatchIndex2].SizeOfJobs.Add(jobSize1);


                                        nonEmptyBatches[selectedBatchIndex1].JobsIndice.RemoveAt(jobIndex1);
                                        nonEmptyBatches[selectedBatchIndex2].JobsIndice.RemoveAt(jobIndex2);

                                        nonEmptyBatches[selectedBatchIndex1].SizeOfJobs.RemoveAt(jobIndex1);
                                        nonEmptyBatches[selectedBatchIndex2].SizeOfJobs.RemoveAt(jobIndex2);

                                        double maxP1 = p1[nonEmptyBatches[selectedBatchIndex1].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex1].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatches[selectedBatchIndex1].Pbs[0] = maxP1;

                                        double maxP2 = p2[nonEmptyBatches[selectedBatchIndex1].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex1].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatches[selectedBatchIndex1].Pbs[1] = maxP2;


                                        double maxP12 = p1[nonEmptyBatches[selectedBatchIndex2].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex2].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        nonEmptyBatches[selectedBatchIndex2].Pbs[0] = maxP12;

                                        double maxP22 = p2[nonEmptyBatches[selectedBatchIndex2].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex2].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        nonEmptyBatches[selectedBatchIndex2].Pbs[1] = maxP22;
                                    }
                                }

                                #endregion

                                break;
                            case 2:

                                #region OP2 Create new batch

                                bool stopOP3Flag = virtualBatches.All(vb => vb.JobsIndice.Count < kMin);

                                if (stopOP3Flag) break;

                                #region New Batch Instance

                                Batch newBatch = new Batch();

                                #endregion

                                int virtualBatchIndex = -1;

                                List<Batch> VirtualBatchesMoreThanKmin = new List<Batch>();


                                foreach (var batch in virtualBatches)
                                    if (batch.JobsIndice.Count >= kMin)
                                        VirtualBatchesMoreThanKmin.Add(batch);

                                bool[] selectedFamilyVirtualBatch = new bool[VirtualBatchesMoreThanKmin.Count];

                                #region New Batch Init/Reset

                                newBatch.JobsIndice = new List<int>();
                                newBatch.SizeOfJobs = new List<int>();
                                newBatch.UrgentMetric = new List<double>();
                                newBatch.Pbs = new double[2]; //batch processing time in step 1 & 2
                                                              //newBatch.SizeOfJobs.Sum() = 0;
                                newBatch.Family = -1;
                                newBatch.machineNumber = new int[] { -1, -1 };
                                newBatch.batchIndex = -1;

                                #endregion

                                //if (selectedFamilyCounter >= virtualBatches.Length)
                                //    break;

                                do
                                {
                                    virtualBatchIndex = r.Next(VirtualBatchesMoreThanKmin.Count);

                                } while (selectedFamilyVirtualBatch[virtualBatchIndex] &&
                                         selectedFamilyVirtualBatch.Count(item => item) <
                                         VirtualBatchesMoreThanKmin.Count);

                                if (!selectedFamilyVirtualBatch[virtualBatchIndex])
                                    selectedFamilyVirtualBatch[virtualBatchIndex] = true;


                                //if (selectedFamilyVirtualBatch.Count(item => item) >= VirtualBatchesMoreThanKmin.Count) break;


                                int selectedVirBatchLength =
                                    VirtualBatchesMoreThanKmin[virtualBatchIndex].JobsIndice.Count;

                                double avgSjOfSelectedVirBatch =
                                    VirtualBatchesMoreThanKmin[virtualBatchIndex].SizeOfJobs.Average();

                                int nm = r.Next(kMin, Convert.ToInt32(Math.Ceiling((kMax) / avgSjOfSelectedVirBatch)));

                                int jobIndexOP2 = -1;

                                int jobOP2 = -1;

                                bool[] selectedJobFromSelectedVirtualBatchOP3 = new bool[selectedVirBatchLength];

                                for (int j = 0; j < Math.Min(nm, selectedVirBatchLength); j++)
                                {
                                    do
                                    {
                                        jobIndexOP2 = r.Next(selectedVirBatchLength);

                                        jobOP2 = VirtualBatchesMoreThanKmin[virtualBatchIndex].JobsIndice[jobIndexOP2];

                                    } while (selectedJobFromSelectedVirtualBatchOP3[jobIndexOP2] &&
                                             selectedJobFromSelectedVirtualBatchOP3.Count(item => item) <
                                             selectedVirBatchLength);

                                    if (!selectedJobFromSelectedVirtualBatchOP3[jobIndexOP2] &&
                                        Sj[jobOP2] + newBatch.SizeOfJobs.Sum() < kMax)
                                    {
                                        newBatch.JobsIndice.Add(jobOP2);
                                        newBatch.SizeOfJobs.Add(Sj[jobOP2]);
                                        //newBatch.SizeOfJobs.Sum() += Sj[jobOP2];
                                        newBatch.UrgentMetric.Add(1);
                                        newBatch.Family = virtualBatchIndex;

                                        if (!selectedJobFromSelectedVirtualBatchOP3[jobIndexOP2])
                                            selectedJobFromSelectedVirtualBatchOP3[jobIndexOP2] = true;

                                        double maxP1 = p1[jobOP2];

                                        foreach (int t in newBatch.JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        newBatch.Pbs[0] = maxP1;

                                        double maxP2 = p2[jobOP2];

                                        foreach (int t in newBatch.JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        newBatch.Pbs[1] = maxP2;

                                    }
                                }


                                List<Batch> nonEmptyBatchesAverageDjToWj = nonEmptyBatches;

                                // if (selectedFamilyVirtualBatch.Count(item => item) >= VirtualBatchesMoreThanKmin.Count) break;

                                if (newBatch.JobsIndice.Count >= kMin)
                                {
                                    newBatch.batchIndex = nonEmptyBatches.Count;

                                    double average = 0;
                                    foreach (var item in newBatch.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        selectedJobs[item] = true;

                                    }
                                    average = average / (double)newBatch.JobsIndice.Count;

                                    newBatch.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWj.Count; j++)
                                    {
                                        if (average > nonEmptyBatchesAverageDjToWj[j].AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            nonEmptyBatchesAverageDjToWj.Insert(j + 1, newBatch);
                                            isSet = true;
                                            for (int l = j + 2; l < nonEmptyBatchesAverageDjToWj.Count; l++)
                                            {
                                                nonEmptyBatchesAverageDjToWj.ToArray()[l].batchIndex++;
                                            }
                                        }
                                    }

                                    if (!isSet)
                                    {
                                        nonEmptyBatchesAverageDjToWj.Insert(0, newBatch);
                                    }

                                }

                                nonEmptyBatches = nonEmptyBatchesAverageDjToWj;

                                #endregion

                                break;
                            case 3:

                                #region OP3 Exchange between virtual batch and nonempty

                                bool stopOP4Flag = virtualBatches.All(vb => vb.JobsIndice.Count == 0) ||
                                                   nonEmptyBatches.Count == 0;

                                if (stopOP4Flag) break;

                                int op4SelectedBatchIndex = -1;

                                int op4selectedBatchFamily = -1;

                                bool changeSelectedBatchOP4Flag = false;

                                do
                                {
                                    op4SelectedBatchIndex = r.Next(nonEmptyBatches.Count);

                                    op4selectedBatchFamily = nonEmptyBatches[op4SelectedBatchIndex].Family;

                                    changeSelectedBatchOP4Flag =
                                        virtualBatches[op4selectedBatchFamily].JobsIndice.Count == 0;

                                } while (changeSelectedBatchOP4Flag);


                                int op4selectedBatchLength =
                                    nonEmptyBatches[op4SelectedBatchIndex].JobsIndice.Count;


                                int op4n = r.Next(1,
                                    Math.Max(1,
                                        Math.Min(op4selectedBatchLength,
                                            virtualBatches[op4selectedBatchFamily].JobsIndice.Count)));

                                for (int j = 0; j < op4n; j++)
                                {
                                    int batchSize1 = 0;

                                    int jobIndex1 = r.Next(op4selectedBatchLength);

                                    int jobIndex2 = -1;

                                    bool[] selectJobIndex2 =
                                        new bool[virtualBatches[op4selectedBatchFamily].JobsIndice.Count];

                                    int numberofSelectedJobfromBatch4 = 0;

                                    do
                                    {
                                        //check if jobindex2 is already chosen
                                        numberofSelectedJobfromBatch4 = selectJobIndex2.Count(item => item);

                                        jobIndex2 = r.Next(virtualBatches[op4selectedBatchFamily].JobsIndice.Count);

                                        int a1 = nonEmptyBatches[op4SelectedBatchIndex].SizeOfJobs.Sum() -
                                                 nonEmptyBatches[op4SelectedBatchIndex].SizeOfJobs[jobIndex1];

                                        batchSize1 = virtualBatches[op4selectedBatchFamily].SizeOfJobs[jobIndex2] +
                                                     a1;

                                        if (!selectJobIndex2[jobIndex2])
                                            selectJobIndex2[jobIndex2] = true;

                                    } while (batchSize1 > kMax && numberofSelectedJobfromBatch4 <
                                             virtualBatches[op4selectedBatchFamily].JobsIndice.Count);

                                    if (numberofSelectedJobfromBatch4 >=
                                        virtualBatches[op4selectedBatchFamily].JobsIndice.Count)
                                        break;

                                    int job1 = nonEmptyBatches[op4SelectedBatchIndex].JobsIndice[jobIndex1];
                                    int job2 = virtualBatches[op4selectedBatchFamily].JobsIndice[jobIndex2];

                                    int jobSize1 = nonEmptyBatches[op4SelectedBatchIndex].SizeOfJobs[jobIndex1];
                                    int jobSize2 = virtualBatches[op4selectedBatchFamily].SizeOfJobs[jobIndex2];


                                    nonEmptyBatches[op4SelectedBatchIndex].JobsIndice.Add(job2);
                                    virtualBatches[op4selectedBatchFamily].JobsIndice.Add(job1);

                                    nonEmptyBatches[op4SelectedBatchIndex].SizeOfJobs.Add(jobSize2);
                                    virtualBatches[op4selectedBatchFamily].SizeOfJobs.Add(jobSize1);

                                    bool flag1 = p1[job1] >=
                                                 nonEmptyBatches[op4SelectedBatchIndex].Pbs[0];

                                    bool flag2 = p2[job1] >=
                                                 nonEmptyBatches[op4SelectedBatchIndex].Pbs[1];


                                    nonEmptyBatches[op4SelectedBatchIndex].JobsIndice.RemoveAt(jobIndex1);
                                    virtualBatches[op4selectedBatchFamily].JobsIndice.RemoveAt(jobIndex2);

                                    nonEmptyBatches[op4SelectedBatchIndex].SizeOfJobs.RemoveAt(jobIndex1);
                                    virtualBatches[op4selectedBatchFamily].SizeOfJobs.RemoveAt(jobIndex2);

                                    if (flag1)
                                    {
                                        double maxP1 = p1[nonEmptyBatches[op4SelectedBatchIndex].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[op4SelectedBatchIndex].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatches[op4SelectedBatchIndex].Pbs[0] = maxP1;
                                    }

                                    if (flag2)
                                    {
                                        double maxP2 = p2[nonEmptyBatches[op4SelectedBatchIndex].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[op4SelectedBatchIndex].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatches[op4SelectedBatchIndex].Pbs[1] = maxP2;
                                    }

                                    selectedJobs[job1] = false;

                                    selectedJobs[job2] = true;

                                }


                                #endregion

                                break;
                            case 4:

                                #region OP4 Add

                                bool stopOP5Flag = virtualBatches.All(vb => vb.JobsIndice.Count == 0) ||
                                                   nonEmptyBatches.Count == 0;

                                if (stopOP5Flag) break;

                                int selectedVirtualBatchIndex5 = -1;

                                int selectVirBatchCounterOP5 = 0;

                                bool[] selectedFamilyVirtualBatchOP5 = new bool[model.NumberOfFamilies];
                                do
                                {
                                    selectedVirtualBatchIndex5 = r.Next(model.NumberOfFamilies);

                                    if (!selectedFamilyVirtualBatchOP5[selectedVirtualBatchIndex5])
                                    {
                                        selectedFamilyVirtualBatchOP5[selectedVirtualBatchIndex5] = true;
                                        selectVirBatchCounterOP5++;
                                    }

                                } while (virtualBatches[selectedVirtualBatchIndex5].JobsIndice.Count == 0 &&
                                         selectVirBatchCounterOP5 < virtualBatches.Length);

                                int op5SelectedBatchIndex = -1;

                                bool[] selectedBatchOP5 = new bool[nonEmptyBatches.Count];

                                do
                                {
                                    op5SelectedBatchIndex = r.Next(nonEmptyBatches.Count);

                                    if (!selectedBatchOP5[op5SelectedBatchIndex])
                                        selectedBatchOP5[op5SelectedBatchIndex] = true;

                                } while (nonEmptyBatches[op5SelectedBatchIndex].Family != selectedVirtualBatchIndex5 &&
                                         selectedBatchOP5.Count(item => item) < nonEmptyBatches.Count);

                                if (nonEmptyBatches[op5SelectedBatchIndex].Family != selectedVirtualBatchIndex5)
                                    break;

                                if (nonEmptyBatches[op5SelectedBatchIndex].JobsIndice.Count == 0)
                                {
                                    break;
                                }
                                double avgSjOfSelectedVirBatchOP5 =
                                    nonEmptyBatches[op5SelectedBatchIndex].SizeOfJobs.Average();

                                int op5n = r.Next(1, Convert.ToInt32(Math.Ceiling(kMax / avgSjOfSelectedVirBatchOP5)));

                                for (int j = 0;
                                    j < Math.Min(virtualBatches[selectedVirtualBatchIndex5].JobsIndice.Count, op5n);
                                    j++)
                                {
                                    int jobIndex = r.Next(virtualBatches[selectedVirtualBatchIndex5].JobsIndice.Count);
                                    int job = virtualBatches[selectedVirtualBatchIndex5].JobsIndice[jobIndex];
                                    int jobSize = virtualBatches[selectedVirtualBatchIndex5].SizeOfJobs[jobIndex];

                                    if (jobSize +
                                        nonEmptyBatches[op5SelectedBatchIndex].SizeOfJobs.Sum() <= kMax)
                                    {
                                        nonEmptyBatches[op5SelectedBatchIndex].JobsIndice.Add(job);
                                        nonEmptyBatches[op5SelectedBatchIndex].SizeOfJobs.Add(Sj[job]);
                                        nonEmptyBatches[op5SelectedBatchIndex].UrgentMetric.Add(1);


                                        double maxP1 = p1[job];

                                        foreach (int t in nonEmptyBatches[op5SelectedBatchIndex].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatches[op5SelectedBatchIndex].Pbs[0] = maxP1;

                                        double maxP2 = p2[job];

                                        foreach (int t in nonEmptyBatches[op5SelectedBatchIndex].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatches[op5SelectedBatchIndex].Pbs[1] = maxP2;

                                        selectedJobs[job] = true;

                                        virtualBatches[selectedVirtualBatchIndex5].JobsIndice.RemoveAt(jobIndex);
                                        virtualBatches[selectedVirtualBatchIndex5].SizeOfJobs.RemoveAt(jobIndex);

                                    }

                                }



                                #endregion Add 

                                break;
                            case 5:

                                #region OP5 Remove Random Batch

                                bool stopOP6Flag = nonEmptyBatches.Count == 0;

                                if (stopOP6Flag) break;

                                int selectedBatchIndex6 = r.Next(nonEmptyBatches.Count);

                                int selectedBatchFamily6 = nonEmptyBatches[selectedBatchIndex6].Family;

                                int selectedBatchLength6 = nonEmptyBatches[selectedBatchIndex6].JobsIndice.Count;

                                foreach (var item in nonEmptyBatches[selectedBatchIndex6].JobsIndice)
                                {
                                    #region Add Selected Jobs to Virtual Batch

                                    virtualBatches[selectedBatchFamily6].JobsIndice.Add(item);

                                    virtualBatches[selectedBatchFamily6].SizeOfJobs.Add(Sj[item]);

                                    //virtualBatches[selectedBatchFamily6].SizeOfJobs.Sum() += Sj[item];

                                    virtualBatches[selectedBatchFamily6].Family = selectedBatchFamily6;

                                    virtualBatches[selectedBatchFamily6].UrgentMetric.Add(1);

                                    #endregion

                                    selectedJobs[item] = false;
                                }

                                nonEmptyBatches.RemoveAt(selectedBatchIndex6);

                                Batch NewNonEmptyBatch = new Batch();

                                List<Batch> myNewNonEmptyBatches = new List<Batch>();

                                for (int j = 0; j < nonEmptyBatches.Count; j++)
                                {
                                    NewNonEmptyBatch = nonEmptyBatches[j];

                                    NewNonEmptyBatch.batchIndex = j;

                                    myNewNonEmptyBatches.Add(NewNonEmptyBatch);
                                }


                                nonEmptyBatches = myNewNonEmptyBatches;

                                #endregion Remove Bandom Batch

                                break;
                            case 6:

                                #region OP6 Transform from one nonEmpty to another

                                bool stopOP6flag = nonEmptyBatches.Count < 2;

                                if (stopOP6flag) break;

                                bool[] selectbatchesOP6 = new bool[nonEmptyBatches.Count];

                                int selectedBatchIndex1OP6 = r.Next(nonEmptyBatches.Count);

                                selectbatchesOP6[selectedBatchIndex1OP6] = true;

                                int selectedFamilyOP6 = nonEmptyBatches[selectedBatchIndex1OP6].Family;

                                int selectedBatchIndex2OP6 = -1;


                                do
                                {
                                    selectedBatchIndex2OP6 = r.Next(nonEmptyBatches.Count);

                                    if (!selectbatchesOP6[selectedBatchIndex2OP6])
                                        selectbatchesOP6[selectedBatchIndex2OP6] = true;

                                } while ((selectedBatchIndex1OP6 == selectedBatchIndex2OP6 ||
                                          selectedFamilyOP6 != nonEmptyBatches[selectedBatchIndex2OP6].Family) &&
                                         selectbatchesOP6.Count(item => item) < nonEmptyBatches.Count);

                                if (selectedFamilyOP6 != nonEmptyBatches[selectedBatchIndex2OP6].Family) break;

                                int selectedBatchLength1OP6 = nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice.Count;

                                int c = r.Next(1, Math.Max(1, selectedBatchLength1OP6 - kMin));

                                for (int j = 0; j < c; j++)
                                {
                                    int jobIndex1 = r.Next(selectedBatchLength1OP6);

                                    if ((nonEmptyBatches[selectedBatchIndex2OP6].SizeOfJobs.Sum() +
                                             nonEmptyBatches[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1] <= kMax))
                                    {
                                        int job1 = nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice[jobIndex1];

                                        int jobSize1 = nonEmptyBatches[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1];

                                        nonEmptyBatches[selectedBatchIndex2OP6].JobsIndice.Add(job1);

                                        nonEmptyBatches[selectedBatchIndex2OP6].SizeOfJobs.Add(jobSize1);

                                        nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice.RemoveAt(jobIndex1);

                                        nonEmptyBatches[selectedBatchIndex1OP6].SizeOfJobs.RemoveAt(jobIndex1);

                                        double maxP1 = p1[nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatches[selectedBatchIndex1OP6].Pbs[0] = maxP1;

                                        double maxP2 = p2[nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex1OP6].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatches[selectedBatchIndex1OP6].Pbs[1] = maxP2;

                                        double maxP12 = p1[nonEmptyBatches[selectedBatchIndex2OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex2OP6].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        nonEmptyBatches[selectedBatchIndex2OP6].Pbs[0] = maxP12;

                                        double maxP22 = p2[nonEmptyBatches[selectedBatchIndex2OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatches[selectedBatchIndex2OP6].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        nonEmptyBatches[selectedBatchIndex2OP6].Pbs[1] = maxP22;

                                    }

                                }

                                #endregion Transform from one nonEmpty to another

                                break;
                            case 7:

                                #region OP7 Change nonEmpty batchindice

                                bool stopOP8flag = nonEmptyBatches.Count < 2;

                                if (stopOP8flag) break;

                                Batch[] batchesForExchangeBatchIndice = myNonEmptyBatches;

                                int OP8selectedBatchIndex1 = r.Next(batchesForExchangeBatchIndice.Length);

                                int OP8selectedBatchIndex2 = -1;

                                do
                                    OP8selectedBatchIndex2 = r.Next(batchesForExchangeBatchIndice.Length);
                                while (OP8selectedBatchIndex1 == OP8selectedBatchIndex2);

                                int temp_BatchIndex = batchesForExchangeBatchIndice[OP8selectedBatchIndex1].batchIndex;

                                batchesForExchangeBatchIndice[OP8selectedBatchIndex1].batchIndex =
                                  batchesForExchangeBatchIndice[OP8selectedBatchIndex2].batchIndex;

                                batchesForExchangeBatchIndice[OP8selectedBatchIndex2].batchIndex = temp_BatchIndex;

                                nonEmptyBatches = batchesForExchangeBatchIndice.ToList();

                                #endregion Change nonEmpty batchindice

                                break;
                        }

                        for (int j = 0; j < selectedJobs.Length; j++)
                            PiJ[j] = !selectedJobs[j] ? 1 : 0;

                        model.IgnoranceBinary = PiJ;

                        model.NumberOfNonEmptyBatches = nonEmptyBatches.Count;

                        t1 = new double[model.NumberOfMachinesInStep1];

                        t2 = new double[model.NumberOfMachinesInStep2];

                        Tj = new double[N];

                        sol = Algorithm1(5, nonEmptyBatches, t1, t2, Tj, d, t_now);

                        model.DelayOfJobs = sol.Tj;

                        double cost = CostFunction(model);

                        if (cost < ant[k].Cost)
                        {
                            ant[k].Cost = cost;

                            ant[k].sol = sol;

                            ant[k].SelectedJobs = selectedJobs;

                            ant[k].SelectedJobsForEmptyBatches = selectedJobsForEmptyBatches;
                        }
                        if (ant[k].Cost < bestAnt.Cost)
                        {
                            bestAnt = ant[k];
                        }

                    }

                    for (int j = 0; j < selectedJobs.Length; j++)
                        R[j] = !selectedJobs[j] ? R[j] + 1 : R[j];

                    for (int j = 0; j < tauJ.Length; j++)
                        tauJ[j] = (double)1 / (double)(R[j] + 1);

                    foreach (var nonemptybatch in nonEmptyBatches)
                    {
                        for (int l = 0; l < nonemptybatch.JobsIndice.Count; l++)
                        {
                            //double sumOfPhiJXs = 0;

                            for (int j = 0; j < nonemptybatch.JobsIndice.Count; j++)
                            {
                                if (l != j)
                                {
                                    //sumOfPhiJXs += phiJX[nonemptybatch.JobsIndice[l], nonemptybatch.JobsIndice[j]];

                                    mJX[nonemptybatch.JobsIndice[l], nonemptybatch.JobsIndice[j]]++;
                                }
                            }

                            // tauJB[nonemptybatch.JobsIndice[l]] = sumOfPhiJXs / nonemptybatch.JobsIndice.Count;
                        }
                    }

                }

                #region Update Phromone

                for (int i = 0; i < phiJX.GetLength(0); i++)
                {
                    for (int j = 0; j < phiJX.GetLength(1); j++)
                    {
                        phiJX[i, j] +=
                            ((double)mJX[i, j] * ((double)Q / (double)bestAnt.Cost));
                    }
                }

                for (int i = 0; i < phiJX.GetLength(0); i++)
                {
                    for (int j = 0; j < phiJX.GetLength(1); j++)
                    {
                        phiJX[i, j] *= (double)(1 - rho);
                    }
                }

                for (int j = 0; j < tauJ.Length; j++)
                    tauJ[j] += (rho * (double)1 / (double)(R[j] + 1));

                for (int j = 0; j < tauJ.Length; j++)
                    tauJ[j] *= (double)(1 - rho);

                #endregion

                #region Show Results

                Console.WriteLine("{0}\t{1}", it, bestAnt.Cost);

                //if (it == maxIteration - 1)
                //{
                // Console.WriteLine("bestAnt.z= {0}", bestAnt.sol.z);

                // foreach (int item in bestAnt.Tour)
                // {
                //    Console.Write("{0}, ", item);
                //  }
                //  Console.WriteLine();

                // }

                streamWriter.WriteLine("{0},{1},", it, bestAnt.Cost);

                #endregion

            }


            #endregion

            streamWriter.Flush();

            streamWriter.Close();

            Console.WriteLine("done!");
            sw.Stop();
            Console.WriteLine("Elapsed Time in Milliseconds: " + sw.ElapsedMilliseconds.ToString());
        }

        static void Main(string[] args)
        {

            Console.Write("Enter the File Path: ");

            string pathToExcelFile = "D:\\129.xls";
            //string pathToExcelFile = Console.ReadLine();


            while ((!File.Exists(pathToExcelFile)))
            {
                Console.WriteLine("File does not exist! ");

                Console.Write("Enter the File Path: ");

                pathToExcelFile = Console.ReadLine();
            }

            Run_ACO(pathToExcelFile);


            Console.ReadLine();
        }
    }
}
