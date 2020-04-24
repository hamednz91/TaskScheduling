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
using Remotion.FunctionalProgramming;

namespace TaskScheduling
{
    class Program
    {
        class Sol
        {

            public List<Batch> BatchesAllocatedToMachines;

            public double[] TimeofMachinesStep1;

            public double[] TimeofMachinesStep2;

            public double[] Tj;
        }

        class Model
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

        class Batch
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

            public List<double> DueTime; //d[j] Due time of jobs

            //public int SumOfJobSizes;

            public int[] machineNumber; //machineNumber[in_step1,in_step2]

            public int[,] Pjbs; //job j proccessing time of batch b in step s (1 and 2 [step1,step2] )

            public double[] Pbs; //batch proccessing time in step 1 and 2 [step1,step2]

            public int Family;

            public int batchIndex;

            public double AverageDueTimeofJobToDelayImportanceFactor; // dj/wj


        }
        class Ant
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
        /// <param name="pathToExcelFile"> the input file directory
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

            int maxIteration = 5000;

            int numberOfAnts = 50;

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

                        batches[i].DueTime = new List<double>();

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

                        virtualBatches[i].DueTime = new List<double>();

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

                            batches[b].DueTime = new List<double>();

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

                                batches[b].DueTime.Add(d[batches[b].JobsIndice[j]]);

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
                        myNonEmptyBatches[x].batchIndex = x;
                    }


                    nonEmptyBatches = myNonEmptyBatches.ToList();

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
                        virtualBatches[j].DueTime = new List<double>();
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

                                    virtualBatches[i].DueTime.Add(d[j]);

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


                    //List<Batch> nonEmptyBatchesAfterOPs = nonEmptyBatches;
                    List<Batch> nonEmptyBatchesAfterOPs = new List<Batch>();

                    foreach (var item in nonEmptyBatches)
                        nonEmptyBatchesAfterOPs.Add(item);


                    Model modelAfterOPs = model;
                    bool[] selectedJobsAfterOPs = selectedJobs;
                    Sol solAfterOPs = sol;

                    for (int i = 0; i < A; i++)
                    {
                        List<Batch> BatchesGreaterThanKmin = new List<Batch>();

                        foreach (var item in nonEmptyBatchesAfterOPs)
                            if (item.JobsIndice.Count > kMin)
                                BatchesGreaterThanKmin.Add(item);

                        int opSelector = RouletteWheelSelection(new[] { .1, .1, .0, .1, .1, .1, .1, .1, .0, .1, .2 });

                        //opSelector = 10;

                        switch (opSelector)
                        {
                            case 0:

                                #region OP0 Drop

                                bool stopflagOP0 = BatchesGreaterThanKmin.Count == 0;

                                if (stopflagOP0) break;

                                int selectedBatchIndex = r.Next(BatchesGreaterThanKmin.Count);

                                //selectedBatchIndex = BatchesGreaterThanKmin[selectedBatchIndex].batchIndex;
                                // real index in nonEmptyBatchesAfterOPs using batch index field

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
                                        //    Sj[nonEmptyBatchesAfterOPs[selectedBatchIndex].JobsIndice[jobIndex]];

                                        virtualBatches[selectedBatchFamily].Family = selectedBatchFamily;

                                        virtualBatches[selectedBatchFamily].UrgentMetric.Add(1);

                                        virtualBatches[selectedBatchFamily].DueTime.Add(
                                            d[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]]);

                                        #endregion

                                        #region Remove Selected Jobs from NonEmptyBatches and Update Batch (size of jobs, Urgent metric, pbs if needed)

                                        bool flag1 =
                                            p1[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] >=
                                            BatchesGreaterThanKmin[selectedBatchIndex].Pbs[0];

                                        bool flag2 =
                                            p2[BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] >=
                                            BatchesGreaterThanKmin[selectedBatchIndex].Pbs[1];

                                        selectedJobsAfterOPs[
                                            BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice[jobIndex]] = false;

                                        BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice.RemoveAt(jobIndex);

                                        BatchesGreaterThanKmin[selectedBatchIndex].SizeOfJobs.RemoveAt(jobIndex);

                                        BatchesGreaterThanKmin[selectedBatchIndex].UrgentMetric.RemoveAt(jobIndex);

                                        BatchesGreaterThanKmin[selectedBatchIndex].DueTime.RemoveAt(jobIndex);

                                        selectedBatchLength =
                                            BatchesGreaterThanKmin[selectedBatchIndex].JobsIndice.Count;

                                        if (selectedBatchLength <= 0) continue;

                                        if (flag1)
                                        {
                                            double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex].JobsIndice)
                                                if (p1[t] > maxP1)
                                                    maxP1 = p1[t];

                                            nonEmptyBatchesAfterOPs[selectedBatchIndex].Pbs[0] = maxP1;
                                        }
                                        if (flag2)
                                        {
                                            double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex].JobsIndice)
                                                if (p2[t] > maxP2)
                                                    maxP2 = p2[t];

                                            nonEmptyBatchesAfterOPs[selectedBatchIndex].Pbs[1] = maxP2;
                                        }

                                        #endregion
                                    }
                                }

                                #endregion

                                break;
                            case 1:

                                #region OP1 Exchange between 2 nonempty

                                bool stopflagOP1 = nonEmptyBatchesAfterOPs.Count < 2;

                                if (stopflagOP1) break;

                                bool[] selectbatches2 = new bool[nonEmptyBatchesAfterOPs.Count];

                                int selectedBatchIndex1 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                selectbatches2[selectedBatchIndex1] = true;

                                int selectedFamily = nonEmptyBatchesAfterOPs[selectedBatchIndex1].Family;

                                int numberOfSameFamilyBatches2 =
                                    nonEmptyBatchesAfterOPs.Count(item => item.Family == selectedFamily);

                                int selectedBatchIndex2 = -1;


                                do
                                {
                                    selectedBatchIndex2 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectbatches2[selectedBatchIndex2])
                                        selectbatches2[selectedBatchIndex2] = true;

                                } while ((selectedBatchIndex1 == selectedBatchIndex2 ||
                                          selectedFamily != nonEmptyBatchesAfterOPs[selectedBatchIndex2].Family) &&
                                         selectbatches2.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (selectedFamily != nonEmptyBatchesAfterOPs[selectedBatchIndex2].Family) break;

                                int selectedBatchLength1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice.Count;
                                int selectedBatchLength2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice.Count;

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

                                                selectedBatchIndex2 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                                // in batch entekhab nashode bashe
                                            } while (selectbatches2[selectedBatchIndex2] ||
                                                     selectedBatchIndex1 == selectedBatchIndex2 ||
                                                     selectedFamily !=
                                                     nonEmptyBatchesAfterOPs[selectedBatchIndex2].Family);

                                            if (!selectbatches2[selectedBatchIndex2])
                                                selectbatches2[selectedBatchIndex2] = true;

                                            selectedBatchLength2 =
                                                nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice.Count;

                                            selectJobIndex2 = new bool[selectedBatchLength2];
                                        }

                                        do jobIndex2 = r.Next(selectedBatchLength2); while (selectJobIndex2[jobIndex2]);

                                        int a1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs[jobIndex1];

                                        batchSize1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs[jobIndex2] +
                                            a1;

                                        int a2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs[jobIndex2];

                                        batchSize2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs[jobIndex1] +
                                            a2;

                                        if (!selectJobIndex2[jobIndex2])
                                            selectJobIndex2[jobIndex2] = true;

                                        numberOfSelectedBatches2 = selectbatches2.Count(item => item);

                                    } while (numberOfSelectedBatches2 < numberOfSameFamilyBatches2 - 1
                                             && (batchSize1 > kMax || batchSize2 > kMax));

                                    if (numberOfSelectedBatches2 < numberOfSameFamilyBatches2 - 1)
                                    {

                                        int job1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice[jobIndex1];
                                        int job2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice[jobIndex2];

                                        int jobSize1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs[jobIndex1];
                                        int jobSize2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs[jobIndex2];

                                        double urgentMetric1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1].UrgentMetric[jobIndex1];
                                        double urgentMetric2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2].UrgentMetric[jobIndex2];

                                        double dueTime1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1].DueTime[jobIndex1];
                                        double dueTime2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2].DueTime[jobIndex2];


                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice.Add(job2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice.Add(job1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs.Add(jobSize2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs.Add(jobSize1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].UrgentMetric.Add(urgentMetric2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].UrgentMetric.Add(urgentMetric1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].DueTime.Add(dueTime2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].DueTime.Add(dueTime1);


                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].SizeOfJobs.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].SizeOfJobs.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].UrgentMetric.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].UrgentMetric.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].DueTime.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].DueTime.RemoveAt(jobIndex2);


                                        double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].Pbs[0] = maxP1;

                                        double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1].Pbs[1] = maxP2;


                                        double maxP12 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].Pbs[0] = maxP12;

                                        double maxP22 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2].Pbs[1] = maxP22;
                                    }
                                }

                                #endregion

                                break;
                            case 2:

                                #region OP2 Create new batch from VirtualBatches

                                bool stopFlagOP2 = virtualBatches.All(vb => vb.JobsIndice.Count < kMin);

                                if (stopFlagOP2) break;

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
                                newBatch.DueTime = new List<double>();
                                newBatch.Pbs = new double[2]; //batch processing time in step 1 & 2
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

                                bool[] selectedJobFromSelectedVirtualBatchOP2 = new bool[selectedVirBatchLength];

                                for (int j = 0; j < Math.Min(nm, selectedVirBatchLength); j++)
                                {
                                    do
                                    {
                                        jobIndexOP2 = r.Next(selectedVirBatchLength);

                                        jobOP2 = VirtualBatchesMoreThanKmin[virtualBatchIndex].JobsIndice[jobIndexOP2];

                                    } while (selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] &&
                                             selectedJobFromSelectedVirtualBatchOP2.Count(item => item) <
                                             selectedVirBatchLength);

                                    if (!selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] &&
                                        Sj[jobOP2] + newBatch.SizeOfJobs.Sum() < kMax)
                                    {
                                        newBatch.JobsIndice.Add(jobOP2);
                                        newBatch.SizeOfJobs.Add(Sj[jobOP2]);
                                        newBatch.UrgentMetric.Add(1);
                                        newBatch.DueTime.Add(d[jobOP2]);
                                        newBatch.Family = virtualBatchIndex;

                                        if (!selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2])
                                            selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] = true;

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

                                        VirtualBatchesMoreThanKmin[virtualBatchIndex].JobsIndice.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndex].SizeOfJobs.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndex].UrgentMetric.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndex].DueTime.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin.ToArray()[virtualBatchIndex].Family = -1;

                                        selectedVirBatchLength =
                                            VirtualBatchesMoreThanKmin[virtualBatchIndex].JobsIndice.Count;
                                    }
                                }


                                //Batch[] nonEmptyBatchesAverageDjToWjOP2 = nonEmptyBatchesAfterOPs.ToArray();
                                //List<Batch> nonEmptyBatchesAverageDjToWjOP2 = nonEmptyBatchesAfterOPs;
                                List<Batch> nonEmptyBatchesAverageDjToWjOP2 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                    nonEmptyBatchesAverageDjToWjOP2.Add(item);


                                // if (selectedFamilyVirtualBatch.Count(item => item) >= VirtualBatchesMoreThanKmin.Count) break;

                                if (newBatch.JobsIndice.Count >= kMin)
                                {
                                    double average = 0;
                                    foreach (var item in newBatch.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        selectedJobsAfterOPs[item] = true;

                                    }
                                    average = average / (double)newBatch.JobsIndice.Count;

                                    newBatch.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    int l = 0;

                                    //for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP2.Length; j++)
                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP2.Count; j++)
                                    {
                                        if (average >
                                            nonEmptyBatchesAverageDjToWjOP2[j].AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            newBatch.batchIndex = j + 1;

                                            //nonEmptyBatchesAverageDjToWjOP2.ToList().Insert(j + 1, newBatch);
                                            nonEmptyBatchesAverageDjToWjOP2.Insert(j + 1, newBatch);

                                            isSet = true;

                                            l = j + 2;

                                            break;
                                        }
                                    }

                                    if (isSet)
                                        //for (; l < nonEmptyBatchesAverageDjToWjOP2.Length; l++)
                                        for (; l < nonEmptyBatchesAverageDjToWjOP2.Count; l++)
                                            //nonEmptyBatchesAverageDjToWjOP2[l].batchIndex++;
                                            nonEmptyBatchesAverageDjToWjOP2[l].batchIndex++;
                                    else
                                        //nonEmptyBatchesAverageDjToWjOP2.ToList().Insert(0, newBatch);
                                        nonEmptyBatchesAverageDjToWjOP2.Insert(0, newBatch);


                                }

                                //nonEmptyBatchesAfterOPs = nonEmptyBatchesAverageDjToWjOP2.ToList();
                                nonEmptyBatchesAfterOPs = nonEmptyBatchesAverageDjToWjOP2;

                                #endregion

                                break;
                            case 3:

                                #region OP3 Exchange between virtual batch and nonempty

                                bool stopFlagOP3 = virtualBatches.All(vb => vb.JobsIndice.Count == 0) ||
                                                   nonEmptyBatchesAfterOPs.Count == 0;

                                if (stopFlagOP3) break;

                                int selectedBatchIndexOP3;

                                int selectedBatchFamilyOP3;

                                bool changeSelectedBatchFlagOP3 = false;

                                bool[] selectedBatchOP3 = new bool[nonEmptyBatchesAfterOPs.Count];

                                do
                                {
                                    selectedBatchIndexOP3 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    selectedBatchFamilyOP3 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Family;

                                    changeSelectedBatchFlagOP3 =
                                        virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count == 0;

                                    if (!selectedBatchOP3[selectedBatchIndexOP3])
                                        selectedBatchOP3[selectedBatchIndexOP3] = true;

                                } while (changeSelectedBatchFlagOP3 &&
                                         selectedBatchOP3.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (selectedBatchOP3.Count(item => item) >= nonEmptyBatchesAfterOPs.Count) break;


                                int selectedBatchLengthOP3 =
                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice.Count;


                                int op4n = r.Next(1,
                                    Math.Max(1,
                                        Math.Min(selectedBatchLengthOP3,
                                            virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count)));

                                for (int j = 0; j < op4n; j++)
                                {
                                    int batchSize1 = 0;

                                    int jobIndex1 = r.Next(selectedBatchLengthOP3);

                                    int jobIndex2 = -1;

                                    bool[] selectJobIndex2 =
                                        new bool[virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count];

                                    int numberofSelectedJobfromBatch4 = 0;

                                    do
                                    {
                                        //check if jobindex2 is already chosen
                                        numberofSelectedJobfromBatch4 = selectJobIndex2.Count(item => item);

                                        jobIndex2 = r.Next(virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count);

                                        int a1 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs[jobIndex1];

                                        batchSize1 = virtualBatches[selectedBatchFamilyOP3].SizeOfJobs[jobIndex2] +
                                                     a1;

                                        if (!selectJobIndex2[jobIndex2])
                                            selectJobIndex2[jobIndex2] = true;

                                    } while (batchSize1 > kMax && numberofSelectedJobfromBatch4 <
                                             virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count);

                                    if (numberofSelectedJobfromBatch4 >=
                                        virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count)
                                        break;

                                    int job1 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[jobIndex1];
                                    int job2 = virtualBatches[selectedBatchFamilyOP3].JobsIndice[jobIndex2];

                                    int jobSize1 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs[jobIndex1];
                                    int jobSize2 = virtualBatches[selectedBatchFamilyOP3].SizeOfJobs[jobIndex2];


                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice.Add(job2);
                                    virtualBatches[selectedBatchFamilyOP3].JobsIndice.Add(job1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.Add(jobSize2);
                                    virtualBatches[selectedBatchFamilyOP3].SizeOfJobs.Add(jobSize1);

                                    bool flag1 = p1[job1] >=
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[0];

                                    bool flag2 = p2[job1] >=
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[1];


                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice.RemoveAt(jobIndex1);
                                    virtualBatches[selectedBatchFamilyOP3].JobsIndice.RemoveAt(jobIndex2);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.RemoveAt(jobIndex1);
                                    virtualBatches[selectedBatchFamilyOP3].SizeOfJobs.RemoveAt(jobIndex2);

                                    if (flag1)
                                    {
                                        double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[0] = maxP1;
                                    }

                                    if (flag2)
                                    {
                                        double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[1] = maxP2;
                                    }

                                    selectedJobsAfterOPs[job1] = false;

                                    selectedJobsAfterOPs[job2] = true;

                                }


                                #endregion

                                break;
                            case 4:

                                #region OP4 Add

                                bool stopOP5Flag = virtualBatches.All(vb => vb.JobsIndice.Count == 0) ||
                                                   nonEmptyBatchesAfterOPs.Count == 0;

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

                                bool[] selectedBatchOP5 = new bool[nonEmptyBatchesAfterOPs.Count];

                                do
                                {
                                    op5SelectedBatchIndex = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectedBatchOP5[op5SelectedBatchIndex])
                                        selectedBatchOP5[op5SelectedBatchIndex] = true;

                                } while (nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].Family !=
                                         selectedVirtualBatchIndex5 &&
                                         selectedBatchOP5.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].Family != selectedVirtualBatchIndex5)
                                    break;

                                if (nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].JobsIndice.Count == 0)
                                {
                                    break;
                                }
                                double avgSjOfSelectedVirBatchOP5 =
                                    nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].SizeOfJobs.Average();

                                int op5n = r.Next(1, Convert.ToInt32(Math.Ceiling(kMax / avgSjOfSelectedVirBatchOP5)));

                                for (int j = 0;
                                    j < Math.Min(virtualBatches[selectedVirtualBatchIndex5].JobsIndice.Count, op5n);
                                    j++)
                                {
                                    int jobIndex = r.Next(virtualBatches[selectedVirtualBatchIndex5].JobsIndice.Count);
                                    int job = virtualBatches[selectedVirtualBatchIndex5].JobsIndice[jobIndex];
                                    int jobSize = virtualBatches[selectedVirtualBatchIndex5].SizeOfJobs[jobIndex];

                                    if (jobSize +
                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].SizeOfJobs.Sum() <= kMax)
                                    {
                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].JobsIndice.Add(job);
                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].SizeOfJobs.Add(Sj[job]);
                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].UrgentMetric.Add(1);
                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].DueTime.Add(d[job]);


                                        double maxP1 = p1[job];

                                        foreach (int t in nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].Pbs[0] = maxP1;

                                        double maxP2 = p2[job];

                                        foreach (int t in nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[op5SelectedBatchIndex].Pbs[1] = maxP2;

                                        selectedJobsAfterOPs[job] = true;

                                        virtualBatches[selectedVirtualBatchIndex5].JobsIndice.RemoveAt(jobIndex);
                                        virtualBatches[selectedVirtualBatchIndex5].SizeOfJobs.RemoveAt(jobIndex);

                                    }

                                }



                                #endregion Add 

                                break;
                            case 5:

                                #region OP5 Remove Random Batch

                                bool stopOP6Flag = nonEmptyBatchesAfterOPs.Count == 0;

                                if (stopOP6Flag) break;

                                int selectedBatchIndex6 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                int selectedBatchFamily6 = nonEmptyBatchesAfterOPs[selectedBatchIndex6].Family;

                                int selectedBatchLength6 = nonEmptyBatchesAfterOPs[selectedBatchIndex6].JobsIndice.Count;

                                foreach (var item in nonEmptyBatchesAfterOPs[selectedBatchIndex6].JobsIndice)
                                {
                                    #region Add Selected Jobs to Virtual Batch

                                    virtualBatches[selectedBatchFamily6].JobsIndice.Add(item);

                                    virtualBatches[selectedBatchFamily6].SizeOfJobs.Add(Sj[item]);

                                    //virtualBatches[selectedBatchFamily6].SizeOfJobs.Sum() += Sj[item];

                                    virtualBatches[selectedBatchFamily6].Family = selectedBatchFamily6;

                                    virtualBatches[selectedBatchFamily6].UrgentMetric.Add(1);

                                    virtualBatches[selectedBatchFamily6].DueTime.Add(d[item]);

                                    #endregion

                                    selectedJobsAfterOPs[item] = false;
                                }

                                nonEmptyBatchesAfterOPs.RemoveAt(selectedBatchIndex6);

                                Batch NewNonEmptyBatch = new Batch();

                                List<Batch> myNewNonEmptyBatches = new List<Batch>();

                                for (int j = 0; j < nonEmptyBatchesAfterOPs.Count; j++)
                                {
                                    NewNonEmptyBatch = nonEmptyBatchesAfterOPs[j];

                                    NewNonEmptyBatch.batchIndex = j;

                                    myNewNonEmptyBatches.Add(NewNonEmptyBatch);
                                }


                                nonEmptyBatchesAfterOPs = myNewNonEmptyBatches;

                                #endregion Remove Bandom Batch

                                break;
                            case 6:

                                #region OP6 Transform from one nonEmpty to another

                                bool stopOP6flag = nonEmptyBatchesAfterOPs.Count < 2 ||
                                                   nonEmptyBatchesAfterOPs.All(item => item.JobsIndice.Count <= kMin);

                                if (stopOP6flag) break;

                                bool[] selectbatchesOP6 = new bool[BatchesGreaterThanKmin.Count];

                                int selectedBatchIndex1OP6;

                                do
                                {
                                    selectedBatchIndex1OP6 = r.Next(BatchesGreaterThanKmin.Count);

                                    if (!selectbatchesOP6[selectedBatchIndex1OP6])
                                    {
                                        selectbatchesOP6[selectedBatchIndex1OP6] = true;
                                    }

                                } while (selectbatchesOP6.Count(item => item) < BatchesGreaterThanKmin.Count);

                                if (selectbatchesOP6.Count(item => item) >= BatchesGreaterThanKmin.Count) break;

                                int selectedFamilyOP6 = BatchesGreaterThanKmin[selectedBatchIndex1OP6].Family;

                                int selectedBatchIndex2OP6;


                                do
                                {
                                    selectedBatchIndex2OP6 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectbatchesOP6[selectedBatchIndex2OP6])
                                        selectbatchesOP6[selectedBatchIndex2OP6] = true;

                                } while ((selectedBatchIndex1OP6 == selectedBatchIndex2OP6 ||
                                          selectedFamilyOP6 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].Family) &&
                                         nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Sum() >= kMax &&
                                         selectbatchesOP6.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (selectedFamilyOP6 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].Family) break;

                                int selectedBatchLength1OP6 =
                                    BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice.Count;

                                double minMetric =
                                    (kMax - nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Sum()) /
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Average();

                                // minMetric wouldn`t be integer, and the Minimum would be double as well. is it correct to convert to integer?
                                int c = r.Next(1, (int)Math.Min(minMetric, (selectedBatchLength1OP6 - kMin)));

                                for (int j = 0; j < c; j++)
                                {
                                    int jobIndex1 = r.Next(selectedBatchLength1OP6);

                                    if ((nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Sum() +
                                         BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1] <= kMax))
                                    {
                                        int job1 = BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice[jobIndex1];

                                        int jobSize1 =
                                            BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice.Add(job1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Add(jobSize1);

                                        BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice.RemoveAt(jobIndex1);

                                        BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs.RemoveAt(jobIndex1);

                                        double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].Pbs[0] = maxP1;

                                        double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP6].Pbs[1] = maxP2;

                                        double maxP12 =
                                            p1[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].Pbs[0] = maxP12;

                                        double maxP22 =
                                            p2[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].Pbs[1] = maxP22;

                                    }

                                }

                                #endregion Transform from one nonEmpty to another

                                break;
                            case 7:

                                #region OP7 Change nonEmpty batchindice

                                bool stopflagOP7 = nonEmptyBatchesAfterOPs.Count < 2;

                                if (stopflagOP7) break;

                                Batch[] batchesForExchangeBatchIndice = myNonEmptyBatches;

                                int selectedBatchIndex1OP7 = r.Next(batchesForExchangeBatchIndice.Length);

                                int selectedBatchIndex2OP7 = -1;

                                do
                                    selectedBatchIndex2OP7 = r.Next(batchesForExchangeBatchIndice.Length); while (
                                    selectedBatchIndex1OP7 == selectedBatchIndex2OP7);

                                int temp_BatchIndex = batchesForExchangeBatchIndice[selectedBatchIndex1OP7].batchIndex;

                                batchesForExchangeBatchIndice[selectedBatchIndex1OP7].batchIndex =
                                    batchesForExchangeBatchIndice[selectedBatchIndex2OP7].batchIndex;

                                batchesForExchangeBatchIndice[selectedBatchIndex2OP7].batchIndex = temp_BatchIndex;

                                nonEmptyBatches = batchesForExchangeBatchIndice.ToList();

                                #endregion Change nonEmpty batchindice

                                break;
                            case 8:

                                #region OP8 Create New Batch from nonEmptyBatchesAfterOPs

                                bool stopOP8Flag = nonEmptyBatchesAfterOPs.All(item => item.JobsIndice.Count <= 2 * kMin);

                                if (stopOP8Flag) break;

                                #region New Batch Instance

                                Batch newBatchOP8 = new Batch();

                                #endregion

                                Batch[] BatchesWithJobsGreaterthan2Kmin =
                                    nonEmptyBatchesAfterOPs.Where(item => item.JobsIndice.Count > 2 * kMin).ToArray();

                                bool[] selectedBatchOP8 = new bool[BatchesWithJobsGreaterthan2Kmin.Length];

                                #region New Batch Init/Reset

                                newBatchOP8.JobsIndice = new List<int>();
                                newBatchOP8.SizeOfJobs = new List<int>();
                                newBatchOP8.UrgentMetric = new List<double>();
                                newBatchOP8.DueTime = new List<double>();
                                newBatchOP8.Pbs = new double[2]; //batch processing time in step 1 & 2
                                newBatchOP8.Family = -1;
                                newBatchOP8.machineNumber = new int[] { -1, -1 };
                                newBatchOP8.batchIndex = -1;

                                #endregion

                                int selectedBatchIndexOP8 = r.Next(BatchesWithJobsGreaterthan2Kmin.Length);

                                if (!selectedBatchOP8[selectedBatchIndexOP8])
                                    selectedBatchOP8[selectedBatchIndexOP8] = true;

                                int selectedBatchLengthOP8 =
                                    BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice.Count;


                                int cnOP8 = r.Next(kMin, selectedBatchLengthOP8 - kMin);


                                bool[] selectedJobFromSelectedBatchOP8 = new bool[selectedBatchLengthOP8];

                                int jobOP8, jobIndexOP8;

                                for (int j = 0; j < cnOP8; j++)
                                {
                                    do
                                    {
                                        jobIndexOP8 = r.Next(selectedBatchLengthOP8);

                                        jobOP8 =
                                            BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice[
                                                jobIndexOP8];

                                    } while (selectedJobFromSelectedBatchOP8[jobIndexOP8] &&
                                             selectedJobFromSelectedBatchOP8.Count(item => item) <
                                             selectedBatchLengthOP8);

                                    if (!selectedJobFromSelectedBatchOP8[jobIndexOP8] &&
                                        Sj[jobOP8] + newBatchOP8.SizeOfJobs.Sum() < kMax)
                                    {
                                        newBatchOP8.JobsIndice.Add(jobOP8);
                                        newBatchOP8.SizeOfJobs.Add(Sj[jobOP8]);
                                        newBatchOP8.UrgentMetric.Add(1);
                                        newBatchOP8.DueTime.Add(d[jobOP8]);
                                        newBatchOP8.Family =
                                            BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].Family;

                                        if (!selectedJobFromSelectedBatchOP8[jobIndexOP8])
                                            selectedJobFromSelectedBatchOP8[jobIndexOP8] = true;

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice.RemoveAt(
                                            jobIndexOP8);

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].SizeOfJobs.RemoveAt(
                                            jobIndexOP8);

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].UrgentMetric.RemoveAt(
                                            jobIndexOP8);

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].DueTime.RemoveAt(
                                            jobIndexOP8);

                                        double maxP12 =
                                            p1[BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice[0]];

                                        foreach (
                                            int t in BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].Pbs[0] = maxP12;

                                        double maxP22 =
                                            p2[BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice[0]];

                                        foreach (
                                            int t in BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].Pbs[1] = maxP22;

                                        selectedBatchLengthOP8 =
                                            BatchesWithJobsGreaterthan2Kmin[selectedBatchIndexOP8].JobsIndice.Count;

                                        double maxP1 = p1[jobOP8];

                                        foreach (int t in newBatchOP8.JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        newBatchOP8.Pbs[0] = maxP1;

                                        double maxP2 = p2[jobOP8];

                                        foreach (int t in newBatchOP8.JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        newBatchOP8.Pbs[1] = maxP2;

                                    }
                                }


                                //Batch[] nonEmptyBatchesAverageDjToWjOP8 = nonEmptyBatchesAfterOPs.ToArray();
                                //List<Batch> nonEmptyBatchesAverageDjToWjOP8 = nonEmptyBatchesAfterOPs;

                                List<Batch> nonEmptyBatchesAverageDjToWjOP8 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                    nonEmptyBatchesAverageDjToWjOP8.Add(item);


                                // if (selectedFamilyVirtualBatch.Count(item => item) >= VirtualBatchesMoreThanKmin.Count) break;

                                if (newBatchOP8.JobsIndice.Count >= kMin)
                                {

                                    double average = 0;
                                    foreach (var item in newBatchOP8.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        // selectedJobsAfterOPs[item] = true;

                                    }
                                    average = average / (double)newBatchOP8.JobsIndice.Count;

                                    newBatchOP8.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    int l = 0;

                                    //for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP8.Length; j++)
                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP8.Count; j++)
                                    {
                                        if (average >
                                            nonEmptyBatchesAverageDjToWjOP8[j]
                                                .AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            newBatchOP8.batchIndex = j + 1;

                                            //nonEmptyBatchesAverageDjToWjOP8.ToList().Insert(j + 1, newBatchOP8);
                                            nonEmptyBatchesAverageDjToWjOP8.Insert(j + 1, newBatchOP8);

                                            isSet = true;

                                            l = j + 2;

                                            break;

                                        }
                                    }



                                    if (isSet)
                                        //for (; l < nonEmptyBatchesAverageDjToWjOP8.Length; l++)
                                        for (; l < nonEmptyBatchesAverageDjToWjOP8.Count; l++)
                                        {
                                            //nonEmptyBatchesAverageDjToWjOP8[l].batchIndex++;
                                            nonEmptyBatchesAverageDjToWjOP8[l].batchIndex++;

                                        }
                                    else
                                        //nonEmptyBatchesAverageDjToWjOP8.ToList().Insert(0, newBatchOP8);
                                        nonEmptyBatchesAverageDjToWjOP8.Insert(0, newBatchOP8);


                                }

                                //nonEmptyBatchesAfterOPs = nonEmptyBatchesAverageDjToWjOP8.ToList();
                                nonEmptyBatchesAfterOPs = nonEmptyBatchesAverageDjToWjOP8;

                                #endregion

                                break;
                            case 9:

                                #region OP9 Exchange between 2 nonempty using D[j]

                                bool stopflagOP9 = nonEmptyBatchesAfterOPs.Count < 2;

                                if (stopflagOP9) break;

                                bool[] selectedBatchesOP9 = new bool[nonEmptyBatchesAfterOPs.Count];

                                int selectedBatchIndex1OP9 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                selectedBatchesOP9[selectedBatchIndex1OP9] = true;

                                int selectedFamilyOP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].Family;

                                int numberOfSameFamilyBatches2OP9 =
                                    nonEmptyBatchesAfterOPs.Count(item => item.Family == selectedFamilyOP9);

                                int selectedBatchIndex2OP9;

                                do
                                {
                                    selectedBatchIndex2OP9 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectedBatchesOP9[selectedBatchIndex2OP9])
                                        selectedBatchesOP9[selectedBatchIndex2OP9] = true;

                                } while ((selectedBatchIndex1OP9 == selectedBatchIndex2OP9 ||
                                          selectedFamilyOP9 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].Family) &&
                                         selectedBatchesOP9.Count(item => item) < nonEmptyBatchesAfterOPs.Count - 1);

                                if (selectedFamilyOP9 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].Family) break;

                                int selectedBatchLength1OP9 =
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice.Count;
                                int selectedBatchLength2OP9 =
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.Count;

                                //
                                //can number 2 be chosen in the random selection below?
                                //
                                int nl = r.Next(1, Math.Min(selectedBatchLength1OP9, selectedBatchLength2OP9));

                                for (int j = 0; j < nl; j++)
                                {
                                    int batchSize1, batchSize2 = 0;

                                    int jobIndex1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.IndexOf(
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.Max());

                                    int jobIndex2 = -1;




                                    do
                                    {
                                        do
                                        {

                                            selectedBatchIndex2OP9 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                            if (!selectedBatchesOP9[selectedBatchIndex2OP9])
                                                selectedBatchesOP9[selectedBatchIndex2OP9] = true;

                                            // in batch entekhab nashode bashe
                                        } while ((selectedBatchesOP9.Count(item => item) <
                                                  nonEmptyBatchesAfterOPs.Count - 1) &&
                                                 (selectedBatchIndex1OP9 == selectedBatchIndex2OP9 ||
                                                  selectedFamilyOP9 !=
                                                  nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].Family));

                                        if (selectedBatchesOP9.Count(item => item) >= nonEmptyBatchesAfterOPs.Count - 1)
                                            break;

                                        selectedBatchLength2OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.Count;


                                        jobIndex2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.IndexOf(
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.Min());

                                        int a1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[jobIndex1];

                                        batchSize1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[jobIndex2] +
                                            a1;

                                        int a2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[jobIndex2];

                                        batchSize2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[jobIndex1] +
                                            a2;

                                    } while (selectedBatchesOP9.Count(item => item) < numberOfSameFamilyBatches2OP9 - 1
                                             && (batchSize1 > kMax || batchSize2 > kMax));

                                    if (selectedBatchesOP9.Count(item => item) >= nonEmptyBatchesAfterOPs.Count - 1)
                                        break;

                                    if (selectedBatchesOP9.Count(item => item) < numberOfSameFamilyBatches2OP9 - 1)
                                    {

                                        int job1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice[jobIndex1];
                                        int job2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice[jobIndex2];

                                        int jobSize1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[jobIndex1];
                                        int jobSize2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[jobIndex2];

                                        double urgentMetric1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric[jobIndex1];
                                        double urgentMetric2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric[jobIndex2];

                                        double dueTime1 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime[jobIndex1];
                                        double dueTime2 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime[jobIndex2];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice.Add(job2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.Add(job1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.Add(jobSize2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.Add(jobSize1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric.Add(urgentMetric2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric.Add(urgentMetric1);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.Add(dueTime2);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.Add(dueTime1);


                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric.RemoveAt(jobIndex2);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.RemoveAt(jobIndex1);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.RemoveAt(jobIndex2);


                                        double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].Pbs[0] = maxP1;

                                        double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].Pbs[1] = maxP2;


                                        double maxP12 =
                                            p1[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice)
                                            if (p1[t] > maxP12)
                                                maxP12 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].Pbs[0] = maxP12;

                                        double maxP22 =
                                            p2[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice)
                                            if (p2[t] > maxP22)
                                                maxP22 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].Pbs[1] = maxP22;
                                    }
                                }

                                #endregion

                                break;
                            case 10:

                                #region OP10 Create new batch from VirtualBatches

                                bool stopFlagOP10 = nonEmptyBatchesAfterOPs.All(
                                    aa => aa.JobsIndice.Count + virtualBatches[aa.Family].JobsIndice.Count < 2 * kMin);

                                if (stopFlagOP10) break;

                                #region New Batch Instance

                                Batch newBatchOP10 = new Batch();

                                #endregion

                                Batch[] BatchesWithJobsGreaterthan2KminOP10 =
                                    nonEmptyBatchesAfterOPs.Where(
                                        item =>
                                            item.JobsIndice.Count + virtualBatches[item.Family].JobsIndice.Count >=
                                            2 * kMin).ToArray();

                                bool[] selectedBatchOP10 = new bool[BatchesWithJobsGreaterthan2KminOP10.Length];

                                #region New Batch Init/Reset

                                newBatchOP10.JobsIndice = new List<int>();
                                newBatchOP10.SizeOfJobs = new List<int>();
                                newBatchOP10.UrgentMetric = new List<double>();
                                newBatchOP10.DueTime = new List<double>();
                                newBatchOP10.Pbs = new double[2]; //batch processing time in step 1 & 2
                                newBatchOP10.Family = -1;
                                newBatchOP10.machineNumber = new int[] { -1, -1 };
                                newBatchOP10.batchIndex = -1;

                                #endregion

                                int selectedBatchIndexOP10 = r.Next(BatchesWithJobsGreaterthan2KminOP10.Length);

                                if (!selectedBatchOP10[selectedBatchIndexOP10])
                                    selectedBatchOP10[selectedBatchIndexOP10] = true;

                                int selectedBatchLengthOP10 =
                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice.Count;

                                int selectedBatchFamilyOP10 =
                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].Family;

                                if (virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count == 0) break;

                                int cnOP10 = r.Next(kMin,
                                    (selectedBatchLengthOP10 - kMin) +
                                    Math.Min((int)(kMax / virtualBatches[selectedBatchFamilyOP10].SizeOfJobs.Average()),
                                        virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count));

                                int selectedVirtualBatchLengthOP10 =
                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count;


                                bool[] selectedJobFromSelectedBatchOP10 = new bool[selectedBatchLengthOP10];

                                bool[] selectedJobFromSelectedVirtualBatchOP10 = new bool[selectedVirtualBatchLengthOP10];

                                int jobOP10, jobIndexOP10, numberOfUnSelectedVirtualBatches=0;

                                for (int j = 0; j < cnOP10;)
                                {
                                    if (selectedVirtualBatchLengthOP10 > 0)
                                    {
                                        do
                                        {
                                            jobIndexOP10 = r.Next(selectedVirtualBatchLengthOP10);

                                            jobOP10 =
                                                virtualBatches[selectedBatchFamilyOP10].JobsIndice[
                                                    jobIndexOP10];

                                        } while (selectedJobFromSelectedVirtualBatchOP10[jobIndexOP10] &&
                                                 selectedJobFromSelectedVirtualBatchOP10.Count(item => item) <
                                                 selectedVirtualBatchLengthOP10);

                                        if (!selectedJobFromSelectedVirtualBatchOP10[jobIndexOP10] &&
                                            Sj[jobOP10] + newBatchOP10.SizeOfJobs.Sum() < kMax)
                                        {
                                            selectedJobFromSelectedVirtualBatchOP10[jobIndexOP10] = true;

                                            newBatchOP10.JobsIndice.Add(jobOP10);
                                            newBatchOP10.SizeOfJobs.Add(Sj[jobOP10]);
                                            newBatchOP10.UrgentMetric.Add(1);
                                            newBatchOP10.DueTime.Add(d[jobOP10]);
                                            newBatchOP10.Family = selectedBatchFamilyOP10;

                                            double maxP1 = p1[jobOP10];

                                            foreach (int t in newBatchOP10.JobsIndice)
                                                if (p1[t] > maxP1)
                                                    maxP1 = p1[t];

                                            newBatchOP10.Pbs[0] = maxP1;

                                            double maxP2 = p2[jobOP10];

                                            foreach (int t in newBatchOP10.JobsIndice)
                                                if (p2[t] > maxP2)
                                                    maxP2 = p2[t];

                                            newBatchOP10.Pbs[1] = maxP2;

                                        }
                                        if (!selectedJobFromSelectedVirtualBatchOP10[jobIndexOP10])
                                        {
                                            selectedJobFromSelectedVirtualBatchOP10[jobIndexOP10] = true;
                                            numberOfUnSelectedVirtualBatches++;
                                        }
                                    }

                                    else if (selectedBatchLengthOP10 > 0 &&
                                        (selectedJobFromSelectedVirtualBatchOP10.All(item => item)))
                                    {

                                        do
                                        {
                                            jobIndexOP10 = r.Next(selectedBatchLengthOP10);

                                            jobOP10 =
                                                BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice[
                                                    jobIndexOP10];

                                        } while (selectedJobFromSelectedBatchOP10[jobIndexOP10] &&
                                                     selectedJobFromSelectedBatchOP10.Count(item => item) <
                                                     selectedBatchLengthOP10);

                                        if (!selectedJobFromSelectedBatchOP10[jobIndexOP10] &&
                                            Sj[jobOP10] + newBatchOP10.SizeOfJobs.Sum() < kMax)
                                        {
                                            selectedJobFromSelectedBatchOP10[jobIndexOP10] = true;

                                            newBatchOP10.JobsIndice.Add(jobOP10);
                                            newBatchOP10.SizeOfJobs.Add(Sj[jobOP10]);
                                            newBatchOP10.UrgentMetric.Add(1);
                                            newBatchOP10.DueTime.Add(d[jobOP10]);
                                            newBatchOP10.Family = selectedBatchFamilyOP10;

                                            double maxP1 = p1[jobOP10];

                                            foreach (int t in newBatchOP10.JobsIndice)
                                                if (p1[t] > maxP1)
                                                    maxP1 = p1[t];

                                            newBatchOP10.Pbs[0] = maxP1;

                                            double maxP2 = p2[jobOP10];

                                            foreach (int t in newBatchOP10.JobsIndice)
                                                if (p2[t] > maxP2)
                                                    maxP2 = p2[t];

                                            newBatchOP10.Pbs[1] = maxP2;

                                        }
                                        if (!selectedJobFromSelectedBatchOP10[jobIndexOP10])
                                            selectedJobFromSelectedBatchOP10[jobIndexOP10] = true;
                                    }

                                    if (numberOfUnSelectedVirtualBatches>0)
                                        j += numberOfUnSelectedVirtualBatches;
                                    else
                                        j++;


                                }


                                List<Batch> nonEmptyBatchesAverageDjToWjOP10 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                {
                                    nonEmptyBatchesAverageDjToWjOP10.Add(item);
                                }

                                if (newBatchOP10.JobsIndice.Count >= kMin)
                                {
                                    double average = 0;
                                    foreach (var item in newBatchOP10.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        selectedJobsAfterOPs[item] = true;

                                    }
                                    average = average / (double)newBatchOP10.JobsIndice.Count;

                                    newBatchOP10.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    int l = 0;

                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP10.Count; j++)
                                    {
                                        if (average >
                                            nonEmptyBatchesAverageDjToWjOP10[j]
                                                .AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            newBatchOP10.batchIndex = j + 1;

                                            nonEmptyBatchesAverageDjToWjOP10.Insert(j + 1, newBatchOP10);

                                            isSet = true;

                                            l = j + 2;

                                            break;

                                        }
                                    }

                                    if (isSet)
                                        for (; l < nonEmptyBatchesAverageDjToWjOP10.Count; l++)
                                            nonEmptyBatchesAverageDjToWjOP10[l].batchIndex++;
                                    else
                                        nonEmptyBatchesAverageDjToWjOP10.Insert(0, newBatchOP10);

                                    nonEmptyBatchesAfterOPs = nonEmptyBatchesAverageDjToWjOP10;

                                    foreach (var job in newBatchOP10.JobsIndice)
                                    {
                                        int index;
                                        double maxP12, maxP22;

                                        #region Remove from VirtualBatches
                                        if (virtualBatches[selectedBatchFamilyOP10].JobsIndice.Any(vj => vj == job))
                                        {
                                            index =
                                                Array.IndexOf(
                                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice.ToArray(), job);

                                            virtualBatches[selectedBatchFamilyOP10].JobsIndice.RemoveAt(index);

                                            virtualBatches[selectedBatchFamilyOP10].SizeOfJobs.RemoveAt(index);

                                            virtualBatches[selectedBatchFamilyOP10].UrgentMetric.RemoveAt(index);

                                            virtualBatches[selectedBatchFamilyOP10].DueTime.RemoveAt(index);

                                            selectedVirtualBatchLengthOP10 =
                                                virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count;

                                            if (selectedVirtualBatchLengthOP10 > 0)
                                            {
                                                maxP12 =
                                                   p1[virtualBatches[selectedBatchFamilyOP10].JobsIndice[0]
                                                   ];

                                                foreach (
                                                    int t in
                                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice)
                                                    if (p1[t] > maxP12)
                                                        maxP12 = p1[t];

                                                virtualBatches[selectedBatchFamilyOP10].Pbs[0] = maxP12;

                                                maxP22 =
                                                   p2[virtualBatches[selectedBatchFamilyOP10].JobsIndice[0]
                                                   ];

                                                foreach (
                                                    int t in
                                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice)
                                                    if (p2[t] > maxP22)
                                                        maxP22 = p2[t];

                                                virtualBatches[selectedBatchFamilyOP10].Pbs[1] = maxP22;
                                            }
                                        }
                                        #endregion

                                        #region Remove from BatchesWithJobsGreaterthan2KminOP10
                                        else
                                        {
                                            index =
                                                Array.IndexOf(
                                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                        .JobsIndice.ToArray(), job);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice
                                                .RemoveAt(
                                                    index);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].SizeOfJobs
                                                .RemoveAt(
                                                    index);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].UrgentMetric
                                                .RemoveAt(
                                                    index);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].DueTime.RemoveAt
                                            (
                                                index);

                                            selectedBatchLengthOP10 =
                                                BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice
                                                    .Count;

                                            if (selectedBatchLengthOP10 > 0)
                                            {

                                                maxP12 =
                                                    p1[
                                                        BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                            .JobsIndice[0]];

                                                foreach (
                                                    int t in
                                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                        .JobsIndice)
                                                    if (p1[t] > maxP12)
                                                        maxP12 = p1[t];

                                                BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].Pbs[0] =
                                                    maxP12;

                                                maxP22 =
                                                    p2[
                                                        BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                            .JobsIndice[0]];

                                                foreach (
                                                    int t in
                                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                        .JobsIndice)
                                                    if (p2[t] > maxP22)
                                                        maxP22 = p2[t];

                                                BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].Pbs[1] =
                                                    maxP22;
                                            }
                                        }
                                        #endregion

                                    }

                                }



                                #endregion

                                break;
                        }
                        for (int j = 0; j < selectedJobsAfterOPs.Length; j++)
                            PiJ[j] = !selectedJobsAfterOPs[j] ? 1 : 0;

                        modelAfterOPs.IgnoranceBinary = PiJ;

                        modelAfterOPs.NumberOfNonEmptyBatches = nonEmptyBatchesAfterOPs.Count;

                        t1 = new double[modelAfterOPs.NumberOfMachinesInStep1];

                        t2 = new double[modelAfterOPs.NumberOfMachinesInStep2];

                        Tj = new double[N];

                        solAfterOPs = Algorithm1(5, nonEmptyBatches, t1, t2, Tj, d, t_now);

                        modelAfterOPs.DelayOfJobs = solAfterOPs.Tj;

                        double cost = CostFunction(modelAfterOPs);

                        if (cost < ant[k].Cost)
                        {
                            ant[k].Cost = cost;

                            ant[k].sol = solAfterOPs;

                            ant[k].SelectedJobs = selectedJobsAfterOPs;

                            ant[k].SelectedJobsForEmptyBatches = selectedJobsForEmptyBatches;

                            model = modelAfterOPs;

                            nonEmptyBatches = nonEmptyBatchesAfterOPs;

                            sol = solAfterOPs;

                            selectedJobs = selectedJobsAfterOPs;

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

            string pathToExcelFile = "D:\\125.xls";
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
