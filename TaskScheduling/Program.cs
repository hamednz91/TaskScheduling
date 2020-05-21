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
        static void Run_ACO(string inputFilePath, long elapsedTime)
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

            Model model = new Model();

            model.CreateModel(inputFilePath);

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

            int maxIteration = 10000;

            int numberOfAnts = 50;

            int Q = 10;

            double alpha = 1;

            double beta1 = 1;
            double beta2 = 1;
            double beta3 = 1;
            double beta4 = 1;

            double rho = 0.05;
            double rhoL = 0.05;

            #endregion

            #region ACO Initialization

            double[] eta1J = new double[N];

            double[] eta2J = new double[N];

            double[] eta3J = new double[N];

            double[] eta4J = new double[N];

            double[,] phiJX = new double[N, N];

            double[,] phiJX0 = new double[N, N];

            int[,] mJX = new int[N, N];

            bool[,] mJXFlag = new bool[N, N];

            double[] R = new double[N];

            double[] tauJB = new double[N];

            double[] tauJ = new double[N];

            double[] tauJ0 = new double[N];

            for (int i = 0; i < N; i++)
            {
                eta1J[i] = wTj[i];

                eta2J[i] = (double)1 / (d[i] - (double)(p1[i] + p2[i]));

                eta3J[i] = 1;

                eta4J[i] = sigmaJ[i];

                tauJB[i] = 1;

                tauJ[i] = 1;

                tauJ0[i] = 1;

                for (int j = 0; j < N; j++)
                {
                    phiJX[i, j] = 0.1;

                    phiJX0[i, j] = 0.1;
                }
            }

            Ant[] ant = new Ant[numberOfAnts];

            Ant bestAnt = new Ant();

            bestAnt.Tour = new List<Batch>();

            bestAnt.Cost = double.MaxValue;

            bestAnt.SelectedJobs = new bool[N];

            bestAnt.SelectedJobsForEmptyBatches = new bool[N];

            bestAnt.R = new double[N];

            bestAnt.mJXFlag = new bool[N, N];

            Sol sol = new Sol();

            sol.BatchesAllocatedToMachines = new List<Batch>();

            sol.TimeofMachinesStep1 = new double[t1.Length];

            sol.TimeofMachinesStep2 = new double[t2.Length];

            #endregion

            #region ACO Main Loop

            for (int it = 0; it < maxIteration && sw.ElapsedMilliseconds < elapsedTime; it++)
            {
                Ant bestAntPerIteration = new Ant();

                bestAntPerIteration.Tour = new List<Batch>();

                bestAntPerIteration.Cost = double.MaxValue;

                bestAntPerIteration.SelectedJobs = new bool[N];

                bestAntPerIteration.SelectedJobsForEmptyBatches = new bool[N];

                bestAntPerIteration.R = new double[N];

                bestAntPerIteration.mJXFlag = new bool[N, N];

                R = new double[N];

                // Ants Movement
                for (int k = 0; k < numberOfAnts; k++)
                {
                    #region Parameters Initialization
                    int t_now = 0;

                    sol = new Sol();

                    ant[k] = new Ant();

                    ant[k].Tour = new List<Batch>();

                    ant[k].sol = new Sol();

                    ant[k].SelectedJobs = new bool[N];

                    ant[k].SelectedJobsForEmptyBatches = new bool[N];

                    ant[k].R = new double[N];

                    ant[k].mJXFlag = new bool[N, N];

                    t1 = new double[model.NumberOfMachinesInStep1];

                    t2 = new double[model.NumberOfMachinesInStep2];

                    selectedJobs = new bool[N];

                    selectedJobsForEmptyBatches = new bool[N];

                    mJXFlag = new bool[N, N];
                    #endregion

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

                        batches[i].Pbs = new double[2]; //batch processing time in step 1 & 2

                        batches[i].Family = -1;

                        batches[i].machineNumber = new int[] { -1, -1 };

                        batches[i].batchIndex = -1;

                        batches[i].BatchID = -1;

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

                        virtualBatches[i].Pbs = new double[2]; //batch processing time in step 1 & 2

                        virtualBatches[i].Family = -1;

                        virtualBatches[i].machineNumber = new int[] { -1, -1 };

                        virtualBatches[i].batchIndex = -1;

                        virtualBatches[i].BatchID = -1;

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
                            if (!Helper.CheckExistInBatch(batches[b], i) &&
                                Helper.CheckInSameFamily(hjf, rand, i,
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

                            double s = 0;

                            counter++;

                            for (int j = 0; j < batches[b].BatchCandidateList.Count; j++)
                            {
                                eta3J[batches[b].BatchCandidateList[j]] = (double)1 /
                                                                          (double)
                                                                          (Math.Abs((kMax - batches[b].SizeOfJobs.Sum()) -
                                                                           Sj[batches[b].BatchCandidateList[j]]) + 1);
                                double sumOfPhiJXs = 0;

                                for (int l = 0; l < batches[b].JobsIndice.Count; l++)
                                    sumOfPhiJXs += phiJX[batches[b].BatchCandidateList[j], batches[b].JobsIndice[l]];

                                tauJB[batches[b].BatchCandidateList[j]] = sumOfPhiJXs / batches[b].JobsIndice.Count;

                            }

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

                                s += (Math.Pow((batches[b].TauJBCandidates[l] + batches[b].TauJCandidates[l]), alpha) *
                                      Math.Pow(batches[b].Eta1jCandidates[l], beta1) *
                                      Math.Pow(batches[b].Eta2jCandidates[l], beta2) *
                                      Math.Pow(batches[b].Eta3jCandidates[l], beta3) *
                                      Math.Pow(batches[b].Eta4jCandidates[l], beta4));
                            }
                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                            {
                                batches[b].JobCandidateSelectionProbability[l] =
                                    (Math.Pow((batches[b].TauJBCandidates[l] + batches[b].TauJCandidates[l]), alpha) *
                                     Math.Pow(batches[b].Eta1jCandidates[l], beta1) *
                                     Math.Pow(batches[b].Eta2jCandidates[l], beta2) *
                                     Math.Pow(batches[b].Eta3jCandidates[l], beta3) *
                                     Math.Pow(batches[b].Eta4jCandidates[l], beta4)
                                    ) / s;
                                sumOfSum += batches[b].JobCandidateSelectionProbability[l];
                            }
                            for (int l = 0; l < batches[b].BatchCandidateList.Count; l++)
                                batches[b].JobCandidateSelectionProbability[l] =
                                    batches[b].JobCandidateSelectionProbability[l] /
                                    sumOfSum;


                            int ind = Helper.RouletteWheelSelection(batches[b].JobCandidateSelectionProbability.ToArray());

                            if (Sj[batches[b].BatchCandidateList[ind]] + batches[b].SizeOfJobs.Sum() > kMax)
                                continue;

                            if (!selectedJobs[batches[b].BatchCandidateList[ind]])
                            {
                                batches[b].JobsIndice.Add(batches[b].BatchCandidateList[ind]);

                                batches[b].SizeOfJobs.Add(Sj[batches[b].BatchCandidateList[ind]]);

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

                            batches[b].Pbs = new double[2];

                            batches[b].Family = -1;

                            batches[b].batchIndex = -1;

                            batches[b].BatchID = -1;

                            batches[b].machineNumber = new int[] { -1, -1 };

                        }

                        #endregion

                        #region Add Pbs

                        batches[b].Pbs[0] = maxP1j;

                        batches[b].Pbs[1] = maxP2j;

                        #endregion


                    }

                    #endregion

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

                            batches[b].batchIndex = count;

                            batches[b].BatchID = count++;

                            nonEmptyBatches.Add(batches[b]);
                        }

                    #endregion

                    for (int j = 0; j < selectedJobs.Length; j++)
                        PiJ[j] = !selectedJobs[j] ? 1 : 0;

                    model.IgnoranceBinary = PiJ;

                    model.NumberOfNonEmptyBatches = nonEmptyBatches.Count;

                    sol = Helper.Algorithm1(4, nonEmptyBatches, t1, t2, Tj, d, t_now);

                    model.DelayOfJobs = sol.Tj;

                    ant[k].Cost = model.CostFunction();

                    ant[k].sol = sol;

                    ant[k].SelectedJobs = selectedJobs;

                    ant[k].SelectedJobsForEmptyBatches = selectedJobsForEmptyBatches;

                    if (ant[k].Cost < bestAntPerIteration.Cost)
                    {
                        bestAntPerIteration = ant[k];
                    }

                    nonEmptyBatches = sol.BatchesAllocatedToMachines;

                    for (int x = 0; x < nonEmptyBatches.Count; x++)
                        nonEmptyBatches[x].batchIndex = x;

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
                        virtualBatches[j].Pbs = new double[2]; //batch processing time in step 1 & 2
                                                               // virtualBatches[j].SumOfJobSizes = 0;
                        virtualBatches[j].Family = -1;
                        virtualBatches[j].machineNumber = new int[] { -1, -1 };
                        virtualBatches[j].batchIndex = -1;
                        virtualBatches[j].BatchID = -1;
                    }


                    #endregion

                    for (int j = 0; j < selectedJobs.Length; j++)
                    {
                        if (!selectedJobs[j])
                        {
                            for (int i = 0; i < model.NumberOfFamilies; i++)
                            {
                                if (hjf[i, j] == 1)
                                {
                                    virtualBatches[i].JobsIndice.Add(j);

                                    virtualBatches[i].SizeOfJobs.Add(Sj[j]);

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


                    Model modelAfterOPs = new Model(model);

                    bool[] selectedJobsAfterOPs = new bool[selectedJobs.Length];

                    for (int i = 0; i < selectedJobsAfterOPs.Length; i++)
                        selectedJobsAfterOPs[i] = selectedJobs[i];

                    Sol solAfterOPs = new Sol(sol);

                    #region Operators

                    for (int i = 0; i < A; i++)
                    {
                        List<Batch> BatchesGreaterThanKmin = new List<Batch>();

                        foreach (var item in nonEmptyBatchesAfterOPs)
                            if (item.JobsIndice.Count > kMin)
                                BatchesGreaterThanKmin.Add(item);

                        int opSelector = Helper.RouletteWheelSelection(new[] { .1, .1, .0, .1, .1, .1, .1, .1, .0, .1, .2 });

                        //opSelector = 10;

                        switch (opSelector)
                        {
                            case 0:

                                #region OP0 Drop

                                bool stopflagOP0 = BatchesGreaterThanKmin.Count == 0;

                                if (stopflagOP0) break;

                                int selectedBatchIndexOP0 = r.Next(BatchesGreaterThanKmin.Count);

                                //selectedBatchIndex = BatchesGreaterThanKmin[selectedBatchIndex].batchIndex;
                                // real index in nonEmptyBatchesAfterOPs using batch index field

                                int selectedBatchFamilyOP0 = BatchesGreaterThanKmin[selectedBatchIndexOP0].Family;

                                int selectedBatchLengthOP0 = BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice.Count;

                                int nOP0 = r.Next(1, Math.Max(1, selectedBatchLengthOP0 - kMin));

                                for (int j = 0; j < nOP0; j++)
                                {

                                    int jobIndexOP0 = r.Next(selectedBatchLengthOP0);

                                    if (selectedBatchLengthOP0 > 0)
                                    {
                                        #region Add Selected Jobs to Virtual Batch

                                        virtualBatches[selectedBatchFamilyOP0].JobsIndice.Add(
                                            BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]);

                                        virtualBatches[selectedBatchFamilyOP0].SizeOfJobs.Add(
                                            Sj[BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]]);

                                        virtualBatches[selectedBatchFamilyOP0].Family = selectedBatchFamilyOP0;

                                        virtualBatches[selectedBatchFamilyOP0].UrgentMetric.Add(1);

                                        virtualBatches[selectedBatchFamilyOP0].DueTime.Add(
                                            d[BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]]);

                                        #endregion

                                        #region Remove Selected Jobs from NonEmptyBatches and Update Batch (size of jobs, Urgent metric, pbs if needed)

                                        bool flag1 =
                                            p1[BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]] >=
                                            BatchesGreaterThanKmin[selectedBatchIndexOP0].Pbs[0];

                                        bool flag2 =
                                            p2[BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]] >=
                                            BatchesGreaterThanKmin[selectedBatchIndexOP0].Pbs[1];

                                        selectedJobsAfterOPs[
                                            BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice[jobIndexOP0]] = false;

                                        BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice.RemoveAt(jobIndexOP0);

                                        BatchesGreaterThanKmin[selectedBatchIndexOP0].SizeOfJobs.RemoveAt(jobIndexOP0);

                                        BatchesGreaterThanKmin[selectedBatchIndexOP0].UrgentMetric.RemoveAt(jobIndexOP0);

                                        BatchesGreaterThanKmin[selectedBatchIndexOP0].DueTime.RemoveAt(jobIndexOP0);

                                        selectedBatchLengthOP0 =
                                            BatchesGreaterThanKmin[selectedBatchIndexOP0].JobsIndice.Count;

                                        if (selectedBatchLengthOP0 <= 0) continue;

                                        if (flag1)
                                        {
                                            double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].JobsIndice)
                                                if (p1[t] > maxP1)
                                                    maxP1 = p1[t];

                                            nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].Pbs[0] = maxP1;
                                        }
                                        if (flag2)
                                        {
                                            double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].JobsIndice[0]];

                                            foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].JobsIndice)
                                                if (p2[t] > maxP2)
                                                    maxP2 = p2[t];

                                            nonEmptyBatchesAfterOPs[selectedBatchIndexOP0].Pbs[1] = maxP2;
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

                                bool[] selectbatchesOP1 = new bool[nonEmptyBatchesAfterOPs.Count];

                                int selectedBatchIndex1OP1 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                selectbatchesOP1[selectedBatchIndex1OP1] = true;

                                int selectedFamilyOP1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].Family;

                                int selectedBatchIndex2OP1 = -1;

                                do
                                {
                                    selectedBatchIndex2OP1 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectbatchesOP1[selectedBatchIndex2OP1])
                                        selectbatchesOP1[selectedBatchIndex2OP1] = true;

                                } while ((selectedBatchIndex1OP1 == selectedBatchIndex2OP1 ||
                                          selectedFamilyOP1 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].Family) &&
                                         selectbatchesOP1.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (selectedFamilyOP1 != nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].Family) break;

                                int selectedBatchLength1OP1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice.Count;
                                int selectedBatchLength2OP1 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice.Count;

                                //
                                //can number 2 be chosen in the random selection below?
                                //
                                int nOP1 = r.Next(1, Math.Min(selectedBatchLength1OP1, selectedBatchLength2OP1));

                                for (int j = 0; j < nOP1; j++)
                                {
                                    int batchSize1OP1 = 0, batchSize2OP1 = 0;

                                    int jobIndex1OP1 = r.Next(selectedBatchLength1OP1);

                                    int jobIndex2OP1 = -1;

                                    // int numberOfSelectedBatches2 = 0;

                                    bool[] selectjobIndex2OP1 = new bool[selectedBatchLength2OP1];

                                    do jobIndex2OP1 = r.Next(selectedBatchLength2OP1); while (selectjobIndex2OP1[jobIndex2OP1] && selectjobIndex2OP1.Count(item => item) < selectedBatchLength2OP1);

                                    int a1OP1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs.Sum() -
                                             nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs[jobIndex1OP1];

                                    batchSize1OP1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs[jobIndex2OP1] +
                                        a1OP1;

                                    int a2OP1 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs.Sum() -
                                             nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs[jobIndex2OP1];

                                    batchSize2OP1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs[jobIndex1OP1] +
                                        a2OP1;

                                    if (!selectjobIndex2OP1[jobIndex2OP1])
                                        selectjobIndex2OP1[jobIndex2OP1] = true;

                                    if (batchSize1OP1 > kMax || batchSize2OP1 > kMax) break;

                                    int job1 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice[jobIndex1OP1];
                                    int job2 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice[jobIndex2OP1];

                                    int jobSize1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs[jobIndex1OP1];
                                    int jobSize2 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs[jobIndex2OP1];

                                    double urgentMetric1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].UrgentMetric[jobIndex1OP1];
                                    double urgentMetric2 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].UrgentMetric[jobIndex2OP1];

                                    double dueTime1 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].DueTime[jobIndex1OP1];
                                    double dueTime2 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].DueTime[jobIndex2OP1];

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice.Add(job2);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice.Add(job1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs.Add(jobSize2);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs.Add(jobSize1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].UrgentMetric.Add(urgentMetric2);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].UrgentMetric.Add(urgentMetric1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].DueTime.Add(dueTime2);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].DueTime.Add(dueTime1);


                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice.RemoveAt(jobIndex1OP1);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice.RemoveAt(jobIndex2OP1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].SizeOfJobs.RemoveAt(jobIndex1OP1);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].SizeOfJobs.RemoveAt(jobIndex2OP1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].UrgentMetric.RemoveAt(jobIndex1OP1);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].UrgentMetric.RemoveAt(jobIndex2OP1);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].DueTime.RemoveAt(jobIndex1OP1);
                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].DueTime.RemoveAt(jobIndex2OP1);


                                    double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice[0]];

                                    foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice)
                                        if (p1[t] > maxP1)
                                            maxP1 = p1[t];

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].Pbs[0] = maxP1;

                                    double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice[0]];

                                    foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].JobsIndice)
                                        if (p2[t] > maxP2)
                                            maxP2 = p2[t];

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP1].Pbs[1] = maxP2;


                                    double maxP12 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice[0]];

                                    foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice)
                                        if (p1[t] > maxP12)
                                            maxP12 = p1[t];

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].Pbs[0] = maxP12;

                                    double maxP22 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice[0]];

                                    foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].JobsIndice)
                                        if (p2[t] > maxP22)
                                            maxP22 = p2[t];

                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP1].Pbs[1] = maxP22;
                                }


                                #endregion

                                break;
                            case 2:

                                #region OP2 Create new batch from VirtualBatches

                                bool stopFlagOP2 = virtualBatches.All(vb => vb.JobsIndice.Count < kMin);

                                if (stopFlagOP2) break;

                                #region New Batch Instance

                                Batch newBatchOP2 = new Batch();

                                #endregion

                                int virtualBatchIndexOP2 = -1;

                                List<Batch> VirtualBatchesMoreThanKmin = new List<Batch>();

                                foreach (var batch in virtualBatches)
                                    if (batch.JobsIndice.Count >= kMin)
                                        VirtualBatchesMoreThanKmin.Add(batch);

                                bool[] selectedFamilyVirtualBatch = new bool[VirtualBatchesMoreThanKmin.Count];

                                #region New Batch Init/Reset

                                newBatchOP2.JobsIndice = new List<int>();
                                newBatchOP2.SizeOfJobs = new List<int>();
                                newBatchOP2.UrgentMetric = new List<double>();
                                newBatchOP2.DueTime = new List<double>();
                                newBatchOP2.Pbs = new double[2]; //batch processing time in step 1 & 2
                                newBatchOP2.Family = -1;
                                newBatchOP2.machineNumber = new int[] { -1, -1 };
                                newBatchOP2.batchIndex = -1;
                                newBatchOP2.BatchID = -1;

                                #endregion

                                do
                                {
                                    virtualBatchIndexOP2 = r.Next(VirtualBatchesMoreThanKmin.Count);

                                } while (selectedFamilyVirtualBatch[virtualBatchIndexOP2] &&
                                         selectedFamilyVirtualBatch.Count(item => item) <
                                         VirtualBatchesMoreThanKmin.Count);

                                if (!selectedFamilyVirtualBatch[virtualBatchIndexOP2])
                                    selectedFamilyVirtualBatch[virtualBatchIndexOP2] = true;

                                int selectedVirBatchLength =
                                    VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].JobsIndice.Count;

                                double avgSjOfSelectedVirBatch =
                                    VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].SizeOfJobs.Average();

                                int nOP2 = r.Next(kMin, Convert.ToInt32(Math.Ceiling((kMax) / avgSjOfSelectedVirBatch)));

                                int jobIndexOP2 = -1;

                                int jobOP2 = -1;

                                bool[] selectedJobFromSelectedVirtualBatchOP2 = new bool[selectedVirBatchLength];

                                for (int j = 0; j < Math.Min(nOP2, selectedVirBatchLength); j++)
                                {
                                    do
                                    {
                                        jobIndexOP2 = r.Next(selectedVirBatchLength);

                                        jobOP2 = VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].JobsIndice[jobIndexOP2];

                                    } while (selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] &&
                                             selectedJobFromSelectedVirtualBatchOP2.Count(item => item) <
                                             selectedVirBatchLength);

                                    if (!selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] &&
                                        Sj[jobOP2] + newBatchOP2.SizeOfJobs.Sum() < kMax)
                                    {
                                        newBatchOP2.JobsIndice.Add(jobOP2);
                                        newBatchOP2.SizeOfJobs.Add(Sj[jobOP2]);
                                        newBatchOP2.UrgentMetric.Add(1);
                                        newBatchOP2.DueTime.Add(d[jobOP2]);
                                        newBatchOP2.Family = virtualBatchIndexOP2;

                                        if (!selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2])
                                            selectedJobFromSelectedVirtualBatchOP2[jobIndexOP2] = true;

                                        double maxP1 = p1[jobOP2];

                                        foreach (int t in newBatchOP2.JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        newBatchOP2.Pbs[0] = maxP1;

                                        double maxP2 = p2[jobOP2];

                                        foreach (int t in newBatchOP2.JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        newBatchOP2.Pbs[1] = maxP2;

                                        VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].JobsIndice.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].SizeOfJobs.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].UrgentMetric.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].DueTime.RemoveAt(jobIndexOP2);
                                        VirtualBatchesMoreThanKmin.ToArray()[virtualBatchIndexOP2].Family = -1;

                                        selectedVirBatchLength =
                                            VirtualBatchesMoreThanKmin[virtualBatchIndexOP2].JobsIndice.Count;
                                    }
                                }

                                List<Batch> nonEmptyBatchesAverageDjToWjOP2 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                    nonEmptyBatchesAverageDjToWjOP2.Add(item);

                                if (newBatchOP2.JobsIndice.Count >= kMin)
                                {
                                    double average = 0;
                                    foreach (var item in newBatchOP2.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        selectedJobsAfterOPs[item] = true;

                                    }
                                    average = average / (double)newBatchOP2.JobsIndice.Count;

                                    newBatchOP2.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    int l = 0;

                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP2.Count; j++)
                                    {
                                        newBatchOP2.BatchID = nonEmptyBatchesAverageDjToWjOP2.Count;
                                        if (average >
                                            nonEmptyBatchesAverageDjToWjOP2[j].AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            newBatchOP2.batchIndex = j + 1;

                                            nonEmptyBatchesAverageDjToWjOP2.Insert(j + 1, newBatchOP2);

                                            isSet = true;

                                            l = j + 2;

                                            break;
                                        }
                                    }

                                    if (isSet)
                                        for (; l < nonEmptyBatchesAverageDjToWjOP2.Count; l++)
                                            nonEmptyBatchesAverageDjToWjOP2[l].batchIndex++;
                                    else
                                        nonEmptyBatchesAverageDjToWjOP2.Insert(0, newBatchOP2);


                                }

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


                                int nOP3 = r.Next(1,
                                        Math.Min(selectedBatchLengthOP3,
                                            virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count));

                                for (int j = 0; j < nOP3; j++)
                                {
                                    int batchSize1OP3 = 0;

                                    int jobIndex1OP3 = r.Next(selectedBatchLengthOP3);

                                    int jobIndex2OP3 = -1;

                                    bool[] selectJobIndex2OP3 =
                                        new bool[virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count];

                                    int numberofSelectedJobfromBatchOP3 = 0;

                                    do
                                    {
                                        //check if jobindex2 is already chosen
                                        numberofSelectedJobfromBatchOP3 = selectJobIndex2OP3.Count(item => item);

                                        jobIndex2OP3 = r.Next(virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count);

                                        int a1OP3 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.Sum() -
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs[jobIndex1OP3];

                                        batchSize1OP3 = virtualBatches[selectedBatchFamilyOP3].SizeOfJobs[jobIndex2OP3] +
                                                     a1OP3;

                                        if (!selectJobIndex2OP3[jobIndex2OP3])
                                            selectJobIndex2OP3[jobIndex2OP3] = true;

                                    } while (batchSize1OP3 > kMax && numberofSelectedJobfromBatchOP3 <
                                             virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count);

                                    if (numberofSelectedJobfromBatchOP3 >=
                                        virtualBatches[selectedBatchFamilyOP3].JobsIndice.Count)
                                        break;

                                    int job1OP3 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[jobIndex1OP3];
                                    int job2OP3 = virtualBatches[selectedBatchFamilyOP3].JobsIndice[jobIndex2OP3];

                                    int jobSize1OP3 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs[jobIndex1OP3];
                                    int jobSize2OP3 = virtualBatches[selectedBatchFamilyOP3].SizeOfJobs[jobIndex2OP3];


                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice.Add(job2OP3);
                                    virtualBatches[selectedBatchFamilyOP3].JobsIndice.Add(job1OP3);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.Add(jobSize2OP3);
                                    virtualBatches[selectedBatchFamilyOP3].SizeOfJobs.Add(jobSize1OP3);

                                    bool flag1OP3 = p1[job1OP3] >=
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[0];

                                    bool flag2OP3 = p2[job1OP3] >=
                                                 nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[1];


                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice.RemoveAt(jobIndex1OP3);
                                    virtualBatches[selectedBatchFamilyOP3].JobsIndice.RemoveAt(jobIndex2OP3);

                                    nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].SizeOfJobs.RemoveAt(jobIndex1OP3);
                                    virtualBatches[selectedBatchFamilyOP3].SizeOfJobs.RemoveAt(jobIndex2OP3);

                                    if (flag1OP3)
                                    {
                                        double maxP1 = p1[nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[0] = maxP1;
                                    }

                                    if (flag2OP3)
                                    {
                                        double maxP2 = p2[nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice[0]];

                                        foreach (int t in nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndexOP3].Pbs[1] = maxP2;
                                    }

                                    selectedJobsAfterOPs[job1OP3] = false;

                                    selectedJobsAfterOPs[job2OP3] = true;

                                }


                                #endregion

                                break;
                            case 4:

                                #region OP4 Add

                                bool stopFlagOP4 = virtualBatches.All(vb => vb.JobsIndice.Count == 0) ||
                                                   nonEmptyBatchesAfterOPs.Count == 0;

                                if (stopFlagOP4) break;

                                int selectedVirtualBatchIndexOP4 = -1;

                                int selectVirBatchCounterOP4 = 0;

                                bool[] selectedFamilyVirtualBatchOP4 = new bool[model.NumberOfFamilies];
                                do
                                {
                                    selectedVirtualBatchIndexOP4 = r.Next(model.NumberOfFamilies);

                                    if (!selectedFamilyVirtualBatchOP4[selectedVirtualBatchIndexOP4])
                                    {
                                        selectedFamilyVirtualBatchOP4[selectedVirtualBatchIndexOP4] = true;
                                        selectVirBatchCounterOP4++;
                                    }

                                } while (virtualBatches[selectedVirtualBatchIndexOP4].JobsIndice.Count == 0 &&
                                         selectVirBatchCounterOP4 < virtualBatches.Length);

                                int SelectedBatchIndexOP4 = -1;

                                bool[] selectedBatchOP4 = new bool[nonEmptyBatchesAfterOPs.Count];

                                do
                                {
                                    SelectedBatchIndexOP4 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                    if (!selectedBatchOP4[SelectedBatchIndexOP4])
                                        selectedBatchOP4[SelectedBatchIndexOP4] = true;

                                } while (nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].Family !=
                                         selectedVirtualBatchIndexOP4 &&
                                         selectedBatchOP4.Count(item => item) < nonEmptyBatchesAfterOPs.Count);

                                if (nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].Family != selectedVirtualBatchIndexOP4)
                                    break;

                                if (nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].JobsIndice.Count == 0)
                                {
                                    break;
                                }
                                double avgSjOfSelectedVirBatchOP4 =
                                    nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].SizeOfJobs.Average();

                                int nOP4 = r.Next(1, Convert.ToInt32(Math.Ceiling(kMax / avgSjOfSelectedVirBatchOP4)));

                                for (int j = 0;
                                    j < Math.Min(virtualBatches[selectedVirtualBatchIndexOP4].JobsIndice.Count, nOP4);
                                    j++)
                                {
                                    int jobIndexOP4 = r.Next(virtualBatches[selectedVirtualBatchIndexOP4].JobsIndice.Count);
                                    int jobOP4 = virtualBatches[selectedVirtualBatchIndexOP4].JobsIndice[jobIndexOP4];
                                    int jobSizeOP4 = virtualBatches[selectedVirtualBatchIndexOP4].SizeOfJobs[jobIndexOP4];

                                    if (jobSizeOP4 +
                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].SizeOfJobs.Sum() <= kMax)
                                    {
                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].JobsIndice.Add(jobOP4);
                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].SizeOfJobs.Add(Sj[jobOP4]);
                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].UrgentMetric.Add(1);
                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].DueTime.Add(d[jobOP4]);


                                        double maxP1 = p1[jobOP4];

                                        foreach (int t in nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].JobsIndice)
                                            if (p1[t] > maxP1)
                                                maxP1 = p1[t];

                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].Pbs[0] = maxP1;

                                        double maxP2 = p2[jobOP4];

                                        foreach (int t in nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].JobsIndice)
                                            if (p2[t] > maxP2)
                                                maxP2 = p2[t];

                                        nonEmptyBatchesAfterOPs[SelectedBatchIndexOP4].Pbs[1] = maxP2;

                                        selectedJobsAfterOPs[jobOP4] = true;

                                        virtualBatches[selectedVirtualBatchIndexOP4].JobsIndice.RemoveAt(jobIndexOP4);
                                        virtualBatches[selectedVirtualBatchIndexOP4].SizeOfJobs.RemoveAt(jobIndexOP4);

                                    }

                                }



                                #endregion Add 

                                break;
                            case 5:

                                #region OP5 Remove Random Batch

                                bool stopFlagOP5 = nonEmptyBatchesAfterOPs.Count == 0;

                                if (stopFlagOP5) break;

                                int selectedBatchIndexOP5 = r.Next(nonEmptyBatchesAfterOPs.Count);

                                int selectedBatchFamilyOP5 = nonEmptyBatchesAfterOPs[selectedBatchIndexOP5].Family;

                                foreach (var item in nonEmptyBatchesAfterOPs[selectedBatchIndexOP5].JobsIndice)
                                {
                                    #region Add Selected Jobs to Virtual Batch

                                    virtualBatches[selectedBatchFamilyOP5].JobsIndice.Add(item);

                                    virtualBatches[selectedBatchFamilyOP5].SizeOfJobs.Add(Sj[item]);

                                    virtualBatches[selectedBatchFamilyOP5].Family = selectedBatchFamilyOP5;

                                    virtualBatches[selectedBatchFamilyOP5].UrgentMetric.Add(1);

                                    virtualBatches[selectedBatchFamilyOP5].DueTime.Add(d[item]);

                                    #endregion

                                    selectedJobsAfterOPs[item] = false;
                                }

                                nonEmptyBatchesAfterOPs.RemoveAt(selectedBatchIndexOP5);

                                Batch NewNonEmptyBatchOP5 = new Batch();

                                List<Batch> myNewNonEmptyBatchesOP5 = new List<Batch>();

                                for (int j = 0; j < nonEmptyBatchesAfterOPs.Count; j++)
                                {
                                    NewNonEmptyBatchOP5 = nonEmptyBatchesAfterOPs[j];

                                    NewNonEmptyBatchOP5.batchIndex = j;

                                    myNewNonEmptyBatchesOP5.Add(NewNonEmptyBatchOP5);
                                }


                                nonEmptyBatchesAfterOPs = myNewNonEmptyBatchesOP5;

                                #endregion Remove Bandom Batch

                                break;
                            case 6:

                                #region OP6 Transform from one nonEmpty to another

                                bool stopflagOP6 = nonEmptyBatchesAfterOPs.Count < 2 || BatchesGreaterThanKmin.Count == 0;

                                if (stopflagOP6) break;

                                bool[] selectbatchesOP6 = new bool[nonEmptyBatchesAfterOPs.Count];

                                int selectedBatchIndex1OP6 = r.Next(BatchesGreaterThanKmin.Count);

                                if (!selectbatchesOP6[selectedBatchIndex1OP6])
                                    selectbatchesOP6[selectedBatchIndex1OP6] = true;

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

                                double minMetricOP6 =
                                    (kMax - nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Sum()) /
                                    BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs.Average();

                                // minMetric wouldn`t be integer, and the Minimum would be double as well. is it correct to convert to integer?
                                int c = r.Next(1, (int)Math.Max(1, Math.Min(minMetricOP6, (selectedBatchLength1OP6 - kMin))));

                                for (int j = 0; j < c; j++)
                                {
                                    selectedBatchLength1OP6 =
                                    BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice.Count;

                                    int jobIndex1OP6 = r.Next(selectedBatchLength1OP6);

                                    if ((nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Sum() +
                                         BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1OP6] <= kMax))
                                    {
                                        int job1OP6 = BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice[jobIndex1OP6];

                                        int jobSize1OP6 =
                                            BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs[jobIndex1OP6];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].JobsIndice.Add(job1OP6);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP6].SizeOfJobs.Add(jobSize1OP6);

                                        BatchesGreaterThanKmin[selectedBatchIndex1OP6].JobsIndice.RemoveAt(jobIndex1OP6);

                                        BatchesGreaterThanKmin[selectedBatchIndex1OP6].SizeOfJobs.RemoveAt(jobIndex1OP6);

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

                                List<Batch> batchesForExchangeBatchIndiceOP7 = nonEmptyBatchesAfterOPs;

                                int selectedBatchIndex1OP7 = r.Next(batchesForExchangeBatchIndiceOP7.Count);

                                int selectedBatchIndex2OP7 = -1;

                                do
                                    selectedBatchIndex2OP7 = r.Next(batchesForExchangeBatchIndiceOP7.Count); while (
                                    selectedBatchIndex1OP7 == selectedBatchIndex2OP7);

                                int temp_BatchIndexOP7 = batchesForExchangeBatchIndiceOP7[selectedBatchIndex1OP7].batchIndex;

                                batchesForExchangeBatchIndiceOP7[selectedBatchIndex1OP7].batchIndex =
                                    batchesForExchangeBatchIndiceOP7[selectedBatchIndex2OP7].batchIndex;

                                batchesForExchangeBatchIndiceOP7[selectedBatchIndex2OP7].batchIndex = temp_BatchIndexOP7;

                                nonEmptyBatchesAfterOPs = batchesForExchangeBatchIndiceOP7;

                                #endregion Change nonEmpty batchindice

                                break;
                            case 8:

                                #region OP8 Create New Batch from nonEmptyBatchesAfterOPs

                                bool stopFlagOP8 = nonEmptyBatchesAfterOPs.All(item => item.JobsIndice.Count <= 2 * kMin);

                                if (stopFlagOP8) break;

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

                                List<Batch> nonEmptyBatchesAverageDjToWjOP8 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                    nonEmptyBatchesAverageDjToWjOP8.Add(item);

                                if (newBatchOP8.JobsIndice.Count >= kMin)
                                {

                                    double average = 0;
                                    foreach (var item in newBatchOP8.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                    }
                                    average = average / (double)newBatchOP8.JobsIndice.Count;

                                    newBatchOP8.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    bool isSet = false;

                                    int l = 0;

                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP8.Count; j++)
                                    {
                                        newBatchOP8.BatchID = nonEmptyBatchesAverageDjToWjOP8.Count;
                                        if (average >
                                            nonEmptyBatchesAverageDjToWjOP8[j]
                                                .AverageDueTimeofJobToDelayImportanceFactor)
                                        {
                                            newBatchOP8.batchIndex = j + 1;

                                            nonEmptyBatchesAverageDjToWjOP8.Insert(j + 1, newBatchOP8);

                                            isSet = true;

                                            l = j + 2;

                                            break;

                                        }
                                    }



                                    if (isSet)
                                        for (; l < nonEmptyBatchesAverageDjToWjOP8.Count; l++)
                                        {
                                            nonEmptyBatchesAverageDjToWjOP8[l].batchIndex++;
                                        }
                                    else
                                        nonEmptyBatchesAverageDjToWjOP8.Insert(0, newBatchOP8);


                                }

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

                                int nOP9 = r.Next(1, Math.Min(selectedBatchLength1OP9, selectedBatchLength2OP9));

                                for (int j = 0; j < nOP9; j++)
                                {
                                    int batchSize1OP9 = 0, batchSize2OP9 = 0;

                                    int jobIndex1OP9 =
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.IndexOf(
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.Max());

                                    int jobIndex2OP9 = -1;

                                    int counterOP9 = 0;

                                    do
                                    {

                                        jobIndex2OP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.IndexOf(
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.OrderBy(dt => dt)
                                                .Skip(counterOP9++)
                                                .First());

                                        int a1OP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.Sum() -
                                                    nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[
                                                        jobIndex1OP9];

                                        batchSize1OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[jobIndex2OP9] +
                                            a1OP9;

                                        int a2OP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.Sum() -
                                                    nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[
                                                        jobIndex2OP9];

                                        batchSize2OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[jobIndex1OP9] +
                                            a2OP9;


                                    } while (counterOP9 < selectedBatchLength2OP9 &&
                                             (batchSize1OP9 > kMax || batchSize2OP9 > kMax));

                                    if (counterOP9 == selectedBatchLength2OP9 || selectedBatchesOP9.Count(item => item) >= nonEmptyBatchesAfterOPs.Count - 1)
                                        break;

                                    if (selectedBatchesOP9.Count(item => item) < numberOfSameFamilyBatches2OP9 - 1)
                                    {

                                        int job1OP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice[jobIndex1OP9];
                                        int job2OP9 = nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice[jobIndex2OP9];

                                        int jobSize1OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs[jobIndex1OP9];
                                        int jobSize2OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs[jobIndex2OP9];

                                        double urgentMetric1OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric[jobIndex1OP9];
                                        double urgentMetric2OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric[jobIndex2OP9];

                                        double dueTime1OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime[jobIndex1OP9];
                                        double dueTime2OP9 =
                                            nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime[jobIndex2OP9];

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice.Add(job2OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.Add(job1OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.Add(jobSize2OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.Add(jobSize1OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric.Add(urgentMetric2OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric.Add(urgentMetric1OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.Add(dueTime2OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.Add(dueTime1OP9);


                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].JobsIndice.RemoveAt(jobIndex1OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].JobsIndice.RemoveAt(jobIndex2OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].SizeOfJobs.RemoveAt(jobIndex1OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].SizeOfJobs.RemoveAt(jobIndex2OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].UrgentMetric.RemoveAt(jobIndex1OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].UrgentMetric.RemoveAt(jobIndex2OP9);

                                        nonEmptyBatchesAfterOPs[selectedBatchIndex1OP9].DueTime.RemoveAt(jobIndex1OP9);
                                        nonEmptyBatchesAfterOPs[selectedBatchIndex2OP9].DueTime.RemoveAt(jobIndex2OP9);


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
                                newBatchOP10.BatchID = -1;

                                #endregion

                                int selectedBatchIndexOP10 = r.Next(BatchesWithJobsGreaterthan2KminOP10.Length);

                                if (!selectedBatchOP10[selectedBatchIndexOP10])
                                    selectedBatchOP10[selectedBatchIndexOP10] = true;

                                int selectedBatchLengthOP10 =
                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice.Count;

                                int selectedBatchFamilyOP10 =
                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].Family;

                                if (virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count == 0) break;

                                int nOP10 = r.Next(kMin,
                                    (selectedBatchLengthOP10 - kMin) +
                                    Math.Min((int)(kMax / virtualBatches[selectedBatchFamilyOP10].SizeOfJobs.Average()),
                                        virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count));

                                int selectedVirtualBatchLengthOP10 =
                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice.Count;


                                bool[] selectedJobFromSelectedBatchOP10 = new bool[selectedBatchLengthOP10];

                                bool[] selectedJobFromSelectedVirtualBatchOP10 = new bool[selectedVirtualBatchLengthOP10];

                                int jobOP10 = 0, jobIndexOP10 = 0, numberOfUnSelectedVirtualBatchesOP10 = 0;

                                for (int j = 0; j < nOP10;)
                                {
                                    if (numberOfUnSelectedVirtualBatchesOP10 > 0)
                                        j += numberOfUnSelectedVirtualBatchesOP10;
                                    else
                                        j++;

                                    if (selectedVirtualBatchLengthOP10 > 0 && selectedJobFromSelectedVirtualBatchOP10.Any(item => !item))
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
                                            numberOfUnSelectedVirtualBatchesOP10++;
                                        }

                                        if (selectedJobFromSelectedVirtualBatchOP10.Any(item => !item)) continue;


                                    }


                                    else if (selectedBatchLengthOP10 > 0)
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


                                }


                                List<Batch> nonEmptyBatchesAverageDjToWjOP10 = new List<Batch>();

                                foreach (var item in nonEmptyBatchesAfterOPs)
                                {
                                    nonEmptyBatchesAverageDjToWjOP10.Add(item);
                                }

                                if (newBatchOP10.JobsIndice.Count >= kMin)
                                {
                                    foreach (var job in newBatchOP10.JobsIndice)
                                    {
                                        int indexOP10;
                                        double maxP12, maxP22;

                                        #region Remove from VirtualBatches
                                        if (virtualBatches[selectedBatchFamilyOP10].JobsIndice.Any(vj => vj == job))
                                        {
                                            indexOP10 =
                                                Array.IndexOf(
                                                    virtualBatches[selectedBatchFamilyOP10].JobsIndice.ToArray(), job);

                                            virtualBatches[selectedBatchFamilyOP10].JobsIndice.RemoveAt(indexOP10);

                                            virtualBatches[selectedBatchFamilyOP10].SizeOfJobs.RemoveAt(indexOP10);

                                            virtualBatches[selectedBatchFamilyOP10].UrgentMetric.RemoveAt(indexOP10);

                                            virtualBatches[selectedBatchFamilyOP10].DueTime.RemoveAt(indexOP10);

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
                                            indexOP10 =
                                                Array.IndexOf(
                                                    BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10]
                                                        .JobsIndice.ToArray(), job);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].JobsIndice
                                                .RemoveAt(
                                                    indexOP10);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].SizeOfJobs
                                                .RemoveAt(
                                                    indexOP10);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].UrgentMetric
                                                .RemoveAt(
                                                    indexOP10);

                                            BatchesWithJobsGreaterthan2KminOP10[selectedBatchIndexOP10].DueTime.RemoveAt
                                            (
                                                indexOP10);

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

                                    double average = 0;
                                    foreach (var item in newBatchOP10.JobsIndice)
                                    {
                                        average += d[item] / wTj[item];

                                        selectedJobsAfterOPs[item] = true;

                                    }
                                    average = average / (double)newBatchOP10.JobsIndice.Count;

                                    newBatchOP10.AverageDueTimeofJobToDelayImportanceFactor = average;

                                    double averageUpdate = 0;

                                    foreach (var item in nonEmptyBatchesAverageDjToWjOP10)
                                    {
                                        for (int j = 0; j < item.JobsIndice.Count; j++)
                                        {
                                            averageUpdate += d[item.JobsIndice[j]] / wTj[item.JobsIndice[j]];
                                        }
                                        averageUpdate = averageUpdate / (double)item.JobsIndice.Count;

                                        item.AverageDueTimeofJobToDelayImportanceFactor = averageUpdate;
                                    }

                                    bool isSet = false;

                                    int l = 0;

                                    for (int j = 0; j < nonEmptyBatchesAverageDjToWjOP10.Count; j++)
                                    {
                                        newBatchOP10.BatchID = nonEmptyBatchesAverageDjToWjOP10.Count;
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

                        solAfterOPs = Helper.Algorithm1(5, nonEmptyBatchesAfterOPs, t1, t2, Tj, d, t_now);

                        modelAfterOPs.DelayOfJobs = solAfterOPs.Tj;

                        double cost = modelAfterOPs.CostFunction();

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


                    }

                    #endregion


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

                                    mJXFlag[nonemptybatch.JobsIndice[l], nonemptybatch.JobsIndice[j]] = true;
                                }
                            }

                            // tauJB[nonemptybatch.JobsIndice[l]] = sumOfPhiJXs / nonemptybatch.JobsIndice.Count;
                        }
                    }

                    ant[k].mJXFlag = mJXFlag;

                    if (ant[k].Cost < bestAntPerIteration.Cost)
                    {
                        bestAntPerIteration = ant[k];

                        for (int j = 0; j < selectedJobs.Length; j++)
                            bestAntPerIteration.R[j] = !bestAntPerIteration.SelectedJobs[j]
                                ? bestAntPerIteration.R[j] + 1
                                : bestAntPerIteration.R[j];

                        //for (int j = 0; j < tauJ.Length; j++)
                        //    tauJ[j] = (double)1 / (double)(bestAntPerIteration.R[j] + 1);

                    }

                    if (bestAntPerIteration.Cost < bestAnt.Cost)
                    {
                        bestAnt = bestAntPerIteration;
                    }

                    for (int j = 0; j < selectedJobs.Length; j++)
                        R[j] = !selectedJobs[j] ? R[j] + 1 : R[j];

                    for (int j = 0; j < tauJ.Length; j++)
                        tauJ[j] = !selectedJobs[j] ? rhoL * ((double)1 / (double)(R[j] + 1)) : ((double)(1 - rhoL) * tauJ[j]) + (rhoL * tauJ0[j]);

                    for (int i = 0; i < phiJX.GetLength(0); i++)
                    {
                        for (int j = 0; j < phiJX.GetLength(1); j++)
                        {
                            phiJX[i, j] = mJXFlag[i, j]
                                ? (double)(1 - rhoL) * phiJX[i, j] + rhoL * phiJX0[i, j]
                                : (double)(1 - rhoL) * phiJX[i, j];

                        }
                    }


                }

                #region Update Phromone

                for (int i = 0; i < phiJX.GetLength(0); i++)
                {
                    for (int j = 0; j < phiJX.GetLength(1); j++)
                    {
                        phiJX[i, j] = bestAntPerIteration.mJXFlag[i, j]
                            ? ((double)(1 - rho) * phiJX[i, j]) + ((double)mJX[i, j] * phiJX0[i, j])
                            : ((double)(1 - rho) * phiJX[i, j]);
                    }
                }

                //for (int i = 0; i < phiJX.GetLength(0); i++)
                //{
                //    for (int j = 0; j < phiJX.GetLength(1); j++)
                //    {
                //        phiJX[i, j] *= (double)(1 - rho);
                //    }
                //}

                for (int j = 0; j < tauJ.Length; j++)
                {
                    //tauJ[j] = (tauJ[j] * (double)(1 - rho)) + (rho * (double)1 / (double)(R[j] + 1));
                    tauJ[j] = !bestAntPerIteration.SelectedJobs[j]
                        ? (tauJ[j] * (double)(1 - rho))
                        : (tauJ[j] * (double)(1 - rho)) + (tauJ0[j] * (double)1 / (double)(R[j] + 1));
                }
                //for (int j = 0; j < tauJ.Length; j++)
                //    tauJ[j] *= (double)(1 - rho);

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

            string pathToExcelFile = "D:\\504.xls";
            //string pathToExcelFile = Console.ReadLine();


            while ((!File.Exists(pathToExcelFile)))
            {
                Console.WriteLine("File does not exist! ");

                Console.Write("Enter the File Path: ");

                pathToExcelFile = Console.ReadLine();
            }

            Console.WriteLine("Enter the stop time:");

            //long stopElapsedTime = Convert.ToInt64(Console.ReadLine());

            long stopElapsedTime = 50000;
            Run_ACO(pathToExcelFile, stopElapsedTime);


            Console.ReadLine();
        }
    }
}
