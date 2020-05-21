using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskScheduling
{
    public class Helper
    {


        public static int RouletteWheelSelection(double[] pr)
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


        public static bool CheckInSameFamily(int[,] belong2family, int job1, int job2, int family)
        {
            int[,] hjf = belong2family;

            if (hjf[family, job1] == 1 && hjf[family, job2] == 1)
            {
                return true;
            }
            return false;
        }



        public static bool CheckExistInBatch(Batch batch, int jobIndex)
        {
            return batch.JobsIndice.Any(t => t == jobIndex);
        }

        public static Sol Algorithm1(int option, List<Batch> noneEmptybatches, double[] t1, double[] t2, double[] Tj, double[] d,
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

                #region Case 6

                case 6:

                    List<Batch> nonEmptyBatchesSortedByMeanUrgentC6 = noneEmptybatches;

                    int count = nonEmptyBatchesSortedByMeanUrgentC6.Count;

                    for (int bb = 0; bb < count; bb++)
                    {
                        for (int b = 0; b < nonEmptyBatchesSortedByMeanUrgentC6.Count; b++)
                        {
                            for (int j = 0; j < nonEmptyBatchesSortedByMeanUrgentC6[b].UrgentMetric.Count; j++)
                            {
                                nonEmptyBatchesSortedByMeanUrgentC6[b].UrgentMetric[j] =
                                    (double)(d[nonEmptyBatchesSortedByMeanUrgentC6[b].JobsIndice[j]] - t_Now) /
                                    (double)(nonEmptyBatchesSortedByMeanUrgentC6[b].Pbs[0] +
                                              nonEmptyBatchesSortedByMeanUrgentC6[b].Pbs[1] +
                                              nonEmptyBatchesSortedByMeanUrgentC6[b].idleTime);
                            }
                        }

                        nonEmptyBatchesSortedByMeanUrgentC6 =
                            nonEmptyBatchesSortedByMeanUrgentC6.OrderBy(a => a.UrgentMetric.Average()).ToList();


                        double T1 = t1.Min(a => a);

                        int minIndexT1 = Array.IndexOf(t1, T1);

                        double T2 = t2.Min(a => a);

                        int minIndexT2 = Array.IndexOf(t2, T2);

                        double T = T2 - T1;

                        if (nonEmptyBatchesSortedByMeanUrgentC6[0].Pbs[0] - T >= 0)
                        {
                            t1[minIndexT1] += nonEmptyBatchesSortedByMeanUrgentC6[0].Pbs[0];

                            t2[minIndexT2] = t1[minIndexT1] + nonEmptyBatchesSortedByMeanUrgentC6[0].Pbs[1];
                        }
                        else
                        {
                            t1[minIndexT1] = t2[minIndexT2];

                            t2[minIndexT2] += nonEmptyBatchesSortedByMeanUrgentC6[0].Pbs[1];
                        }

                        t_Now = t1[minIndexT1];

                        nonEmptyBatchesSortedByMeanUrgentC6[0].machineNumber[0] = minIndexT1;

                        nonEmptyBatchesSortedByMeanUrgentC6[0].machineNumber[1] = minIndexT2;

                        foreach (int j in nonEmptyBatchesSortedByMeanUrgentC6[0].JobsIndice)
                            Tj[j] = Math.Max((double)t2[minIndexT2] - d[j], 0.0);

                        T1 = t1.Min(a => a);

                        T2 = t2.Min(a => a);

                        T = T2 - T1;

                        for (int i = 1; i < nonEmptyBatchesSortedByMeanUrgentC6.Count; i++)
                        {

                            if (nonEmptyBatchesSortedByMeanUrgentC6[i].Pbs[0] - T >= 0)
                            {
                                nonEmptyBatchesSortedByMeanUrgentC6[i].idleTime = 0;
                            }
                            else
                            {
                                nonEmptyBatchesSortedByMeanUrgentC6[i].idleTime =
                                    T - nonEmptyBatchesSortedByMeanUrgentC6[i].Pbs[0];
                            }

                        }


                        sol.BatchesAllocatedToMachines.Add(nonEmptyBatchesSortedByMeanUrgentC6[0]);

                        nonEmptyBatchesSortedByMeanUrgentC6.RemoveAt(0);
                    }

                    sol.TimeofMachinesStep1 = t1;

                    sol.TimeofMachinesStep2 = t2;

                    break;

                    #endregion
            }

            sol.Tj = Tj;

            return sol;
        }
    }
}
