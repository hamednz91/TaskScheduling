using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskScheduling
{

    public class Batch
    {
        public List<int> BatchCandidateList { get; set; }

        public List<double> Eta1jCandidates { get; set; }

        public List<double> Eta2jCandidates { get; set; }

        public List<double> Eta3jCandidates { get; set; }

        public List<double> Eta4jCandidates { get; set; }

        public List<double> TauJBCandidates { get; set; }

        public List<double> TauJCandidates { get; set; }

        public List<double> JobCandidateSelectionProbability { get; set; }

        public List<int> JobsIndice { get; set; }

        public List<int> SizeOfJobs { get; set; }

        public List<double> UrgentMetric { get; set; }

        public List<double> DueTime { get; set; } //d[j] Due time of 

        public int[] machineNumber { get; set; } //machineNumber[in_step1,in_step2

        public double[] Pbs { get; set; } //batch proccessing time in step 1 and 2 [step1,step2]

        public int Family { get; set; }

        public int batchIndex { get; set; }

        public double AverageDueTimeofJobToDelayImportanceFactor { get; set; } // dj/wj

        public double idleTime { get; set; }

        public int BatchID { get; set; }

    }

}
