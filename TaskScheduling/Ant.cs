using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskScheduling
{
    public class Ant
    {
        public List<Batch> Tour { get; set; }

        public double Cost { get; set; }

        public Sol sol { get; set; }

        public bool[] SelectedJobs { get; set; }

        public bool[] SelectedJobsForEmptyBatches { get; set; }

        public double[] R { get; set; }

        public bool[,] mJXFlag { get; set; }

    }

}
