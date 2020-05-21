using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TaskScheduling
{
    public class Sol
    {

        public List<Batch> BatchesAllocatedToMachines { get; set; }

        public double[] TimeofMachinesStep1 { get; set; }

        public double[] TimeofMachinesStep2 { get; set; }

        public double[] Tj { get; set; }

        public Sol()
        {
            BatchesAllocatedToMachines = new List<Batch>();
        }

        public Sol(Sol sol)
        {
            this.BatchesAllocatedToMachines = sol.BatchesAllocatedToMachines;
            this.TimeofMachinesStep1 = sol.TimeofMachinesStep1;
            this.TimeofMachinesStep2 = sol.TimeofMachinesStep2;
            this.Tj = sol.Tj;
        }

    }
}
