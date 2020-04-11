using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TestTaskScheduling
{
    class Program
    {
        static void Main(string[] args)
        {
            List<double> a=new List<double>();
            a.Add(14.2);
            a.Add(112.2);
            a.Add(111.5);
            a.Add(1.223);
            a.Add(4.2);


            double x = 13.2;
            bool isSet = false;
            for (int i = 0; i < a.Count; i++)
            {
                if (x>a[i])
                {
                    a.Insert(i+1,x);
                    isSet = true;
                    for (int j = i+2; j < a.Count; j++)
                    {
                        a[j]++;
                    }
                }
                
            }
            if (!isSet)
            {
                a.Insert(0,x);
            }
            Console.ReadKey();

        }
    }
}
