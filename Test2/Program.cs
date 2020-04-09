using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Test2
{
    class Program
    {
  
        static void Main(string[] args)
        {
            int n, a = 0, b = 1, c, i;
            Console.Write("enter numer of series");
            n = int.Parse(Console.ReadLine());
            Console.Write("{0,5}{1,5} ", a, b);
            for (i = 3; i <= n; ++i)
            {
                c = a + b;
                Console.Write("{0,5}", c);
                a = b;
                b = c;

            }
            Console.ReadKey();

        }
    }
}
