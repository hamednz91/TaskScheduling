using System;

namespace Test
{
    class Program
    {
        static void Main(string[] args)
        {
            int n = 1;
            int sum = 0;
            n = int.Parse(Console.ReadLine());

            for (int i = 1; i <= n; i = i + sum)
            {
                Console.Write("," + i);
                sum = sum + i;
                Console.Write("," + sum);

            }

            Console.ReadKey();



        }
    }
}
