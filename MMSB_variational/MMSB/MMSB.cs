using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Models;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Factors;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Models
{
    class MMSB
    {
        static void Main(string[] args)
        {
            // Observed interaction matrix
            var YObs = new bool[20][];
            YObs[0] = new bool[] { true, false, true, false, true, false, false, true, false, false, false, false, true, false, false, false, false, false, false, true };
            YObs[1] = new bool[] { true, true, false, true, false, true, false, false, true, true, true, true, false, true, false, true, true, true, true, true };
            YObs[2] = new bool[] { true, false, true, false, true, false, false, true, false, false, true, false, true, false, false, true, false, false, false, false };
            YObs[3] = new bool[] { true, true, false, true, false, true, true, true, false, false, true, true, true, true, true, true, true, true, false, true };
            YObs[4] = new bool[] { false, false, true, false, true, false, true, false, true, false, true, true, true, true, true, false, false, true, false, false };
            YObs[5] = new bool[] { false, true, false, false, false, true, true, false, false, true, false, true, true, true, false, true, true, false, true, true };
            YObs[6] = new bool[] { false, false, false, true, false, true, true, true, true, true, true, true, false, true, false, true, false, true, true, true };
            YObs[7] = new bool[] { false, true, false, false, false, true, true, true, false, true, true, true, false, true, false, true, true, true, true, false };
            YObs[8] = new bool[] { false, true, false, true, false, true, true, true, true, false, true, true, true, true, false, true, true, true, true, true };
            YObs[9] = new bool[] { false, false, true, true, true, false, false, true, true, true, true, true, false, true, false, true, false, false, true, true };
            YObs[10] = new bool[] { false, true, false, true, false, true, true, true, true, true, true, false, false, true, false, false, false, false, true, true };
            YObs[11] = new bool[] { false, true, false, true, false, true, true, true, true, true, true, true, false, false, true, true, false, false, true, true };
            YObs[12] = new bool[] { false, false, true, false, true, true, false, true, false, false, false, false, true, false, true, false, false, true, false, false };
            YObs[13] = new bool[] { false, true, false, true, false, true, false, false, false, true, true, true, false, true, true, true, true, true, true, true };
            YObs[14] = new bool[] { true, false, false, false, true, false, false, false, false, false, false, true, true, false, false, false, true, false, false, true };
            YObs[15] = new bool[] { true, true, false, false, false, false, false, true, false, true, true, true, false, false, true, true, true, false, true, true };
            YObs[16] = new bool[] { true, false, false, true, false, false, true, false, false, true, true, true, false, true, true, true, true, true, false, true };
            YObs[17] = new bool[] { false, true, true, true, true, true, true, true, true, true, false, false, false, true, false, true, true, false, false, true };
            YObs[18] = new bool[] { false, true, false, true, false, true, false, true, true, false, true, true, false, true, false, false, true, true, false, true };
            YObs[19] = new bool[] { false, true, false, false, false, false, true, false, true, true, true, true, false, false, false, true, false, true, true, true };
            int K = 2;           // Number of blocks
            int N = YObs.Length; // Number of nodes

            // Ranges
            Range p = new Range(N).Named("p");   // Range for initiator
            Range q = p.Clone().Named("q");      // Range for receiver
            Range kp = new Range(K).Named("kp"); // Range for initiator block membership
            Range kq = kp.Clone().Named("kq");   // Range for receiver block membership

            // The model
            var Y = Variable.Array(Variable.Array<bool>(q), p); // Interaction matrix
            var pi = Variable.Array<Vector>(p).Named("pi");     // Block-membership probability vector
            pi[p] = Variable.DirichletUniform(kp).ForEach(p);
            var B = Variable.Array<double>(kp, kq).Named("B");  // Link probability matrix
            B[kp, kq] = Variable.Beta(1, 1).ForEach(kp, kq);

            using (Variable.ForEach(p))
            {
                using (Variable.ForEach(q))
                {
                    var z1 = Variable.Discrete(pi[p]).Named("z1"); // Draw initiator membership indicator
                    var z2 = Variable.Discrete(pi[q]).Named("z2"); // Draw receiver membership indicator
                    z2.SetValueRange(kq);
                    using (Variable.Switch(z1))
                    using (Variable.Switch(z2))
                        Y[p][q] = Variable.Bernoulli(B[z1, z2]);   // Sample interaction value
                }
            }

            // Initialise to break symmetry
            var piInit = new Dirichlet[N];
            for (int i = 0; i < N; i++)
            {
                Vector v = Vector.Zero(K);
                for (int j = 0; j < K; j++) v[j] = Rand.Double(); //10 * 
                piInit[i] = new Dirichlet(v);
            }

            // Hook up the data
            Y.ObservedValue = YObs;

            // Infer
            var engine = new InferenceEngine(new VariationalMessagePassing());
            pi.InitialiseTo(Distribution<Vector>.Array(piInit));
            var posteriorPi = engine.Infer<Dirichlet[]>(pi);
            var posteriorB = engine.Infer<Beta[,]>(B);

            Console.WriteLine("\n\n************************************");
            Console.WriteLine("The Value of the Pi's is This:");

            for (int i = 0; i < posteriorPi.Length; i++)
            {
                Console.WriteLine(posteriorPi[i].GetMean());
            }


            Console.WriteLine("\n\n************************************");
            Console.WriteLine("\nThe value of Beta is \n");
 
            for (int i = 0; i < K; i++)
			{
				for (int j = 0; j < K; j++)
				{
					Console.Write(posteriorB[i,j] + "  ");
				}
 
				Console.WriteLine("\n");
			}

            //Console.WriteLine(
            //    String.Format("\nA random double is here: {0}\n",
            //    Rand.Double()));


            Console.ReadLine();
        }
    }
}
