// (C) Copyright 2009 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Collections;

using GaussianArray = MicrosoftResearch.Infer.Distributions.DistributionStructArray<MicrosoftResearch.Infer.Distributions.Gaussian, double>;



namespace MicrosoftResearch.Infer.Factors
{
    public class IndexOfMaximumBuffer
    {
        public IList<Gaussian> MessagesToMax;
        public IList<Gaussian> to_list;
    }


    /// <summary>
    /// TODO: comments
    /// </summary>
    [FactorMethod(typeof(MMath), "IndexOfMaximumDouble")]
    [Quality(QualityBand.Experimental)]
    [Buffers("Buffer")]
    public static class IndexOfMaximumOp
    {

        public static IndexOfMaximumBuffer BufferInit<GaussianList>([IgnoreDependency] GaussianList list)
             where GaussianList : IList<Gaussian>
        {
            return new IndexOfMaximumBuffer
            {
                MessagesToMax = new DistributionStructArray<Gaussian, double>(list.Select(o => Gaussian.Uniform()).ToArray()),
                to_list = new DistributionStructArray<Gaussian, double>(list.Select(o => Gaussian.Uniform()).ToArray())
            }; 
        }

        // redundant parameters required for correct dependency graph
        public static IndexOfMaximumBuffer Buffer<GaussianList>(IndexOfMaximumBuffer Buffer, GaussianList list, int IndexOfMaximumDouble)
            where GaussianList : IList<Gaussian>
        {
            var max_marginal = Buffer.to_list[IndexOfMaximumDouble] * list[IndexOfMaximumDouble];
            //var order = Rand.Perm(list.Count); 
            for (int i = 0; i < list.Count; i++)
            {
                //int c = order[i]; 
                int c = i;
                if (c != IndexOfMaximumDouble)
                {
                    var msg_to_sum = max_marginal / Buffer.MessagesToMax[c];

                    var msg_to_positiveop = DoublePlusOp.AAverageConditional(Sum: msg_to_sum, b: list[c]);
                    var msgFromPositiveOp = IsPositiveOp.XAverageConditional(true, msg_to_positiveop);
                    Buffer.MessagesToMax[c] = DoublePlusOp.SumAverageConditional(list[c], msgFromPositiveOp);
                    Buffer.to_list[c] = DoublePlusOp.AAverageConditional(Sum: msg_to_sum, b: msgFromPositiveOp);
                    max_marginal = msg_to_sum * Buffer.MessagesToMax[c];
                }
            }
            Buffer.to_list[IndexOfMaximumDouble] = max_marginal / list[IndexOfMaximumDouble];
            return Buffer;
        }

        public static GaussianList listAverageConditional<GaussianList>(IndexOfMaximumBuffer Buffer, GaussianList to_list, int IndexOfMaximumDouble)
            where GaussianList : IList<Gaussian>
        {
            to_list.SetTo(Buffer.to_list);
            return to_list;
        }

        public static double LogAverageFactor<GaussianList>(IndexOfMaximumBuffer Buffer, GaussianList list, int IndexOfMaximumDouble)
            where GaussianList : IList<Gaussian>
        {
            double evidence = 0;
            var max_marginal = list[IndexOfMaximumDouble] * Buffer.to_list[IndexOfMaximumDouble];
            for (int c = 0; c < list.Count; c++)
            {
                if (c != IndexOfMaximumDouble)
                {
                    var msg_to_sum = max_marginal / Buffer.MessagesToMax[c];
                    var msg_to_positiveop = DoublePlusOp.AAverageConditional(Sum: msg_to_sum, b: list[c]);
                    evidence += IsPositiveOp.LogEvidenceRatio(true, msg_to_positiveop);
                    // sum operator does not contribute because no projection is involved
                    // the x[index]-x[c] variable does not contribute because it only connects to two factors
                }
            }
            //evidence += ReplicateOp.LogEvidenceRatio<Gaussian>(MessagesToMax, list[IndexOfMaximumDouble], MessagesToMax.Select(o => max_marginal / o).ToArray());
            evidence += max_marginal.GetLogNormalizer() - Buffer.MessagesToMax.Sum(o => o.GetLogNormalizer()) - list[IndexOfMaximumDouble].GetLogNormalizer();
            evidence -= Buffer.MessagesToMax.Sum(o => (max_marginal / o).GetLogAverageOf(o));
            return evidence;
        }

        public static double LogEvidenceRatio<GaussianList>(IndexOfMaximumBuffer Buffer, GaussianList list, int IndexOfMaximumDouble)
    where GaussianList : IList<Gaussian>
        {
            return LogAverageFactor(Buffer, list, IndexOfMaximumDouble); 
        }

    }


    /// <summary>
    /// TODO: comments
    /// </summary>
    [FactorMethod(typeof(MMath), "IndexOfMaximumDouble")]
    [Quality(QualityBand.Experimental)]
    [Buffers("Buffers")]
    public static class IndexOfMaximumStochasticOp
    {
        public static IndexOfMaximumBuffer[] BuffersInit<GaussianList>([IgnoreDependency] GaussianList list, [IgnoreDependency] Discrete IndexOfMaximumDouble)
             where GaussianList : IList<Gaussian>
        {
            return list.Select(o => IndexOfMaximumOp.BufferInit(list)).ToArray(); 
        }

        public static IndexOfMaximumBuffer[] Buffers<GaussianList>(IndexOfMaximumBuffer[] Buffers, [SkipIfUniform] GaussianList list, Discrete IndexOfMaximumDouble)
            where GaussianList : IList<Gaussian>
        {
            for (int i = 0; i < list.Count; i++)
            {
                Buffers[i].to_list = new DistributionStructArray<Gaussian, double>(list.Select(o => Gaussian.Uniform()).ToArray());
                Buffers[i].to_list[i] = Buffers[i].MessagesToMax.Aggregate((p, q) => p * q);
                Buffers[i]=IndexOfMaximumOp.Buffer(Buffers[i], list, i);
            }
            return Buffers;
        }

        public static GaussianList listAverageConditional<GaussianList>([Fresh, SkipIfUniform] IndexOfMaximumBuffer[] Buffers, GaussianList list, [SkipIfUniform] Discrete IndexOfMaximumDouble, GaussianList result)
            where GaussianList : DistributionStructArray<Gaussian, double> // IList<Gaussian>
        {
            // TODO: check if Index is a point mass
            var results = list.Select(o => list.Select(p => Gaussian.Uniform()).ToList()).ToArray();
            //var evidences = new Bernoulli[list.Count];
            for (int i = 0; i < list.Count; i++)
            {
                for (int j = 0; j < list.Count; j++)
                {
                    results[j][i] = Buffers[i].to_list[j];
                }
            }
            for (int i = 0; i < list.Count; i++)
            {
                result[i] = GateEnterPartialOp<double>.ValueAverageConditional<Gaussian>(results[i], IndexOfMaximumDouble, list[i], System.Linq.Enumerable.Range(0, list.Count).ToArray(), result[i]);
            }
            return result;
        }

        public static Discrete IndexOfMaximumDoubleAverageConditional<GaussianList>(GaussianList list, IndexOfMaximumBuffer[] Buffers)
            where GaussianList : DistributionStructArray<Gaussian, double>
        {
            // var results = list.Select(o => list.Select(p => Gaussian.Uniform()).ToList()).ToArray();
            // TODO: if IndexOfMaximumDouble is uniform we will never call this routine so buffers will not get set, so messages to IndexOfMaximumDouble will be incorrect
            var evidences = new double[list.Count];
            for (int i = 0; i < list.Count; i++)
            {
                //var res = new DistributionStructArray<Gaussian, double>(list.Select(o => Gaussian.Uniform()).ToArray());
                //res[i] = Buffer[i].MessagesToMax.Aggregate((p, q) => p * q);
                evidences[i] = IndexOfMaximumOp.LogAverageFactor(Buffers[i], list, i); 
            }
            return new Discrete(MMath.Softmax(evidences));
        }

        public static double LogAverageFactor<GaussianList>([SkipIfUniform] GaussianList list, GaussianList to_list, IndexOfMaximumBuffer[] Buffers, Discrete IndexOfMaximumDouble)
            where GaussianList : DistributionStructArray<Gaussian, double>
        {
            var evidences = new double[list.Count];
            var tempBuffer = new IndexOfMaximumBuffer(); 
            for (int i = 0; i < list.Count; i++)
            {
                tempBuffer.to_list = new DistributionStructArray<Gaussian, double>(list.Select(o => Gaussian.Uniform()).ToArray());
                tempBuffer.to_list[i] = Buffers[i].MessagesToMax.Aggregate((p, q) => p * q);
                tempBuffer.MessagesToMax = new DistributionStructArray<Gaussian, double>(Buffers[i].MessagesToMax.Select(o => (Gaussian)o.Clone()).ToArray());
                evidences[i] = IndexOfMaximumOp.LogAverageFactor(tempBuffer, list, i) + IndexOfMaximumDouble.GetLogProb(i); ; 
            }
            return MMath.LogSumExp(evidences);
        }

        public static double LogEvidenceRatio<GaussianList>(GaussianList list, GaussianList to_list, Discrete IndexOfMaximumDouble, Discrete to_IndexOfMaximumDouble, IndexOfMaximumBuffer[] Buffers)
            where GaussianList : DistributionStructArray<Gaussian, double>
        {
            return LogAverageFactor(list, to_list, Buffers, IndexOfMaximumDouble) - IndexOfMaximumDouble.GetLogAverageOf(to_IndexOfMaximumDouble);
        }
    }

}
