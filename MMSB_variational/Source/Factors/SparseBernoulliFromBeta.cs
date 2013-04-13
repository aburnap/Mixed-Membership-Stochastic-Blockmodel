// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="SparseBernoulliList.Sample(SparseVector)"/>, given random arguments to the function.
	/// </summary>
    [FactorMethod(typeof(SparseBernoulliList), "Sample", typeof(SparseVector))]
	[Buffers("MeanLog","MeanLogOneMinus")]
	[Quality(QualityBand.Preview)]
	public class SparseBernoulliFromBetaOp
    {
        //-- VMP -------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanLog">Buffer 'MeanLog'.</param>
		/// <param name="MeanLogOneMinus">Buffer 'MeanLogOneMinus'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,probsTrue) p(sample,probsTrue) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseBernoulliListBase sample,
			[Proper] SparseBetaList probsTrue, SparseVector MeanLog, SparseVector MeanLogOneMinus)
        {
			//var MeanLogOneMinus = probsTrue.GetMeanLogOneMinus();
            var p = sample.GetProbTrueVector();
			var res = p * MeanLog + (1 - p) * MeanLogOneMinus;
            return res.Sum();
        }

		/// <summary>
		/// Update the buffer 'MeanLogOneMinus'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'MeanLogOneMinus'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static SparseVector MeanLogOneMinus([Proper] SparseBetaList probsTrue)
		{
			return probsTrue.GetMeanLogOneMinus();
		}

		/// <summary>
		/// Update the buffer 'MeanLog'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'MeanLog'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static SparseVector MeanLog([Proper] SparseBetaList probsTrue)
		{
			return probsTrue.GetMeanLog();
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="probsTrue">Constant value for 'probsTrue'.</param>
		/// <param name="MeanLog">Buffer 'MeanLog'.</param>
		/// <param name="MeanLogOneMinus">Buffer 'MeanLogOneMinus'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseBernoulliListBase sample, SparseVector probsTrue, SparseVector MeanLog, SparseVector MeanLogOneMinus)
        {
            return AverageLogFactor(sample, SparseBetaList.PointMass(probsTrue),MeanLog, MeanLogOneMinus);
        }

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanLogOneMinus">Buffer 'MeanLogOneMinus'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor(IList<bool> sample, [Proper] SparseBetaList probsTrue, SparseVector MeanLogOneMinus)
        {
			//var MeanLogOneMinus = probsTrue.GetMeanLogOneMinus();
			double ev = MeanLogOneMinus.Sum();
			int count = sample.Count;
			// assume 'sample' is very sparse
			bool[] sampleAsArray = (bool[])sample;
            for (int i = 0; i < count; i++)
            {
                if (!sampleAsArray[i]) continue;
				ev -= MeanLogOneMinus[i];
                ev += SparseBetaList.ComputeMeanLog(probsTrue.TrueCounts[i],probsTrue.FalseCounts[i]);
            }
            return ev;
        }

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probsTrue">Constant value for 'probsTrue'.</param>
		/// <param name="MeanLogOneMinus">Buffer 'MeanLogOneMinus'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(IList<bool> sample, SparseVector probsTrue, SparseVector MeanLogOneMinus)
        {
			SparseVector tv = SparseVector.Zero(probsTrue.Count);
			tv.SetToFunction(probsTrue, x => Math.Log(1-x));
            double ev = tv.Sum();
			int count = sample.Count;
			bool[] sampleAsArray = (bool[])sample;
			for (int i = 0; i < sampleAsArray.Length; i++)
            {
				if (!sampleAsArray[i]) continue;
                double pTrue = probsTrue[i];
                ev -= Math.Log(1 - pTrue);
                ev += Math.Log(pTrue);
            }
            return ev;
        }

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probsTrue">Constant value for 'probsTrue'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static SparseBernoulliList SampleAverageLogarithm(SparseVector probsTrue)
        {
            var logOdds = probsTrue.Clone();
			logOdds.SetToFunction(probsTrue, MMath.Logit);
            return new SparseBernoulliList(logOdds);
        }

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static SparseBernoulliList SampleAverageLogarithm([SkipIfUniform] SparseBetaList probsTrue)
        {
			var logOdds = SparseVector.Zero(probsTrue.Count);
			logOdds.SetToFunction(probsTrue.TrueCounts, probsTrue.FalseCounts, ComputeLogOdds);
            return new SparseBernoulliList(logOdds);
        }

        /// <summary>
        /// Used to compute log odds in the above operator
        /// </summary>
        /// <param name="trueCount"></param>
        /// <param name="falseCount"></param>
        /// <returns></returns>
        internal static double ComputeLogOdds(double trueCount, double falseCount)
        {
            if (falseCount == Double.PositiveInfinity)
            {
                // compute log odds from prob true
                return MMath.Logit(trueCount);
            }
            else if ((trueCount == 0) || (falseCount == 0))
            {
                throw new ImproperMessageException(new Beta(trueCount, falseCount));
            }
            return MMath.Digamma(trueCount) - MMath.Digamma(falseCount);
        }

		/// <summary>
		/// VMP message to 'probsTrue'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing VMP message to the 'probsTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'probsTrue' conditioned on the given values.
		/// </para></remarks>
		public static SparseBetaList ProbsTrueAverageLogarithm(IList<bool> sample)
        {
			int count = sample.Count;
			var trueCounts = SparseVector.Constant(count, 1.0);
			var falseCounts = SparseVector.Constant(count, 2.0);
			bool[] sampleAsArray = (bool[])sample;
            for (int i = 0; i < sampleAsArray.Length; i++)
            {
				if (!sampleAsArray[i]) continue;
                trueCounts[i] = 2;
                falseCounts[i] = 1;
            }            
            return new SparseBetaList(trueCounts,falseCounts);
        }

		/// <summary>
		/// VMP message to 'probsTrue'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'probsTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'probsTrue'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,probsTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static SparseBetaList ProbsTrueAverageLogarithm([Proper] SparseBernoulliListBase sample)
        {
            // E[x*log(p) + (1-x)*log(1-p)] = E[x]*log(p) + (1-E[x])*log(1-p)
            SparseVector ex = sample.GetProbTrueVector();
			return new SparseBetaList((SparseVector)(ex + 1), (SparseVector)(2 - ex));
        }
    }

	/// <summary>
	/// Provides outgoing messages for <see cref="BernoulliIntegerSubset.Sample(SparseVector)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(BernoulliIntegerSubset), "Sample", typeof(SparseVector))]
	[Buffers("SumMeanLogOneMinus", "MeanLogOneMinus", "MeanLogMinusMeanLogOneMinus")]
	[Quality(QualityBand.Experimental)]
	public class BernoulliIntegerSubsetFromBeta
	{
		//-- VMP -------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanLog"></param>
		/// <param name="MeanLogOneMinus"></param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] BernoulliIntegerSubset sample, [Proper] SparseBetaList probsTrue,SparseVector MeanLog, SparseVector MeanLogOneMinus)
		{
			return SparseBernoulliFromBetaOp.AverageLogFactor(sample, probsTrue,MeanLog,MeanLogOneMinus);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="probsTrue">Constant value 'probsTrue'.</param>
		/// <param name="MeanLog"></param>
		/// <param name="MeanLogOneMinus"></param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] BernoulliIntegerSubset sample, SparseVector probsTrue, SparseVector MeanLog, SparseVector MeanLogOneMinus)
		{
			return SparseBernoulliFromBetaOp.AverageLogFactor(sample, probsTrue,MeanLog, MeanLogOneMinus);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SumMeanLogOneMinus"></param>
		/// <param name="MeanLogMinusMeanLogOneMinus">Buffer 'MeanLogMinusMeanLogOneMinus'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double AverageLogFactor(IList<int> sample, [Proper] SparseBetaList probsTrue, double SumMeanLogOneMinus, SparseVector MeanLogMinusMeanLogOneMinus)
		{
			//var MeanLogOneMinus = probsTrue.GetMeanLogOneMinus();
			double ev = SumMeanLogOneMinus;
			// assume 'sample' is very sparse
			foreach(int i in sample)
			{
				ev += MeanLogMinusMeanLogOneMinus[i];
			}
			return ev;
		}

		/// <summary>
		/// Update the buffer 'MeanLogOneMinus'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'MeanLogOneMinus'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static SparseVector MeanLogOneMinus([Proper] SparseBetaList probsTrue)
		{
			return probsTrue.GetMeanLogOneMinus();
		}

		/// <summary>
		/// Update the buffer 'SumMeanLogOneMinus'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'SumMeanLogOneMinus'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static double SumMeanLogOneMinus([Proper] SparseBetaList probsTrue)
		{
			return probsTrue.GetMeanLogOneMinus().Sum();
		}

		/// <summary>
		/// Update the buffer 'MeanLogMinusMeanLogOneMinus'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'MeanLogMinusMeanLogOneMinus'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static SparseVector MeanLogMinusMeanLogOneMinus([Proper] SparseBetaList probsTrue)
		{
			return (SparseVector)(probsTrue.GetMeanLog() - probsTrue.GetMeanLogOneMinus());
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="probsTrue">Constant value for 'probsTrue'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,probsTrue))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(IList<int> sample, SparseVector probsTrue)
		{
			SparseVector tv = SparseVector.Zero(probsTrue.Count);
			tv.SetToFunction(probsTrue, x => Math.Log(1-x));
			double ev = tv.Sum();
			int count = sample.Count;
			foreach(int i in sample)
			{
				double pTrue = probsTrue[i];
				ev -= Math.Log(1 - pTrue);
				ev += Math.Log(pTrue);
			}
			return ev;
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probsTrue">Constant value for 'probsTrue'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static BernoulliIntegerSubset SampleAverageLogarithm(SparseVector probsTrue)
		{
			var logOdds = probsTrue.Clone();
			logOdds.SetToFunction(probsTrue, MMath.Logit);
			return new BernoulliIntegerSubset(logOdds);
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="probsTrue">Incoming message from 'probsTrue'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(probsTrue) p(probsTrue) log(factor(sample,probsTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="probsTrue"/> is not a proper distribution</exception>
		public static BernoulliIntegerSubset SampleAverageLogarithm([SkipIfUniform] SparseBetaList probsTrue)
		{
			var logOdds = SparseVector.Zero(probsTrue.Count);
			logOdds.SetToFunction(probsTrue.TrueCounts, probsTrue.FalseCounts, SparseBernoulliFromBetaOp.ComputeLogOdds);
			return new BernoulliIntegerSubset(logOdds);
		}

		/// <summary>
		/// VMP message to 'probsTrue'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'probsTrue' conditioned on the given values.
		/// </para></remarks>
		public static SparseBetaList ProbsTrueAverageLogarithm(IList<int> sample, SparseBetaList result)
		{
			var trueCounts =result.TrueCounts;
			var falseCounts =result.FalseCounts;
			trueCounts.SetAllElementsTo(1.0);
			falseCounts.SetAllElementsTo(2.0);
			foreach(int i in sample)
			{
				trueCounts[i] = 2;
				falseCounts[i] = 1;
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'probsTrue'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'probsTrue' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'probsTrue'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,probsTrue)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static SparseBetaList ProbsTrueAverageLogarithm([Proper] BernoulliIntegerSubset sample)
		{
			return SparseBernoulliFromBetaOp.ProbsTrueAverageLogarithm(sample);
		}
	}
}
