// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Linq;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Multinomial"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>The factor is f(sample,p,n) = n!/prod_k sample[k]!  prod_k p[k]^sample[k]</remarks>
	[FactorMethod(typeof(Factor), "Multinomial", typeof(int), typeof(Vector))]
	[FactorMethod(typeof(Factor), "MultinomialList", typeof(int), typeof(Vector))]
	[Buffers("MeanLog")]
	[Quality(QualityBand.Preview)]
	public static class MultinomialOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,trialCount,p) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(IList<Discrete> sample) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,p) p(sample,p) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(IList<int> sample, int trialCount, IList<double> p)
		{
			double result = MMath.GammaLn(trialCount+1);
			for (int i = 0; i < sample.Count; i++) {
				result += sample[i]*Math.Log(p[i]) - MMath.GammaLn(sample[i]+1);
			}
			return result;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,p) p(sample,p) factor(sample,trialCount,p))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(IList<int> sample, int trialCount, Dirichlet p)
		{
			double result = MMath.GammaLn(trialCount+1);
			for (int i = 0; i < sample.Count; i++) {
				result += MMath.GammaLn(sample[i]+p.PseudoCount[i]) + MMath.GammaLn(p.PseudoCount[i])
					-MMath.GammaLn(sample[i]+1);
			}
			result += MMath.GammaLn(p.TotalCount) - MMath.GammaLn(p.TotalCount + trialCount);
			return result;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,p) p(sample,p) factor(sample,trialCount,p) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(IList<int> sample, int trialCount, Dirichlet p)
        {
            return LogAverageFactor(sample, trialCount, p);
        }
		public static double LogEvidenceRatio(IList<int> sample, int trialCount, IList<double> p)
		{
			return LogAverageFactor(sample, trialCount, p);
		}

		/// <summary>
		/// EP message to 'p'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>The outgoing EP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'p' as the random arguments are varied.
		/// The formula is <c>proj[p(p) sum_(sample) p(sample) factor(sample,trialCount,p)]/p(p)</c>.
		/// </para></remarks>
		public static Dirichlet PAverageConditional(IList<int> sample, int trialCount)
        {
            return PAverageLogarithm(sample, trialCount); 
        }

		//-- VMP ------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,p) p(sample,p) log(factor(sample,trialCount,p))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(IList<int> sample, int trialCount, IList<double> p)
		{
			return LogAverageFactor(sample, trialCount, p);
		}

		/// <summary>
		/// Update the buffer 'MeanLog'
		/// </summary>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <returns>New value of buffer 'MeanLog'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static Vector MeanLog(Dirichlet p)
		{
			return p.GetMeanLog();
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <param name="p">Incoming message from 'p'.</param>
		/// <param name="MeanLog">Buffer 'MeanLog'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,p) p(sample,p) log(factor(sample,trialCount,p))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(IList<int> sample, int trialCount, Dirichlet p, [Fresh] Vector MeanLog)
		{
			double result = MMath.GammaLn(trialCount+1);
			if (true) {
				result += sample.Inner(MeanLog);
				result -= sample.ListSum(x => MMath.GammaLn(x+1));
			} else {
				for (int i = 0; i < sample.Count; i++) {
					result += sample[i]*MeanLog[i] - MMath.GammaLn(sample[i]+1);
				}
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'p'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trialCount">Constant value for 'trialCount'.</param>
		/// <returns>The outgoing VMP message to the 'p' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'p'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,trialCount,p)))</c>.
		/// </para></remarks>
		public static Dirichlet PAverageLogarithm(IList<int> sample, int trialCount)
		{
			// This method demonstrates how to write an operator method using extension methods,
			// so that it is efficient for sparse and dense lists.
			//
			// The vector returned here will be sparse if the list is sparse.
			var counts = sample.ListSelect(x => x+1.0).ToVector();
			return new Dirichlet(counts);
		}
	}
}
