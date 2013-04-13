// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.DiscreteFromLogProbs"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "DiscreteFromLogProbs")]
	[Quality(QualityBand.Experimental)]
	[Obsolete("Use SoftmaxOp", true)]
	public static class DiscreteFromLogProbsOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logProbs">Incoming message from 'logProbs'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(logProbs) p(logProbs) factor(sample,logProbs))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, IList<double> logProbs)
		{
			if(logProbs.Count <= 1) return 0;
			return logProbs[sample] - MMath.LogSumExp(logProbs);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logProbs">Incoming message from 'logProbs'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(logProbs) p(logProbs) log(factor(sample,logProbs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(int sample, IList<double> logProbs)
		{
			return LogAverageFactor(sample, logProbs);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logProbs">Incoming message from 'logProbs'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(logProbs) p(logProbs) log(factor(sample,logProbs))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="logProbs"/> is not a proper distribution</exception>
		public static double AverageLogFactor(int sample, [Proper] IList<Gaussian> logProbs)
		{
			double result = 0;
			double ms,vs;
			logProbs[sample].GetMeanAndVariance(out ms, out vs);
			for (int k = 0; k < logProbs.Count; k++)
			{
				if(k == sample) continue;
				double m,v;
				logProbs[k].GetMeanAndVariance(out m, out v);
				Gaussian logProb = new Gaussian(ms-m,vs+v);
				result += BernoulliFromLogOddsOp.AverageLogFactor(true, logProb);
			}
			return result;
		}

		// the factor is exp(logProbs[sample])/(sum_k exp(logProbs[k]))
		// which is bounded by prod_{k!=sample} sigma(logProbs[sample]-logProbs[k])
		/// <summary>
		/// VMP message to 'logProbs'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logProbs">Incoming message from 'logProbs'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'logProbs' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="logProbs"/> is not a proper distribution</exception>
		public static GaussianList LogProbsAverageLogarithm<GaussianList>(int sample, [Proper] IList<Gaussian> logProbs, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			double sumPrecision = 0;
			double sumMeanPrecision = 0;
			double ms,vs;
			logProbs[sample].GetMeanAndVariance(out ms, out vs);
			for (int k = 0; k < logProbs.Count; k++)
			{
				if(k == sample) continue;
				double m,v;
				logProbs[k].GetMeanAndVariance(out m, out v);
				Gaussian logProb = new Gaussian(ms-m,v+vs);
				Gaussian toLogProb = BernoulliFromLogOddsOp_JJ96.LogOddsAverageLogarithm(false, logProb);
				toLogProb.MeanTimesPrecision += ms*toLogProb.Precision;
				result[k] = toLogProb;
				sumPrecision += toLogProb.Precision;
				sumMeanPrecision += m*toLogProb.Precision;
			}
			result[sample] = Gaussian.FromNatural(0.5 + sumMeanPrecision, sumPrecision);
			return result;
		}
	}
}
