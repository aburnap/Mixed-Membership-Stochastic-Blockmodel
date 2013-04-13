// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Poisson"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Poisson", typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class PoissonOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sample, double mean)
		{
			return SampleAverageConditional(mean).GetLogProb(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sample, double mean) { return LogAverageFactor(sample, mean); }
		public static double AverageLogFactor(int sample, double mean) { return LogAverageFactor(sample, mean); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Poisson sample, [Fresh] Poisson to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Poisson sample) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_mean">Outgoing message to 'mean'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean) p(mean) factor(sample,mean))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double LogAverageFactor(int sample, [SkipIfUniform] Gamma mean, [Fresh] Gamma to_mean)
		{
			return to_mean.GetLogAverageOf(mean);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="to_mean">Outgoing message to 'mean'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean) p(mean) factor(sample,mean))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sample, Gamma mean, [Fresh] Gamma to_mean) { return LogAverageFactor(sample, mean, to_mean); }

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static Gamma MeanAverageConditional(int sample)
		{
			// p(mean) = mean^sample exp(-mean)/Gamma(sample+1)
			return new Gamma(sample + 1, 1);
		}

		const string NotSupportedMessage= "A Poisson factor with unobserved output is not yet implemented for Expectation Propagation.  Try using Variational Message Passing.";

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(sample) p(sample) factor(sample,mean)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Gamma MeanAverageConditional([SkipIfUniform] Poisson sample)
		{
			// Z = sum_sample sum_mean (mean*rate)^sample exp(-mean)/Gamma(sample+1)^(nu+1)/Z(rate,nu)
			if (sample.IsUniform()) return Gamma.Uniform();
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Poisson SampleAverageConditional(double mean)
		{
			return new Poisson(mean);
		}
		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(mean) p(mean) factor(sample,mean)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static Poisson SampleAverageConditional([SkipIfUniform] Gamma mean)
		{
			return new Poisson(mean.GetMean());
		}

		//-- VMP -------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean) p(mean) log(factor(sample,mean))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(int sample, [Proper] Gamma mean)
		{
			return sample * mean.GetMeanLog() - MMath.GammaLn(sample+1) - mean.GetMean();
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,mean) p(sample,mean) log(factor(sample,mean))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Poisson sample, [Proper] Gamma mean)
		{
			return sample.GetMean() * mean.GetMeanLog() - sample.GetMeanLogFactorial() - mean.GetMean();
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,mean))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Poisson sample, double mean)
		{
			return sample.GetMean() * mean - sample.GetMeanLogFactorial() - mean;
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static Gamma MeanAverageLogarithm(int sample)
		{
			return MeanAverageConditional(sample);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,mean)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static Gamma MeanAverageLogarithm([Proper] Poisson sample)
		{
			// p(mean) = exp(E[sample]*log(mean) - mean - E[log(sample!)])
			return new Gamma(sample.GetMean()+1, 1);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(sample,mean)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static Poisson SampleAverageLogarithm([Proper] Gamma mean)
		{
			return new Poisson(Math.Exp(mean.GetMeanLog()));
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Poisson SampleAverageLogarithm(double mean)
		{
			return new Poisson(Math.Exp(mean));
		}
	}
}
