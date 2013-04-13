// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.BetaFromMeanAndTotalCount"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "BetaFromMeanAndTotalCount")]
	[Quality(QualityBand.Experimental)]
	public static class BetaOp
	{
		/// <summary>
		/// How much damping to use to avoid improper messages. A higher value implies more damping. 
		/// </summary>
		public static double damping = 0.8;

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (prob,mean,totalCount)</c>.
		/// </para></remarks>
		public static double AverageLogFactor(double prob, double mean, double totalCount)
		{
			return LogAverageFactor(prob, mean, totalCount);
		}

		// TODO: VMP evidence messages for stochastic inputs (see DirichletOp)

		/// <summary>
		/// VMP message to 'Prob'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'Prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Y'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (mean,totalCount)</c>.
		/// </para></remarks>
		public static Beta ProbAverageLogarithm([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount)
		{
			double meanMean = mean.GetMean();
			double totalCountMean = totalCount.GetMean();
			return (new Beta(meanMean * totalCountMean, (1 - meanMean) * totalCountMean));
		}

		/// <summary>
		/// VMP message to 'Prob'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'Prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Y'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (mean,totalCount)</c>.
		/// </para></remarks>
		public static Beta ProbAverageLogarithm([SkipIfUniform] Beta mean, double totalCount)
		{
			double meanMean = mean.GetMean();
			return (new Beta(meanMean * totalCount, (1 - meanMean) * totalCount));
		}

		/// <summary>
		/// VMP message to 'Prob'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'.</param>
		/// <returns>The outgoing VMP message to the 'Prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'Y'.
		/// The formula is <c>int log(f(sample,x)) q(x) dx</c> where <c>x = (mean,totalCount)</c>.
		/// </para></remarks>
		public static Beta ProbAverageLogarithm(double mean, [SkipIfUniform] Gamma totalCount)
		{
			double totalCountMean = totalCount.GetMean();
			return (new Beta(mean * totalCountMean, (1 - mean) * totalCountMean));
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_mean">Previous outgoing message to 'Mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(totalCount) p(totalCount) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Beta MeanAverageLogarithm([Proper] Beta mean, [Proper] Gamma totalCount, double prob, Beta to_mean)
		{
			return MeanAverageLogarithm(mean, totalCount, Beta.PointMass(prob), to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_mean">Previous outgoing message to 'Mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static Beta MeanAverageLogarithm([Proper] Beta mean, double totalCount, double prob, Beta to_mean)
		{
			return MeanAverageLogarithm(mean, Gamma.PointMass(totalCount), Beta.PointMass(prob), to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_mean">Previous outgoing message to 'Mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(prob) p(prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Beta MeanAverageLogarithm([Proper] Beta mean, double totalCount, [SkipIfUniform] Beta prob, Beta to_mean)
		{
			return MeanAverageLogarithm(mean, Gamma.PointMass(totalCount), prob, to_mean);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_mean">Previous outgoing message to 'Mean'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(totalCount,prob) p(totalCount,prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Beta MeanAverageLogarithm([Proper] Beta mean, [Proper] Gamma totalCount, [SkipIfUniform] Beta prob, Beta to_mean)
		{
			// Calculate gradient using method for DirichletOp
			double ELogP, ELogOneMinusP;
			prob.GetMeanLogs(out ELogP, out ELogOneMinusP);
			Vector gradS = DirichletOp.CalculateGradientForMean(
				 Vector.FromArray(new double[] { mean.TrueCount, mean.FalseCount }),
				 totalCount,
				 Vector.FromArray(new double[] { ELogP, ELogOneMinusP }));
			// Project onto a Beta distribution 
			Matrix A = new Matrix(2, 2);
			double c = MMath.Trigamma(mean.TotalCount);
			A[0, 0] = MMath.Trigamma(mean.TrueCount) - c;
			A[1, 0] = A[0, 1] = -c;
			A[1, 1] = MMath.Trigamma(mean.FalseCount) - c;
			Vector theta = GammaFromShapeAndRateOp.twoByTwoInverse(A)*gradS;
			Beta approximateFactor = new Beta(theta[0] + 1, theta[1] + 1);
			if (damping == 0.0)
				return approximateFactor;
			else
				return (approximateFactor^(1-damping)) * (to_mean ^ damping);
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'TotalCount'.</param>
		/// <returns>The outgoing VMP message to the 'totalCount' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'totalCount'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm([Proper] Beta mean, [Proper] Gamma totalCount, double prob, Gamma to_totalCount)
		{
			return TotalCountAverageLogarithm(mean, totalCount, Beta.PointMass(prob), to_totalCount);
		}

		/// <summary>
		/// VMP message to 'totalCount'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_totalCount">Previous outgoing message to 'TotalCount'.</param>
		/// <returns>The outgoing VMP message to the 'totalCount' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'totalCount'.
		/// The formula is <c>exp(sum_(mean,prob) p(mean,prob) log(factor(prob,mean,totalCount)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		public static Gamma TotalCountAverageLogarithm([Proper] Beta mean, [Proper] Gamma totalCount, [SkipIfUniform] Beta prob, Gamma to_totalCount)
		{
			double ELogP, ELogOneMinusP;
			prob.GetMeanLogs(out ELogP, out ELogOneMinusP);
			Gamma approximateFactor = DirichletOp.TotalCountAverageLogarithmHelper(
				Vector.FromArray(new double[] { mean.TrueCount, mean.FalseCount }),
				totalCount,
				Vector.FromArray(new double[] { ELogP, ELogOneMinusP }));
			if (damping == 0.0)
				return approximateFactor;
			else
				return (approximateFactor^(1-damping)) * (to_totalCount ^ damping);
		}

		//---------------------------- EP -----------------------------

		const string NotSupportedMessage = "Expectation Propagation does not currently support beta distributions with stochastic arguments.";

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_() p() factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double prob, double mean, double totalCount)
		{
			var g = new Beta(mean * totalCount, (1 - mean) * totalCount);
			return g.GetLogProb(prob);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="prob">Incoming message from 'prob'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(prob) p(prob) factor(prob,mean,totalCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Beta prob, double mean, double totalCount)
		{
			var g = new Beta(mean*totalCount, (1-mean)*totalCount);
			return g.GetLogAverageOf(prob);
		}

		/// <summary>
		/// EP message to 'prob'.
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing EP message to the 'prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'prob' conditioned on the given values.
		/// </para></remarks>
		public static Beta ProbAverageConditional(double mean, double totalCount)
		{
			return new Beta(mean * totalCount, (1 - mean) * totalCount);
		}

		/// <summary>
		/// EP message to 'prob'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(mean,totalCount) p(mean,totalCount) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta ProbAverageConditional([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <returns>The outgoing EP message to the 'prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(mean) p(mean) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta ProbAverageConditional([SkipIfUniform] Beta mean, double totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'prob'.
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'prob' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'prob' as the random arguments are varied.
		/// The formula is <c>proj[p(prob) sum_(totalCount) p(totalCount) factor(prob,mean,totalCount)]/p(prob)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta ProbAverageConditional(double mean, [SkipIfUniform] Gamma totalCount)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta MeanAverageConditional([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount, double prob, [SkipIfUniform] Beta result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta MeanAverageConditional([SkipIfUniform] Beta mean, double totalCount, double prob, [SkipIfUniform] Beta result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Constant value for 'totalCount'.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(prob) p(prob) factor(prob,mean,totalCount)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta MeanAverageConditional([SkipIfUniform] Beta mean, double totalCount, [SkipIfUniform] Beta prob, [SkipIfUniform] Beta result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'mean'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(totalCount,prob) p(totalCount,prob) factor(prob,mean,totalCount)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Beta MeanAverageConditional([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount, [SkipIfUniform] Beta prob, [SkipIfUniform] Beta result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Constant value for 'prob'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'totalCount' as the random arguments are varied.
		/// The formula is <c>proj[p(totalCount) sum_(mean) p(mean) factor(prob,mean,totalCount)]/p(totalCount)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount, double prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'totalCount'.
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="totalCount">Incoming message from 'totalCount'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="prob">Incoming message from 'prob'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'totalCount' as the random arguments are varied.
		/// The formula is <c>proj[p(totalCount) sum_(mean,prob) p(mean,prob) factor(prob,mean,totalCount)]/p(totalCount)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="totalCount"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="prob"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma TotalCountAverageConditional([SkipIfUniform] Beta mean, [SkipIfUniform] Gamma totalCount, [SkipIfUniform] Beta prob, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Beta.Sample(double,double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Beta), "Sample", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class BetaFromTrueAndFalseCountsOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,trueCount,falseCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Beta sample, double trueCount, double falseCount, [Fresh] Beta to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,trueCount,falseCount) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Beta sample, double trueCount, double falseCount) { return 0.0; }

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,trueCount,falseCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Beta sample, double trueCount, double falseCount, [Fresh] Beta to_sample)
		{
			return to_sample.GetAverageLog(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,trueCount,falseCount))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double sample, double trueCount, double falseCount)
		{
			return SampleAverageConditional(trueCount, falseCount).GetLogProb(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,trueCount,falseCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double sample, double trueCount, double falseCount)
		{
			return LogAverageFactor(sample, trueCount, falseCount);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,trueCount,falseCount))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(double sample, double trueCount, double falseCount)
		{
			return LogAverageFactor(sample, trueCount, falseCount);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Beta SampleAverageLogarithm(double trueCount, double falseCount)
		{
			return new Beta(trueCount, falseCount);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="trueCount">Constant value for 'trueCount'.</param>
		/// <param name="falseCount">Constant value for 'falseCount'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Beta SampleAverageConditional(double trueCount, double falseCount)
		{
			return new Beta(trueCount, falseCount);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Beta.SampleFromMeanAndVariance"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "sample", "mean", "variance" }, typeof(Beta), "SampleFromMeanAndVariance")]
	[Quality(QualityBand.Stable)]
	public static class BetaFromMeanAndVarianceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Beta sample, double mean, double variance, [Fresh] Beta to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Beta sample, double mean, double variance) { return 0.0; }

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,mean,variance))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Beta sample, double mean, double variance, [Fresh] Beta to_sample)
		{
			return to_sample.GetAverageLog(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double sample, double mean, double variance)
		{
			return SampleAverageConditional(mean, variance).GetLogProb(sample);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double sample, double mean, double variance)
		{
			return LogAverageFactor(sample, mean, variance);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(double sample, double mean, double variance)
		{
			return LogAverageFactor(sample, mean, variance);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Beta SampleAverageLogarithm(double mean, double variance)
		{
			return Beta.FromMeanAndVariance(mean, variance);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Beta SampleAverageConditional(double mean, double variance)
		{
			return Beta.FromMeanAndVariance(mean, variance);
		}
	}
}
