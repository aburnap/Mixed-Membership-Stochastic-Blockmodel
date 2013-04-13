// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.GammaFromShapeAndRate"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "GammaFromShapeAndRate")]
	[Quality(QualityBand.Mature)]
	public static class GammaFromShapeAndRateOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(y,shape,rate))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double y, double shape, double rate)
		{
			return Gamma.GetLogProb(y, shape, rate);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double y, double shape, double rate) { return LogAverageFactor(y, shape, rate); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'.</param>
		/// <param name="to_y">Outgoing message to 'y'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y) p(y) factor(y,shape,rate))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma y, [Fresh] Gamma to_y)
		{
			return to_y.GetLogAverageOf(y);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y) p(y) factor(y,shape,rate) / sum_y p(y) messageTo(y))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gamma y, double shape, double rate)
		{
			return 0.0;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(rate) p(rate) factor(y,shape,rate))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double y, double shape, Gamma rate)
		{
			// f(y,a,b) = y^(a-1) b^a/Gamma(a) exp(-by) = y^(-2) Gamma(a+1)/Gamma(a) Ga(b; a+1, y)
			Gamma to_rate = RateAverageConditional(y, shape);
			return rate.GetLogAverageOf(to_rate) - 2*Math.Log(y) + Math.Log(shape);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(rate) p(rate) factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double y, double shape, Gamma rate)
		{
			return LogAverageFactor(y, shape, rate);
		}

		/// <summary>
		/// EP message to 'y'
		/// </summary>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>The outgoing EP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'y' conditioned on the given values.
		/// </para></remarks>
		public static Gamma YAverageConditional(double shape, double rate)
		{
			return Gamma.FromShapeAndRate(shape, rate);
		}

		/// <summary>
		/// EP message to 'rate'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <returns>The outgoing EP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'rate' conditioned on the given values.
		/// </para></remarks>
		public static Gamma RateAverageConditional(double y, double shape)
		{
			return Gamma.FromShapeAndRate(shape+1, y);
		}

		const string NotSupportedMessage = "Expectation Propagation does not support Gamma variables with stochastic shape or rate";

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y,shape,rate) p(y,shape,rate) factor(y,shape,rate))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y,shape) p(y,shape) factor(y,shape,rate))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape, double rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y,rate) p(y,rate) factor(y,shape,rate))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor([SkipIfUniform] Gamma y, double shape, [SkipIfUniform] Gamma rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(shape,rate) p(shape,rate) factor(y,shape,rate))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor(double y, [SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(shape) p(shape) factor(y,shape,rate))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static double LogAverageFactor(double y, [SkipIfUniform] Gamma shape, double rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'y'
		/// </summary>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'y' as the random arguments are varied.
		/// The formula is <c>proj[p(y) sum_(shape,rate) p(shape,rate) factor(y,shape,rate)]/p(y)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma YAverageConditional([SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		[Skip]
		public static Gamma YAverageConditionalInit()
		{
			return Gamma.Uniform();
		}

		/// <summary>
		/// EP message to 'y'
		/// </summary>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>The outgoing EP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'y' as the random arguments are varied.
		/// The formula is <c>proj[p(y) sum_(shape) p(shape) factor(y,shape,rate)]/p(y)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma YAverageConditional([SkipIfUniform] Gamma shape, double rate)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'y'
		/// </summary>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'y' as the random arguments are varied.
		/// The formula is <c>proj[p(y) sum_(rate) p(rate) factor(y,shape,rate)]/p(y)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		//[NotSupported(NotSupportedMessage)]
		[Quality(QualityBand.Experimental)]
		public static Gamma YAverageConditional(Gamma y, double shape, [SkipIfUniform] Gamma rate)
		{
			if (rate.IsPointMass) return YAverageConditional(shape, rate.Point);
			// q(b) = b^(a2-1) exp(-b2 b)
			// f(y,b) = y^(a-1) b^a exp(-by)
			// q(b) f(y,b) = y^(a-1) b^(a+a2-1) exp(-b (y+b2))
			// int over b = y^(a-1) / (y + b2)^(a+a2) = ratio of Ga(a,1)/Ga(a2,1/b2)
			double max = shape/rate.GetMean()*10;
			int n = 1000;
			double inc = max/n;
			GammaEstimator est = new GammaEstimator();
			for (int i = 0; i < n; i++) {
				double x = (i+1)*inc;
				double logp = (shape + y.Shape - 2)*Math.Log(x) - x*y.Rate - (shape + rate.Shape)*Math.Log(x + rate.Rate);
				est.Add(x, Math.Exp(logp));
			}
			Gamma post = est.GetDistribution(new Gamma());
			return post / rate;
			//throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'rate'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'rate' as the random arguments are varied.
		/// The formula is <c>proj[p(rate) sum_(y,shape) p(y,shape) factor(y,shape,rate)]/p(rate)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma RateAverageConditional([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'rate'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'rate' as the random arguments are varied.
		/// The formula is <c>proj[p(rate) sum_(shape) p(shape) factor(y,shape,rate)]/p(rate)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma RateAverageConditional(double y, [SkipIfUniform] Gamma shape)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'rate'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <returns>The outgoing EP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'rate' as the random arguments are varied.
		/// The formula is <c>proj[p(rate) sum_(y) p(y) factor(y,shape,rate)]/p(rate)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		//[NotSupported(NotSupportedMessage)]
		[Quality(QualityBand.Experimental)]
		public static Gamma RateAverageConditional([SkipIfUniform] Gamma y, double shape, Gamma rate)
		{
			if (y.IsPointMass) return RateAverageConditional(y.Point, shape);
			// q(y) = y^(a1-1) exp(-b1 y)
			// q(b) = b^(a2-1) exp(-b2 b)
			// f(y,b) = y^(a-1) b^a exp(-by)
			// q(y) q(b) f(y,b) = y^(a+a1-2) b^(a+a2-1) exp(-b2 b - b1 y - b y)
			// int over y = b^(a+a2-1) exp(-b2 b) / (b + b1)^(a+a1-1)
			// this is dist of ratio Ga(a+1,1)/Ga(a1-2,1) times Ga(a2,b2) pdf
			double max = shape/y.GetMean()*10;
			int n = 1000;
			double inc = max/n;
			GammaEstimator est = new GammaEstimator();
			for (int i = 0; i < n; i++) {
				double b = (i+1)*inc;
				double logp = (shape + rate.Shape - 1)*Math.Log(b) - b*rate.Rate - (shape + y.Shape - 1)*Math.Log(b + y.Rate);
				est.Add(b, Math.Exp(logp));
			}
			Gamma post = est.GetDistribution(new Gamma());
			//Console.WriteLine("y = {0}, shape = {1}, rate = {2}, post = {3}", y, shape, rate, post);
			return post / rate;
			//throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'shape'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'shape' as the random arguments are varied.
		/// The formula is <c>proj[p(shape) sum_(y) p(y) factor(y,shape,rate)]/p(shape)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma ShapeAverageConditional([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape, double rate, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'shape'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'shape' as the random arguments are varied.
		/// The formula is <c>proj[p(shape) sum_(rate) p(rate) factor(y,shape,rate)]/p(shape)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma ShapeAverageConditional(double y, [SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'shape'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'shape' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma ShapeAverageConditional(double y, [SkipIfUniform] Gamma shape, double rate, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'shape'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'shape' as the random arguments are varied.
		/// The formula is <c>proj[p(shape) sum_(y,rate) p(y,rate) factor(y,shape,rate)]/p(shape)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gamma ShapeAverageConditional([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate, [SkipIfUniform] Gamma result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		// VMP ----------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(double y, double shape, double rate)
		{
			return LogAverageFactor(y, shape, rate);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(y,shape,rate) p(y,shape,rate) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor([Proper] Gamma y, [Proper] Gamma shape, [Proper] Gamma rate)
		{
			return (shape.GetMean()-1)*y.GetMeanLog()+shape.GetMeanLog()*rate.GetMeanLog()
				-rate.GetMean()*y.GetMean()-ELogGamma(shape);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(y,shape) p(y,shape) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor([Proper] Gamma y, [Proper] Gamma shape, double rate)
		{
			return AverageLogFactor(y, shape, Gamma.PointMass(rate));
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(y,rate) p(y,rate) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] Gamma y, double shape, [SkipIfUniform] Gamma rate)
		{
			return AverageLogFactor(y, Gamma.PointMass(shape), rate);
		}

        /// <summary>
        /// Evidence message for VMP
        /// </summary>
        /// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="rate">Constant value for 'rate'.</param>
        /// <returns>Average of the factor's log-value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>sum_(y) p(y) log(factor(y,shape,rate))</c>.
        /// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
        /// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
        public static double AverageLogFactor([Proper] Gamma y, double shape, double	rate)
		{
			Gamma to_y = YAverageLogarithm(shape, rate);
			return y.GetAverageLog(to_y);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(shape,rate) p(shape,rate) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor(double y, [SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate)
		{
			return AverageLogFactor(Gamma.PointMass(y), shape, rate);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(shape) p(shape) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static double AverageLogFactor(double y, [SkipIfUniform] Gamma shape, double rate)
		{
			return AverageLogFactor(Gamma.PointMass(y), shape, Gamma.PointMass(rate));
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(rate) p(rate) log(factor(y,shape,rate))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		public static double AverageLogFactor(double y, double shape, [SkipIfUniform] Gamma rate)
		{
			return AverageLogFactor(Gamma.PointMass(y), Gamma.PointMass(shape), rate);
		}

		/// <summary>
		/// VMP message to 'y'
		/// </summary>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'y'.
		/// The formula is <c>exp(sum_(shape,rate) p(shape,rate) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma YAverageLogarithm([SkipIfUniform] Gamma shape, [SkipIfUniform] Gamma rate)
		{
			return Gamma.FromShapeAndRate(shape.GetMean(), rate.GetMean());
		}

		/// <summary>
		/// VMP message to 'y'
		/// </summary>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <returns>The outgoing VMP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'y'.
		/// The formula is <c>exp(sum_(shape) p(shape) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma YAverageLogarithm([SkipIfUniform] Gamma shape, double rate)
		{
			return Gamma.FromShapeAndRate(shape.GetMean(), rate);
		}

		/// <summary>
		/// VMP message to 'y'
		/// </summary>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'y'.
		/// The formula is <c>exp(sum_(rate) p(rate) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		public static Gamma YAverageLogarithm(double shape, [SkipIfUniform] Gamma rate)
		{
			return Gamma.FromShapeAndRate(shape, rate.GetMean());
		}

        /// <summary>
        /// VMP message to 'y'
        /// </summary>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="rate">Constant value for 'rate'..</param>
        /// <returns>The outgoing VMP message to the 'y' argument</returns>
        /// <remarks><para>
        /// The message is simply the factor itself.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		public static Gamma YAverageLogarithm(double shape, double rate)
		{
			return Gamma.FromShapeAndRate(shape, rate);
		}

		/// <summary>
		/// VMP message to 'rate'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'rate'.
		/// The formula is <c>exp(sum_(y,shape) p(y,shape) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		public static Gamma RateAverageLogarithm([SkipIfUniform] Gamma y, [SkipIfUniform] Gamma shape)
		{
			// factor = rate^shape exp(-y*rate)
			// log(factor) = shape*log(rate) - y*rate
			// E[log(factor)] = E[shape]*log(rate) - E[y]*rate
			return Gamma.FromShapeAndRate(shape.GetMean() + 1, y.GetMean());
		}

		/// <summary>
		/// VMP message to 'rate'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'rate'.
		/// The formula is <c>exp(sum_(shape) p(shape) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		public static Gamma RateAverageLogarithm(double y, [SkipIfUniform] Gamma shape)
		{
			return Gamma.FromShapeAndRate(shape.GetMean() + 1, y);
		}

		/// <summary>
		/// VMP message to 'rate'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Constant value for 'shape'.</param>
		/// <returns>The outgoing VMP message to the 'rate' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'rate'.
		/// The formula is <c>exp(sum_(y) p(y) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		public static Gamma RateAverageLogarithm([SkipIfUniform] Gamma y, double shape)
		{
			return Gamma.FromShapeAndRate(shape + 1, y.GetMean());
		}

        /// <summary>
        /// VMP message to 'rate'
        /// </summary>
        /// <param name="y">Constant value for 'y'..</param>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <returns>The outgoing VMP message to the 'rate' argument</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'rate'.
        /// The formula is <c>exp(sum_(y) p(y) log(factor(y,shape,rate)))</c>.
        /// </para></remarks>
		public static Gamma RateAverageLogarithm(double y, double shape)
		{
			return Gamma.FromShapeAndRate(shape + 1, y);
		}

		/// <summary>
		/// VMP message to 'shape'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <param name="to_shape">Previous outgoing message to 'shape'.</param>
		/// <returns>The outgoing VMP message to the 'shape' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'shape'.
		/// The formula is <c>exp(sum_(y) p(y) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma ShapeAverageLogarithm([SkipIfUniform] Gamma y, [Proper] Gamma shape, double rate, Gamma to_shape)
		{
			return ShapeAverageLogarithm(y, shape, Gamma.PointMass(rate), to_shape);
		}

		/// <summary>
		/// VMP message to 'shape'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_shape">Previous outgoing message to 'shape'.</param>
		/// <returns>The outgoing VMP message to the 'shape' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'shape'.
		/// The formula is <c>exp(sum_(rate) p(rate) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma ShapeAverageLogarithm(double y, [Proper] Gamma shape, [Proper] Gamma rate, Gamma to_shape)
		{
			return ShapeAverageLogarithm(Gamma.PointMass(y), shape, rate, to_shape);
		}

		/// <summary>
		/// VMP message to 'shape'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Constant value for 'rate'.</param>
		/// <param name="to_shape">Previous outgoing message to 'shape'.</param>
		/// <returns>The outgoing VMP message to the 'shape' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'shape' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma ShapeAverageLogarithm(double y, [Proper] Gamma shape, double rate, Gamma to_shape)
		{
			return ShapeAverageLogarithm(Gamma.PointMass(y), shape, Gamma.PointMass(rate), to_shape);
		}

		/// <summary>
		/// VMP message to 'shape'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="shape">Incoming message from 'shape'. Must be a proper distribution.  If uniform, the result will be uniform. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="rate">Incoming message from 'rate'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_shape">Previous outgoing message to 'shape'.</param>
		/// <returns>The outgoing VMP message to the 'shape' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'shape'.
		/// The formula is <c>exp(sum_(y,rate) p(y,rate) log(factor(y,shape,rate)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="shape"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="rate"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gamma ShapeAverageLogarithm([SkipIfUniform] Gamma y, [Proper] Gamma shape, [Proper] Gamma rate, Gamma to_shape)
		{
			//if (to_shape.IsUniform())
			//	to_shape.SetShapeAndRate(1, 1);
			//if (y.Rate == 0 || Double.IsInfinity(y.Rate))
			//	y.SetShapeAndRate(1, 1);
			//if (rate.Rate==0 || Double.IsInfinity(rate.Rate))
			//	rate.SetShapeAndRate(1, 1);
			double ElogYR = y.GetMeanLog() + rate.GetMeanLog();
			double a, b;
			a = shape.Shape;
			b = shape.Rate;

			// Find required expectations using quadrature
			Vector gradElogGamma = CalculateDerivatives(shape);

			// Calculate gradients
			Vector gradS = -gradElogGamma;
			gradS[0] += ElogYR / b;
			gradS[1] += -a * ElogYR / (b * b);

			// Calculate the required message to match these gradients
			Gamma approximateFactor = NonconjugateProjection(shape, gradS);

			double damping = 0.0;
			if (damping == 0.0)
				return approximateFactor;
			else
				return (approximateFactor^(1-damping)) * (to_shape ^ damping);
		}

        /// <summary>
        /// Project the standard VMP message onto a gradient matched Gamma message. 
        /// </summary>
        /// <param name="context">Incoming message.</param>
        /// <param name="gradS">Gradient of S=int factor(x) p(x) dx</param>
        /// <returns>Projected gamma message</returns>
		internal static Gamma NonconjugateProjection(Gamma context, Vector gradS)
		{
			Matrix mat = new Matrix(2, 2);
			mat[0, 0] = MMath.Trigamma(context.Shape);
			mat[1, 0] = mat[0, 1] = -1 / context.Rate;
			mat[1, 1] = context.Shape / (context.Rate * context.Rate);
			Vector v = twoByTwoInverse(mat) * gradS;
			return Gamma.FromShapeAndRate(v[0] + 1, v[1]);
		}

        /// <summary>
        /// Two by two matrix inversion. 
        /// </summary>
        /// <param name="a">Matrix to invert</param>
        /// <returns>Inverted matrix</returns>
		internal static Matrix twoByTwoInverse(Matrix a)
		{
			Matrix result = new Matrix(2, 2);
			double det = a[0, 0] * a[1, 1] - a[0, 1] * a[1, 0];
			result[0, 0] = a[1, 1] / det;
			result[0, 1] = -a[0, 1] / det;
			result[1, 0] = -a[1, 0] / det;
			result[1, 1] = a[0, 0] / det;
			return result;
		}

		/// <summary>
		/// Calculate derivatives of \int G(x;a,b) LogGamma(x) dx wrt (a,b)
		/// </summary>
		/// <param name="q">The Gamma distribution G(x;a,b).</param>
		/// <returns>A 2-vector containing derivatives of \int G(x;a,b) LogGamma(x) dx wrt (a,b).</returns>
		/// <remarks><para>
		///  Calculates expectations in x=log(s) space using Gauss-Hermite quadrature. 
		///  For each integral the behaviour as x->0 is subtracted from the integrand 
		///  before performing quadrature to remove the singularity there. 
		/// </para></remarks>
		public static Vector CalculateDerivatives(Gamma q)
		{
			Vector gradElogGamma = Vector.Zero(2);
			// Get shape and scale of the distribution
			double a = q.Shape;
			double b = q.Rate;
			// Mean in the transformed domain
			double proposalMean = q.GetMeanLog();
			// Laplace approximation of variance in transformed domain 
			double proposalVariance = 1 / a;
            //double proposalVariance = Math.Exp(-proposalMean) / b; 

			// Quadrature coefficient
			int n = 11;
			Vector nodes = Vector.Zero(n);
			Vector weights = Vector.Zero(n);
			Quadrature.GaussianNodesAndWeights(proposalMean, proposalVariance, nodes, weights);

			double EXDigamma = 0;
			double ELogGam = 0;
			double ELogXLogGamma = 0;

			// Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
			double logZ = MMath.GammaLn(a) - a * Math.Log(b);
			for (int i = 0; i < n; i++) {
				double x = nodes[i];
				double expx = Math.Exp(x);
				double p = a * x - b * expx - logZ - Gaussian.GetLogProb(x, proposalMean, proposalVariance);
				p = Math.Exp(p) * weights[i];
				EXDigamma += p * (expx * MMath.Digamma(expx) + 1);
				ELogGam += p * (MMath.GammaLn(expx) + x);
				ELogXLogGamma += p * (x * MMath.GammaLn(expx) + x * x);
			}

			// Normalise and add removed components
			ELogGam = ELogGam - proposalMean;
			ELogXLogGamma = ELogXLogGamma - (MMath.Trigamma(a) + proposalMean * proposalMean);
			EXDigamma = EXDigamma - 1;

			// Calculate derivatives
			gradElogGamma[0] = ELogXLogGamma - proposalMean * ELogGam;
			gradElogGamma[1] = -1.0 / b * EXDigamma;
			//gradElogGamma[1] = (ELogGamma(q) - ELogGamma(new Gamma(a + 1, b))) * a / b; 
			return gradElogGamma;
		}

		public static Vector CalculateDerivativesTrapezoid(Gamma q)
		{
			Vector gradElogGamma = Vector.Zero(2);
			// Get shape and scale of the distribution
			double a = q.Shape;
			double b = q.Rate;
			double mean, variance;
			q.GetMeanAndVariance(out mean, out variance);
			double upperBound = 10;

			int n = 10000;
			double ELogGamma = 0, ELogXLogGamma = 0, ExDigamma = 0;
			double inc = upperBound/n;
			for (int i = 0; i < n; i++) {
				double x = inc * (i+1);
				double logp = q.GetLogProb(x);
				double p = Math.Exp(logp);
				double f = p * MMath.GammaLn(x);
				ELogGamma += f;
				ELogXLogGamma += Math.Log(x)*f;
				ExDigamma += x*MMath.Digamma(x)*p;
			}
			ELogGamma *= inc;
			ELogXLogGamma *= inc;
			ExDigamma *= inc;
			gradElogGamma[0] = ELogXLogGamma + (Math.Log(b) - MMath.Digamma(a))*ELogGamma;
			gradElogGamma[1] = -ExDigamma/b;
			return gradElogGamma;
		}

        public static Vector CalculateDerivativesNaive(Gamma q)
        {
            Vector gradElogGamma = Vector.Zero(2);
            // Get shape and scale of the distribution
            double a = q.Shape;
            double b = q.Rate;
            // Mean in the transformed domain
            double proposalMean = q.GetMeanLog();
            // Laplace approximation of variance in transformed domain 
            double proposalVariance = 1 / a;

            // Quadrature coefficient
            int n = 11;
            Vector nodes = Vector.Zero(n);
            Vector weights = Vector.Zero(n);
            Quadrature.GaussianNodesAndWeights(proposalMean, proposalVariance, nodes, weights);

            double EXDigamma = 0;
            double ELogGam = 0;
            double ELogXLogGamma = 0;

            // Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
            double logZ = MMath.GammaLn(a) - a * Math.Log(b);
            for (int i = 0; i < n; i++)
            {
                double x = nodes[i];
                double expx = Math.Exp(x);
                double p = a * x - b * expx - logZ - Gaussian.GetLogProb(x, proposalMean, proposalVariance);
                p = Math.Exp(p) * weights[i];
                EXDigamma += p * (expx * MMath.Digamma(expx));
                ELogGam += p * (MMath.GammaLn(expx));
                ELogXLogGamma += p * (x * MMath.GammaLn(expx));
            }

            // Calculate derivatives
            gradElogGamma[0] = ELogXLogGamma - proposalMean * ELogGam;
            gradElogGamma[1] = -1.0 / b * EXDigamma;
            //gradElogGamma[1] = (ELogGamma(q) - ELogGamma(new Gamma(a + 1, b))) * a / b; 
            return gradElogGamma;
        }

		/// <summary>
		/// Calculates \int G(x;a,b) LogGamma(x) dx
		/// </summary>
		/// <param name="q">G(x;a,b)</param>
		/// <returns>\int G(x;a,b) LogGamma(x) dx</returns>
		public static double ELogGamma(Gamma q)
		{
			if (q.IsPointMass)
				return MMath.GammaLn(q.Point);
			double a = q.Shape;
			double b = q.Rate;
			// Mean in the transformed domain
			double proposalMean = q.GetMeanLog();
			// Laplace approximation of variance in transformed domain 
			double proposalVariance = 1 / a;

			// Quadrature coefficient
			int n = 11;
			Vector nodes = Vector.Zero(n);
			Vector weights = Vector.Zero(n);
			Quadrature.GaussianNodesAndWeights(proposalMean, proposalVariance, nodes, weights);

			double ELogGamma = 0;
			double logZ = -a * Math.Log(b) + MMath.GammaLn(a);
			// Calculate expectations in x=log(s) space using Gauss-Hermite quadrature
			for (int i = 0; i < n; i++) {
				double x = nodes[i];
				double expx = Math.Exp(x);
				double p = a * x - b * expx - Gaussian.GetLogProb(x, proposalMean, proposalVariance) - logZ;
				p = Math.Exp(p + Math.Log(weights[i]));
				ELogGamma += p * (MMath.GammaLn(expx) + x);
			}

			// Add removed components
			return ELogGamma - proposalMean;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Gamma.SampleFromMeanAndVariance"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "sample", "mean", "variance" }, typeof(Gamma), "SampleFromMeanAndVariance")]
	[Quality(QualityBand.Stable)]
	public static class GammaFromMeanAndVarianceOp
	{
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
			return Gamma.FromShapeAndScale(mean, variance).GetLogProb(sample);
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
		public static double LogEvidenceRatio(double sample, double mean, double variance) { return LogAverageFactor(sample, mean, variance); }
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
		public static double AverageLogFactor(double sample, double mean, double variance) { return LogAverageFactor(sample, mean, variance); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma sample, [Fresh] Gamma to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
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
		public static Gamma SampleAverageLogarithm(double mean, double variance)
		{
			return Gamma.FromMeanAndVariance(mean, variance);
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
		public static Gamma SampleAverageConditional(double mean, double variance)
		{
			return Gamma.FromMeanAndVariance(mean, variance);
		}
	}

	/// <summary>
	/// Provides just the operators for constructing Gamma distributions with this parameterisation.
	/// </summary>
	[FactorMethod(typeof(Gamma), "Sample", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class GammaFromShapeAndScaleOp
	{
        /// <summary>
        /// Evidence message for EP
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="scale">Constant value for 'scale'.</param>
        /// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
        /// </para></remarks>
		public static double LogAverageFactor(double sample, double shape, double scale)
		{
			return Gamma.FromShapeAndScale(shape, scale).GetLogProb(sample);
		}

        /// <summary>
        /// Evidence message for EP
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="scale">Constant value for 'scale'.</param>
        /// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
        /// </para></remarks>
		public static double LogEvidenceRatio(double sample, double shape, double scale) { return LogAverageFactor(sample, shape, scale); }

		[Skip]
		public static double LogEvidenceRatio(Gamma sample, double shape, double scale) { return 0.0; }

        /// <summary>
        /// Evidence message for VMP
        /// </summary>
        /// <param name="sample">Constant value for 'sample'.</param>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="scale">Constant value for 'scale'.</param>
        /// <returns>Average of the factor's log-value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
        /// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
        /// </para></remarks>
		public static double AverageLogFactor(double sample, double shape, double scale) { return LogAverageFactor(sample, shape, scale); }

		public static double AverageLogFactor([Proper] Gamma sample, double shape, double scale)
		{
			return GammaFromShapeAndRateOp.AverageLogFactor(sample, shape, 1/scale);
		}

        /// <summary>
        /// Evidence message for EP
        /// </summary>
        /// <param name="sample">Incoming message to 'sample'.</param>
        /// <param name="to_sample">Message sent to 'sample'</param>
        /// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
        /// <remarks><para>
        /// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
        /// </para></remarks>
		public static double LogAverageFactor(Gamma sample, [Fresh] Gamma to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

        /// <summary>
        /// VMP message to 'sample'
        /// </summary>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="scale">Constant value for 'scale'.</param>
        /// <returns>The outgoing VMP message to the 'sample' argument</returns>
        /// <remarks><para>
        /// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
        /// </para></remarks>
		public static Gamma SampleAverageLogarithm(double shape, double scale)
		{
			return Gamma.FromShapeAndScale(shape, scale);
		}

        /// <summary>
        /// EP message to 'sample'
        /// </summary>
        /// <param name="shape">Constant value for 'shape'.</param>
        /// <param name="scale">Constant value for 'scale'.</param>
        /// <returns>The outgoing EP message to the 'sample' argument</returns>
        /// <remarks><para>
        /// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
        /// </para></remarks>
		public static Gamma SampleAverageConditional(double shape, double scale)
		{
			return Gamma.FromShapeAndScale(shape, scale);
		}
	}
}
