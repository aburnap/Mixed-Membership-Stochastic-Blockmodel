// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.IsPositive"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "IsPositive", typeof(double), Default=true)]
	[Quality(QualityBand.Mature)]
	public static class IsPositiveOp
	{
		public static double LogAverageFactor(bool isPositive, double x)
		{
			return (isPositive == Factor.IsPositive(x)) ? 0.0 : Double.NegativeInfinity;
		}
		public static double LogEvidenceRatio(bool isPositive, double x) { return LogAverageFactor(isPositive, x); }
		public static double AverageLogFactor(bool isPositive, double x) { return LogAverageFactor(isPositive, x); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(x) p(x) factor(isPositive,x))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool isPositive, Gaussian x)
		{
			return LogAverageFactor(Bernoulli.PointMass(isPositive), x);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isPositive">Incoming message from 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isPositive,x) p(isPositive,x) factor(isPositive,x))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli isPositive, Gaussian x)
		{
			Bernoulli to_isPositive = IsPositiveAverageConditional(x);
			return isPositive.GetLogAverageOf(to_isPositive);
#if false
			// Z = p(b=T) p(x > 0) + p(b=F) p(x <= 0)
			//   = p(b=F) + (p(b=T) - p(b=F)) p(x > 0)
			if (x.IsPointMass) {
				return Factor.IsPositive(x.Point) ? isPositive.GetLogProbTrue() : isPositive.GetLogProbFalse();
			} else if(x.IsUniform()) {
				return Bernoulli.LogProbEqual(isPositive.LogOdds,0.0);
			} else {
				// m/sqrt(v) = (m/v)/sqrt(1/v)
				double z = x.MeanTimesPrecision / Math.Sqrt(x.Precision);
				if (isPositive.IsPointMass) {
					return isPositive.Point ? MMath.NormalCdfLn(z) : MMath.NormalCdfLn(-z);
				} else {
					return MMath.LogSumExp(isPositive.GetLogProbTrue() + MMath.NormalCdfLn(z), isPositive.GetLogProbFalse() + MMath.NormalCdfLn(-z));
				}
			}
#endif
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isPositive">Incoming message from 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isPositive,x) p(isPositive,x) factor(isPositive,x) / sum_isPositive p(isPositive) messageTo(isPositive))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli isPositive, Gaussian x) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(x) p(x) factor(isPositive,x))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool isPositive, Gaussian x) { return LogAverageFactor(isPositive, x); }

		/// <summary>
		/// EP message to 'isPositive'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'isPositive' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isPositive' as the random arguments are varied.
		/// The formula is <c>proj[p(isPositive) sum_(x) p(x) factor(isPositive,x)]/p(isPositive)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Bernoulli IsPositiveAverageConditional([SkipIfUniform, Proper] Gaussian x)
		{
			Bernoulli result = new Bernoulli();
			if (x.IsPointMass) {
				result.Point = Factor.IsPositive(x.Point);
			} else if (x.IsUniform()) {
				result.LogOdds = 0.0;
			} else if (!x.IsProper()) {
				throw new ImproperMessageException(x);
			} else {
				// m/sqrt(v) = (m/v)/sqrt(1/v)
				double z = x.MeanTimesPrecision / Math.Sqrt(x.Precision);
				// p(true) = NormalCdf(z)
				// log(p(true)/p(false)) = NormalCdfLogit(z)
				result.LogOdds = MMath.NormalCdfLogit(z);
			}
			return result;
		}
		[Skip]
		public static Bernoulli IsPositiveAverageConditionalInit()
		{
			return new Bernoulli();
		}

		/// <summary>
		/// EP message to 'x'
		/// </summary>
		/// <param name="isPositive">Incoming message from 'isPositive'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'x' as the random arguments are varied.
		/// The formula is <c>proj[p(x) sum_(isPositive) p(isPositive) factor(isPositive,x)]/p(x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isPositive"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Gaussian XAverageConditional([SkipIfUniform] Bernoulli isPositive, [SkipIfUniform, Proper] Gaussian x)
		{
			return XAverageConditional_Helper(isPositive, x, false);
		}
		public static Gaussian XAverageConditional_Helper([SkipIfUniform] Bernoulli isPositive, [SkipIfUniform, Proper] Gaussian x, bool forceProper)
		{
			Gaussian result = new Gaussian();
			if (x.IsPointMass) {
				result.SetToUniform();
				return result;
			}
			double prec = x.Precision;
			if (prec == 0.0) {
				result.SetToUniform();
				return result;
			} else if (prec < 0) {
				throw new ImproperMessageException(x);
			}
			double sqrtPrec = Math.Sqrt(prec);
			double tau = x.MeanTimesPrecision;
			// m/sqrt(v) = (m/v)/sqrt(1/v)
			double z = tau / sqrtPrec;
			// epsilon = p(b=F)
			// eq (51) in EP quickref
			double alpha;
			if (isPositive.IsPointMass) {
				if (isPositive.Point) {
					if (z < -100) {
						double Y = MMath.NormalCdfRatio(z);
						double d2Y = 2*MMath.NormalCdfMomentRatio(2, z);
						double m2TimesPrec = d2Y/Y;
						Assert.IsTrue(tau != 0);
						double mp = (m2TimesPrec - 1)/tau;
						double vp = m2TimesPrec/prec - mp*mp;
						return (new Gaussian(mp, vp))/x;
					}
					alpha = sqrtPrec / MMath.NormalCdfRatio(z);
				} else {
					if (z > 100) {
						double Y = MMath.NormalCdfRatio(-z);
						// dY = -(d2Y/prec - Y)/(-z)*sqrtPrec
						// dY/Y/prec = -(d2Y/Y/prec/prec - 1/prec)/(-z)*sqrtPrec
						// dY/Y/prec = -(d2Y/Y/prec - 1)/(-tau)
						//double dY = -MMath.NormalCdfMomentRatio(1,-z)*sqrtPrec;
						double d2Y = 2*MMath.NormalCdfMomentRatio(2, -z);
						double m2TimesPrec = d2Y/Y;
						Assert.IsTrue(tau != 0);
						double mp = (m2TimesPrec - 1)/tau;
						double vp = m2TimesPrec/prec - mp*mp;
						return (new Gaussian(mp, vp))/x;
					}
					alpha = -sqrtPrec / MMath.NormalCdfRatio(-z);
				}
			} else {
				//double v = MMath.LogSumExp(isPositive.LogProbTrue + MMath.NormalCdfLn(z), isPositive.LogProbFalse + MMath.NormalCdfLn(-z));
				double v = LogAverageFactor(isPositive, x);
				alpha = sqrtPrec * Math.Exp(-z * z * 0.5 - MMath.LnSqrt2PI - v) * (2 * isPositive.GetProbTrue() - 1);
			}
			// eq (52) in EP quickref (where tau = mnoti/Vnoti)
			double beta = alpha * (alpha + tau);
			double weight = beta / (prec - beta);
			if (forceProper && weight < 0) weight = 0;
			// eq (31) in EP quickref; same as inv(inv(beta)-inv(prec))
			result.Precision = prec * weight;
			// eq (30) in EP quickref times above and simplified
			result.MeanTimesPrecision = weight * (tau + alpha) + alpha;
			if (double.IsNaN(result.Precision) || double.IsNaN(result.MeanTimesPrecision)) throw new ApplicationException("result is nan");
			return result;
		}

		/// <summary>
		/// EP message to 'x'
		/// </summary>
		/// <param name="isPositive">Incoming message from 'isPositive'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <returns>The outgoing EP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'x' as the random arguments are varied.
		/// The formula is <c>proj[p(x) sum_(isPositive) p(isPositive) factor(isPositive,x)]/p(x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isPositive"/> is not a proper distribution</exception>
		public static Gaussian XAverageConditional([SkipIfUniform] Bernoulli isPositive, double x)
		{
			return XAverageConditional(isPositive, Gaussian.PointMass(x));
		}

		/// <summary>
		/// EP message to 'x'
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		public static Gaussian XAverageConditional(bool isPositive, [SkipIfUniform] Gaussian x)
		{
			return XAverageConditional(Bernoulli.PointMass(isPositive), x);
		}

		[Skip]
		public static Gaussian XAverageConditionalInit()
		{
			return Gaussian.Uniform();
		}

		//-- VMP ---------------------------------------------------------------------------------------------
		const string NotSupportedMessage = "Variational Message Passing does not support the IsPositive factor with Gaussian distributions, since the factor is not conjugate to the Gaussian.";
		const string NotSupportedMessage2 = "Variational Message Passing does not support the IsPositive factor with stochastic output and Gaussian distributions, since the factor is not conjugate to the Gaussian.";

#if true
		//-----------------------------------------------------------------------------------------------------------
		// We can calculate VMP messages by implicitly modelling X as a TruncatedGaussian. This factor hides this
		// details of the implementation from the user. Note we cannot have multiple uses of this variable
		// in this case. To do this requires using truncated Gaussians explicitly. 

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_X">Previous outgoing message to 'X'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static Gaussian XAverageLogarithm(bool isPositive, [SkipIfUniform] Gaussian x, Gaussian to_X)
		{
			var prior = x / to_X;
			return XAverageConditional(isPositive, prior);
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_X">Previous outgoing message to 'X'.</param>
		/// <returns>The outgoing VMP message to the 'x' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static Gaussian XAverageLogarithm([SkipIfUniform] Bernoulli isPositive, [SkipIfUniform] Gaussian x, Gaussian to_X)
		{
			if (isPositive.IsPointMass) return XAverageLogarithm(isPositive.Point, x, to_X);
			if (isPositive.IsUniform()) return Gaussian.Uniform();
			throw new NotSupportedException(NotSupportedMessage2);
			var prior = x / to_X;
			return XAverageConditional(isPositive, prior);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isPositive">Message from isPositive</param>
		/// <param name="X">Message from X</param>
		/// <param name="to_X">Message to X</param>
		/// <remarks>The formula for the result is <c>sum_(sample,precision) p(sample,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.</remarks>
		/// <returns>Evidence</returns>
		[Quality(QualityBand.Preview)]
		public static double AverageLogFactor(Bernoulli isPositive, [SkipIfUniform] Gaussian X, Gaussian to_X)
		{
			if (isPositive.IsPointMass)
				return AverageLogFactor(isPositive.Point, X, to_X);
			var prior = X / to_X;
			//var tg = new TruncatedGaussian(prior);
			//tg.LowerBound = 0;
			// Remove the incorrect Gaussian entropy contribution and add the correct
			// truncated Gaussian entropy. Log(1)=0 so the factor itself does not contribute. 
			return X.GetAverageLog(X) - X.GetAverageLog(prior) + LogAverageFactor(isPositive, prior);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_X">Previous outgoing message to 'X'.</param>
		/// <returns>Evidence</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static double AverageLogFactor(bool isPositive, [SkipIfUniform] Gaussian x, Gaussian to_X)
		{
			if (isPositive)
				return AverageLogFactor_helper(x, to_X);
			else {
				x.MeanTimesPrecision *= -1.0;
				to_X.MeanTimesPrecision *= -1.0;
				return AverageLogFactor_helper(x, to_X);
			}
		}
#else
        
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[Skip]
		public static double AverageLogFactor(bool isPositive, [SkipIfUniform] TruncatedGaussian x)
		{
			return 0.0;
		}

		[NotSupported(NotSupportedMessage)]
		public static Gaussian XAverageLogarithm(bool isPositive, [SkipIfUniform] Gaussian x)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		[NotSupported(NotSupportedMessage)]
		public static double AverageLogFactor(bool isPositive, [SkipIfUniform] Gaussian x)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="isPositive">Incoming message from 'isPositive'.</param>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' with 'isPositive' integrated out.
		/// The formula is <c>sum_isPositive p(isPositive) factor(isPositive,x)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage2)]
		public static Gaussian XAverageLogarithm(Bernoulli isPositive, [SkipIfUniform] Gaussian x)
		{
			throw new NotSupportedException(NotSupportedMessage2);
		}

#endif

		/// <summary>
		/// VMP message to 'x'
		/// </summary>
		/// <param name="isPositive">Constant value for 'isPositive'.</param>
		/// <returns></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'x' conditioned on the given values.
		/// </para></remarks>
		[Quality(QualityBand.Preview)]
		public static TruncatedGaussian XAverageLogarithm(bool isPositive)
		{
			return XAverageConditional(isPositive);
		}
		[Quality(QualityBand.Preview)]
		public static TruncatedGaussian XAverageConditional(bool isPositive)
		{
			if (isPositive) {
				return new TruncatedGaussian(0, double.PositiveInfinity, 0, double.PositiveInfinity);
			} else {
				return new TruncatedGaussian(0, double.PositiveInfinity, double.NegativeInfinity, 0);
			}
		}
		[Quality(QualityBand.Preview)]
		public static TruncatedGaussian XAverageConditional([SkipIfUniform] Bernoulli isPositive)
		{
			if (isPositive.IsUniform()) return TruncatedGaussian.Uniform();
			if (isPositive.IsPointMass) return XAverageConditional(isPositive.Point);
			throw new NotSupportedException("Cannot return a TruncatedGaussian when isPositive is random");
		}

		/// <summary>
		/// VMP message to 'isPositive'
		/// </summary>
		/// <param name="x">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'isPositive' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isPositive' as the random arguments are varied.
		/// The formula is <c>proj[sum_(x) p(x) factor(isPositive,x)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="x"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static Bernoulli IsPositiveAverageLogarithm([SkipIfUniform] Gaussian x)
		{
			// same as BP if you use John Winn's rule.
			return IsPositiveAverageConditional(x);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="X">Incoming message from 'x'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_X">Previous outgoing message to 'X'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="X"/> is not a proper distribution</exception>
		[Quality(QualityBand.Preview)]
		public static double AverageLogFactor_helper([SkipIfUniform] Gaussian X, Gaussian to_X)
		{
			//if (!isPositive) throw new ArgumentException("VariationalMessagePassing requires isPositive=true", "isPositive");
			var prior = X / to_X;
			if (!prior.IsProper()) throw new ImproperMessageException(prior);
			var tg = new TruncatedGaussian(prior);
			tg.LowerBound = 0;
			// Remove the incorrect Gaussian entropy contribution and add the correct
			// truncated Gaussian entropy. Log(1)=0 so the factor itself does not contribute. 
			return X.GetAverageLog(X) - tg.GetAverageLog(tg);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="X">Incoming message from 'x'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		[Quality(QualityBand.Preview)]
		public static double AverageLogFactor(TruncatedGaussian X)
		{
			return 0;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.IsPositive"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "IsPositive", typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class IsPositiveOp_Proper
	{
		public static Gaussian XAverageConditional([SkipIfUniform] Bernoulli isPositive, [SkipIfUniform, Proper] Gaussian x)
		{
			return IsPositiveOp.XAverageConditional_Helper(isPositive, x, forceProper: true);
		}
	}
}
