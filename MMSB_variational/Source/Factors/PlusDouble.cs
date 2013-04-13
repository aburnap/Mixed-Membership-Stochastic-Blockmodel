// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Plus(double,double)"/> and <see cref="Factor.Difference(double,double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Plus", typeof(double), typeof(double))]
	[FactorMethod(new string[] { "A", "Sum", "B" }, typeof(Factor), "Difference", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class DoublePlusOp
	{
		// ----------------------------------------------------------------------------------------------------------------------
		// WrappedGaussian
		// ----------------------------------------------------------------------------------------------------------------------

		public static WrappedGaussian SumAverageConditional([SkipIfUniform] WrappedGaussian a, [SkipIfUniform] WrappedGaussian b)
		{
			if (a.Period != b.Period) throw new ArgumentException("a.Period ("+a.Period+") != b.Period ("+b.Period+")");
			WrappedGaussian result = WrappedGaussian.Uniform(a.Period);
			result.Gaussian = SumAverageConditional(a.Gaussian, b.Gaussian);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian SumAverageConditional([SkipIfUniform] WrappedGaussian a, double b)
		{
			WrappedGaussian result = WrappedGaussian.Uniform(a.Period);
			result.Gaussian = SumAverageConditional(a.Gaussian, b);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian SumAverageConditional(double a, [SkipIfUniform] WrappedGaussian b)
		{
			return SumAverageConditional(b, a);
		}
		public static WrappedGaussian AAverageConditional([SkipIfUniform] WrappedGaussian sum, [SkipIfUniform] WrappedGaussian b)
		{
			if (sum.Period != b.Period) throw new ArgumentException("sum.Period ("+sum.Period+") != b.Period ("+b.Period+")");
			WrappedGaussian result = WrappedGaussian.Uniform(sum.Period);
			result.Gaussian = AAverageConditional(sum.Gaussian, b.Gaussian);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian AAverageConditional([SkipIfUniform] WrappedGaussian sum, double b)
		{
			WrappedGaussian result = WrappedGaussian.Uniform(sum.Period);
			result.Gaussian = AAverageConditional(sum.Gaussian, b);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian BAverageConditional([SkipIfUniform] WrappedGaussian sum, [SkipIfUniform] WrappedGaussian a)
		{
			return AAverageConditional(sum, a);
		}
		public static WrappedGaussian BAverageConditional([SkipIfUniform] WrappedGaussian sum, double a)
		{
			return AAverageConditional(sum, a);
		}

		// ----------------------------------------------------------------------------------------------------------------------
		// TruncatedGaussian
		// ----------------------------------------------------------------------------------------------------------------------

		public static TruncatedGaussian SumAverageConditional(double a, [SkipIfUniform] TruncatedGaussian b)
		{
			Gaussian prior = b.Gaussian;
			Gaussian post = SumAverageConditional(a, prior);
			TruncatedGaussian result = b;
			result.Gaussian = post;
			result.LowerBound += a;
			result.UpperBound += a;
			return result;
		}
		public static TruncatedGaussian SumAverageConditional([SkipIfUniform] TruncatedGaussian a, double b)
		{
			return SumAverageConditional(b, a);
		}
		public static TruncatedGaussian AAverageConditional([SkipIfUniform] TruncatedGaussian sum, double b)
		{
			Gaussian prior = sum.Gaussian;
			Gaussian post = AAverageConditional(prior, b);
			TruncatedGaussian result = sum;
			result.Gaussian = post;
			result.LowerBound -= b;
			result.UpperBound -= b;
			return result;
		}
		public static TruncatedGaussian AAverageConditional(double sum, [SkipIfUniform] TruncatedGaussian b)
		{
			Gaussian prior = b.Gaussian;
			Gaussian post = AAverageConditional(sum, prior);
			return new TruncatedGaussian(post, sum - b.UpperBound, sum - b.LowerBound);
		}
		public static TruncatedGaussian BAverageConditional([SkipIfUniform] TruncatedGaussian sum, double a)
		{
			return AAverageConditional(sum, a);
		}
		public static TruncatedGaussian BAverageConditional(double sum, [SkipIfUniform] TruncatedGaussian a)
		{
			return AAverageConditional(sum, a);
		}
		public static double LogAverageFactor(double sum, [SkipIfUniform] TruncatedGaussian a, double b)
		{
			TruncatedGaussian to_sum = SumAverageConditional(a, b);
			return to_sum.GetLogProb(sum);
		}
		public static double LogAverageFactor(double sum, double a, [SkipIfUniform] TruncatedGaussian b)
		{
			TruncatedGaussian to_sum = SumAverageConditional(a, b);
			return to_sum.GetLogProb(sum);
		}
		public static double LogAverageFactor([SkipIfUniform] TruncatedGaussian sum, double a, double b)
		{
			return sum.GetLogProb(Factor.Plus(a, b));
		}

		// ----------------------------------------------------------------------------------------------------------------------
		// Gaussian 
		// ----------------------------------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum) p(Sum) factor(Sum,A,B))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Gaussian Sum, [Fresh] Gaussian to_sum)
		{
			return to_sum.GetLogAverageOf(Sum);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_a">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(Sum,A,B))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double Sum, [SkipIfUniform] Gaussian a, [Fresh] Gaussian to_a)
		{
			return to_a.GetLogAverageOf(a);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="to_b">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double Sum, double a, Gaussian b, [Fresh] Gaussian to_b)
		{
			return LogAverageFactor(Sum, b, to_b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double Sum, double a, double b)
		{
			return (Sum == Factor.Plus(a, b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(A,B) p(A,B) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian SumAverageConditional([SkipIfUniform] Gaussian a, [SkipIfUniform] Gaussian b)
		{
			if (a.IsPointMass) return SumAverageConditional(a.Point, b);
			if (b.IsPointMass) return SumAverageConditional(a, b.Point);
			if (a.Precision == 0) return Gaussian.FromNatural(a.MeanTimesPrecision, 0);
			double prec = a.Precision + b.Precision;
			if (prec <= 0) throw new ImproperDistributionException(a.IsProper() ? b : a);
			return Gaussian.FromNatural((a.MeanTimesPrecision*b.Precision + b.MeanTimesPrecision*a.Precision)/prec, a.Precision*b.Precision/prec);
		}

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(B) p(B) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian SumAverageConditional(double a, [SkipIfUniform] Gaussian b)
		{
			if (b.IsPointMass) return SumAverageConditional(a, b.Point);
			return Gaussian.FromNatural(b.MeanTimesPrecision + a*b.Precision, b.Precision);
		}

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing EP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(A) p(A) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian SumAverageConditional([SkipIfUniform] Gaussian a, double b)
		{
			return SumAverageConditional(b, a);
		}
		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing EP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Sum' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian SumAverageConditional(double a, double b)
		{
			return Gaussian.PointMass(a + b);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(Sum,B) p(Sum,B) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Sum, [SkipIfUniform] Gaussian b)
		{
			if (Sum.IsPointMass) return AAverageConditional(Sum.Point, b);
			if (b.IsPointMass) return AAverageConditional(Sum, b.Point);
			if (Sum.Precision == 0) return Gaussian.FromNatural(Sum.MeanTimesPrecision, 0);
			double prec = Sum.Precision + b.Precision;
			if (prec <= 0) throw new ImproperDistributionException(Sum.IsProper() ? b : Sum);
			return Gaussian.FromNatural((Sum.MeanTimesPrecision*b.Precision - b.MeanTimesPrecision*Sum.Precision)/prec, Sum.Precision*b.Precision/prec);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(B) p(B) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional(double Sum, [SkipIfUniform] Gaussian b)
		{
			if (b.IsPointMass) return AAverageConditional(Sum, b.Point);
			return Gaussian.FromNatural(Sum*b.Precision - b.MeanTimesPrecision, b.Precision);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing EP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(Sum) p(Sum) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Sum, double b)
		{
			if (Sum.IsPointMass) return AAverageConditional(Sum, b);
			return Gaussian.FromNatural(Sum.MeanTimesPrecision - b*Sum.Precision, Sum.Precision);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing EP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageConditional(double Sum, double b)
		{
			return Gaussian.PointMass(Sum - b);
		}

		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(Sum,A) p(Sum,A) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Sum, [SkipIfUniform] Gaussian a)
		{
			return AAverageConditional(Sum, a);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(A) p(A) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional(double Sum, [SkipIfUniform] Gaussian a)
		{
			return AAverageConditional(Sum, a);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <returns>The outgoing EP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(Sum) p(Sum) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Sum, double a)
		{
			return AAverageConditional(Sum, a);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <returns>The outgoing EP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian BAverageConditional(double Sum, double a)
		{
			return AAverageConditional(Sum, a);
		}
	}

	/// <summary>
	/// Evidence methods for Factor.Plus only
	/// </summary>
	[FactorMethod(typeof(Factor), "Plus", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class DoublePlusEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum) p(Sum) factor(Sum,A,B) / sum_Sum p(Sum) messageTo(Sum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian Sum) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double Sum, double a, double b) { return DoublePlusOp.LogAverageFactor(Sum, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="to_b">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double Sum, double a, Gaussian b, [Fresh] Gaussian to_b) { return DoublePlusOp.LogAverageFactor(b, to_b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="to_a">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double Sum, Gaussian a, [Fresh] Gaussian to_a) { return DoublePlusOp.LogAverageFactor(a, to_a); }
	}

	/// <summary>
	/// Evidence methods for Factor.Difference only
	/// </summary>
	[FactorMethod(typeof(Factor), "Difference", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class DoubleMinusEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="difference">Incoming message from 'difference'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(difference) p(difference) factor(difference,a,b) / sum_difference p(difference) messageTo(difference))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian difference) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="difference">Constant value for 'difference'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(difference,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double difference, double a, double b) { return DoublePlusOp.LogAverageFactor(a, difference, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="difference">Constant value for 'difference'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="to_a">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(difference,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double difference, Gaussian a, [Fresh] Gaussian to_a) { return DoublePlusOp.LogAverageFactor(a, to_a); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="difference">Constant value for 'difference'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="to_b">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(difference,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double difference, double a, Gaussian b, [Fresh] Gaussian to_b) { return DoublePlusOp.LogAverageFactor(b, to_b); }
	}

	/// <summary>
	/// These operator methods are only valid for plus, not minus.
	/// </summary>
	[FactorMethod(typeof(Factor), "Plus", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class DoublePlusVmpOp
	{
		// ----------------------------------------------------------------------------------------------------------------------
		// WrappedGaussian
		// ----------------------------------------------------------------------------------------------------------------------

		[Skip]
		public static double AverageLogFactor(WrappedGaussian sum) { return 0.0; }
		public static WrappedGaussian SumAverageLogarithm([SkipIfUniform] WrappedGaussian a, [SkipIfUniform] WrappedGaussian b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		public static WrappedGaussian SumAverageLogarithm([SkipIfUniform] WrappedGaussian a, double b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		public static WrappedGaussian SumAverageLogarithm(double a, [SkipIfUniform] WrappedGaussian b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		public static WrappedGaussian AAverageLogarithm([SkipIfUniform] WrappedGaussian sum, [SkipIfUniform] WrappedGaussian b)
		{
			if (sum.Period != b.Period) throw new ArgumentException("sum.Period ("+sum.Period+") != b.Period ("+b.Period+")");
			WrappedGaussian result = WrappedGaussian.Uniform(sum.Period);
			result.Gaussian = AAverageLogarithm(sum.Gaussian, b.Gaussian);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian AAverageLogarithm([SkipIfUniform] WrappedGaussian sum, double b)
		{
			WrappedGaussian result = WrappedGaussian.Uniform(sum.Period);
			result.Gaussian = AAverageLogarithm(sum.Gaussian, b);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian BAverageLogarithm([SkipIfUniform] WrappedGaussian sum, [SkipIfUniform] WrappedGaussian a)
		{
			return AAverageLogarithm(sum, a);
		}
		public static WrappedGaussian BAverageLogarithm([SkipIfUniform] WrappedGaussian sum, double a)
		{
			return AAverageLogarithm(sum, a);
		}

		// ----------------------------------------------------------------------------------------------------------------------
		// Gaussian
		// ----------------------------------------------------------------------------------------------------------------------

		internal const string NotSupportedMessage = "Variational Message Passing does not support a Plus factor with fixed output and two random inputs.";

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Gaussian sum) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(double sum, double a, double b) { return DoublePlusOp.LogAverageFactor(sum, a, b); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(double sum, Gaussian a) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(double sum, double a, Gaussian b) { return 0.0; }

		/// <summary>
		/// VMP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[sum_(A,B) p(A,B) factor(Sum,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian SumAverageLogarithm([SkipIfUniform] Gaussian a, [SkipIfUniform] Gaussian b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'Sum'
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[sum_(B) p(B) factor(Sum,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian SumAverageLogarithm(double a, [SkipIfUniform] Gaussian b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing VMP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[sum_(A) p(A) factor(Sum,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian SumAverageLogarithm([SkipIfUniform] Gaussian a, double b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'Sum'
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing VMP message to the 'Sum' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'Sum' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian SumAverageLogarithm(double a, double b)
		{
			return DoublePlusOp.SumAverageConditional(a, b);
		}

		[Skip]
		public static Gaussian SumDeriv(Gaussian a, double b) { return Gaussian.Uniform(); }
		[Skip]
		public static Gaussian SumDeriv(double a, Gaussian b) { return Gaussian.Uniform(); }
		public static Gaussian SumDeriv([SkipIfUniform, Proper] Gaussian Sum, Gaussian a, Gaussian a_deriv, Gaussian to_a, [Proper] Gaussian b, Gaussian b_deriv, Gaussian to_b)
		{
			double sa = a_deriv.MeanTimesPrecision;
			double ta = a.MeanTimesPrecision - sa*to_a.MeanTimesPrecision;
			double sb = b_deriv.MeanTimesPrecision;
			double tb = b.MeanTimesPrecision - sb*to_b.MeanTimesPrecision;
			double va = 1/a.Precision;
			double vb = 1/b.Precision;
			double aa = va*sa*Sum.Precision;
			double ab = vb*sb*Sum.Precision;
			double ba = va*sa;
			double bb = vb*sb;
			double deriv = (ba*(1-ab)+bb*(1-aa))/(1-aa*ab);
			return Gaussian.FromNatural(deriv-1, 0);
		}
		public static Gaussian AAverageLogarithm([SkipIfUniform, Proper] Gaussian Sum, Gaussian a, Gaussian a_deriv, Gaussian to_a, [Proper] Gaussian b, Gaussian b_deriv, Gaussian to_b)
		{
			double sa = a_deriv.MeanTimesPrecision+1;
			double ta = a.MeanTimesPrecision - sa*to_a.MeanTimesPrecision;
			double sb = b_deriv.MeanTimesPrecision+1;
			double tb = b.MeanTimesPrecision - sb*to_b.MeanTimesPrecision;
			double va = 1/a.Precision;
			double vb = 1/b.Precision;
			double aa = va*sa*Sum.Precision;
			double ab = vb*sb*Sum.Precision;
			double ba = va*(ta + sa*Sum.MeanTimesPrecision);
			double bb = vb*(tb + sb*Sum.MeanTimesPrecision);
			double ma = (ba - bb*aa)/(1 - aa*ab);
			double mu = (ma/va - ta)/sa;
			return new Gaussian(mu/Sum.Precision, 1/Sum.Precision);
		}

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'A'.
		/// Because the factor is deterministic, 'Sum' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(B) p(B) log(sum_Sum p(Sum) factor(Sum,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform, Proper] Gaussian Sum, [Proper] Gaussian b)
		{
			Gaussian result = new Gaussian();
			if (Sum.IsUniform()) {
				result.SetToUniform();
			} else if (Sum.Precision < 0) {
				throw new ImproperMessageException(Sum);
			} else {
				// p(a|sum,b) = N(E[sum] - E[b], var(sum) )
				double ms, vs;
				double mb = b.GetMean();
				Sum.GetMeanAndVariance(out ms, out vs);
				result.SetMeanAndVariance(ms - mb, vs);
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'A'.
		/// The formula is <c>exp(sum_(B) p(B) log(factor(Sum,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gaussian AAverageLogarithm(double Sum, [Proper] Gaussian b)
		{
			// Throw an exception rather than return a meaningless point mass.
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' with 'Sum' integrated out.
		/// The formula is <c>sum_Sum p(Sum) factor(Sum,A,B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform, Proper] Gaussian Sum, double b)
		{
			return AAverageLogarithm(Sum, Gaussian.PointMass(b));
		}
		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageLogarithm(double Sum, double b)
		{
			return DoublePlusOp.AAverageConditional(Sum, b);
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'B'.
		/// Because the factor is deterministic, 'Sum' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(A) p(A) log(sum_Sum p(Sum) factor(Sum,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform, Proper] Gaussian Sum, [Proper] Gaussian a)
		{
			return AAverageLogarithm(Sum, a);
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'B'.
		/// The formula is <c>exp(sum_(A) p(A) log(factor(Sum,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gaussian BAverageLogarithm(double Sum, [Proper] Gaussian a)
		{
			return AAverageLogarithm(Sum, a);
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="Sum">Incoming message from 'Sum'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' with 'Sum' integrated out.
		/// The formula is <c>sum_Sum p(Sum) factor(Sum,A,B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sum"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform, Proper] Gaussian Sum, double a)
		{
			return AAverageLogarithm(Sum, a);
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="Sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian BAverageLogarithm(double Sum, double a)
		{
			return AAverageLogarithm(Sum, a);
		}

	}
	/// <summary>
	/// Provides outgoing VMP messages for <see cref="Factor.Difference(double, double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Difference", typeof(double), typeof(double))]
	[Quality(QualityBand.Stable)]
	public static class DoubleMinusVmpOp
	{
		//-- VMP ----------------------------------------------------------------------------------------------

		internal const string NotSupportedMessage = "Variational Message Passing does not support a Minus factor with fixed output and two random inputs.";

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(difference,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'difference'
		/// </summary>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'difference' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'difference' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(difference,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian DifferenceAverageLogarithm([SkipIfUniform] Gaussian a, [SkipIfUniform] Gaussian b)
		{
			return DoublePlusOp.AAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'difference'
		/// </summary>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'difference' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'difference' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(difference,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian DifferenceAverageLogarithm(double a, [SkipIfUniform] Gaussian b)
		{
			return DoublePlusOp.AAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'difference'
		/// </summary>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'difference' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'difference' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(difference,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian DifferenceAverageLogarithm([SkipIfUniform] Gaussian a, double b)
		{
			return DoublePlusOp.AAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'difference'
		/// </summary>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'difference' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'difference' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian DifferenceAverageLogarithm(double a, double b)
		{
			return DoublePlusOp.AAverageConditional(a, b);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Difference">Incoming message from 'difference'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'difference' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_difference p(difference) factor(difference,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Difference"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform, Proper] Gaussian Difference, [Proper] Gaussian b)
		{
			Gaussian result = new Gaussian();
			if (Difference.IsUniform()) {
				result.SetToUniform();
			} else if (Difference.Precision < 0) {
				throw new ImproperMessageException(Difference);
			} else {
				// p(a|diff,b) = N(E[diff] + E[b], var(diff) )
				double ms, vs;
				double mb = b.GetMean();
				Difference.GetMeanAndVariance(out ms, out vs);
				result.SetMeanAndVariance(ms + mb, vs);
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Difference">Incoming message from 'difference'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'difference' integrated out.
		/// The formula is <c>sum_difference p(difference) factor(difference,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Difference"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform, Proper] Gaussian Difference, double b)
		{
			return AAverageLogarithm(Difference, Gaussian.PointMass(b));
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Difference">Constant value for 'difference'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(difference,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gaussian AAverageLogarithm(double Difference, [Proper] Gaussian b)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Difference">Constant value for 'difference'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageLogarithm(double Difference, double b)
		{
			return Gaussian.PointMass(Difference + b);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Difference">Incoming message from 'difference'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'difference' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_difference p(difference) factor(difference,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Difference"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform, Proper] Gaussian Difference, [Proper] Gaussian a)
		{
			Gaussian result = new Gaussian();
			if (Difference.IsUniform()) {
				result.SetToUniform();
			} else if (Difference.Precision < 0) {
				throw new ImproperMessageException(Difference);
			} else {
				// p(b|diff,a) = N(E[a] - E[diff], var(diff) )
				double ms, vs;
				double ma = a.GetMean();
				Difference.GetMeanAndVariance(out ms, out vs);
				result.SetMeanAndVariance(ma - ms, vs);
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Difference">Incoming message from 'difference'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'difference' integrated out.
		/// The formula is <c>sum_difference p(difference) factor(difference,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Difference"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform, Proper] Gaussian Difference, double a)
		{
			return BAverageLogarithm(Difference, Gaussian.PointMass(a));
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Difference">Constant value for 'difference'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(difference,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		[NotSupported(NotSupportedMessage)]
		public static Gaussian BAverageLogarithm(double Difference, [Proper] Gaussian a)
		{
			// Throw an exception rather than return a meaningless point mass.
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Difference">Constant value for 'difference'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian BAverageLogarithm(double Difference, double a)
		{
			return Gaussian.PointMass(a - Difference);
		}
	}
}
