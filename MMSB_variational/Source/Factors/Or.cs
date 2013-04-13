// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Or"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Or")]
	[Quality(QualityBand.Mature)]
	public static class BooleanOrOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool or, bool a, bool b)
		{
			return (or == Factor.Or(a, b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool or, bool a, bool b)
		{
			return LogAverageFactor(or, a, b);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool or, bool a, bool b)
		{
			return LogAverageFactor(or, a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Incoming message from 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(or) p(or) factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli or, bool a, bool b)
		{
			return or.GetLogProb(Factor.Or(a, b));
		}

		/// <summary>
		/// EP message to 'or'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[p(or) sum_(a,b) p(a,b) factor(or,a,b)]/p(or)</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageConditional(Bernoulli A, Bernoulli B)
		{
			return Bernoulli.FromLogOdds(Bernoulli.Or(A.LogOdds, B.LogOdds));
		}

		/// <summary>
		/// EP message to 'or'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[p(or) sum_(b) p(b) factor(or,a,b)]/p(or)</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageConditional(bool A, Bernoulli B)
		{
			if (A) return Bernoulli.PointMass(true);
			else return B;
		}

		/// <summary>
		/// EP message to 'or'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[p(or) sum_(a) p(a) factor(or,a,b)]/p(or)</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageConditional(Bernoulli A, bool B)
		{
			return OrAverageConditional(B, A);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(or,b) p(or,b) factor(or,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli or, Bernoulli B)
		{
			if (B.IsPointMass) return AAverageConditional(or, B.Point);
			return Bernoulli.FromLogOdds(Bernoulli.Gate(or.LogOdds, B.LogOdds));
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(b) p(b) factor(or,a,b)]/p(a)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageConditional(bool or, Bernoulli B)
		{
			if (B.IsPointMass) return AAverageConditional(or, B.Point);
			return AAverageConditional(Bernoulli.PointMass(or), B);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(or) p(or) factor(or,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli or, bool B)
		{
			if (or.IsPointMass) return AAverageConditional(or.Point, B);
			return Bernoulli.FromLogOdds(B ? 0 : or.LogOdds);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli AAverageConditional(bool or, bool B)
		{
			if (!B) return Bernoulli.PointMass(or);
			else if (or) return Bernoulli.Uniform();
			else throw new AllZeroException();
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(or,a) p(or,a) factor(or,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli or, Bernoulli A)
		{
			return AAverageConditional(or, A);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(or) p(or) factor(or,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli or, bool A)
		{
			return AAverageConditional(or, A);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(a) p(a) factor(or,a,b)]/p(b)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool or, Bernoulli A)
		{
			return AAverageConditional(or, A);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool or, bool A)
		{
			return AAverageConditional(or, A);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Incoming message from 'or'.</param>
		/// <param name="to_or">Outgoing message to 'or'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(or) p(or) factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli or, [Fresh] Bernoulli to_or)
		{
			return to_or.GetLogAverageOf(or);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool or, Bernoulli a, Bernoulli b)
		{
			Bernoulli to_or = OrAverageConditional(a, b);
			return to_or.GetLogProb(or);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool or, Bernoulli a, bool b)
		{
			Bernoulli to_or = OrAverageConditional(a, b);
			return to_or.GetLogProb(or);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(or,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool or, bool a, Bernoulli b)
		{
			return LogAverageFactor(or, b, a);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Incoming message from 'or'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(or) p(or) factor(or,a,b) / sum_or p(or) messageTo(or))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli or) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool or, Bernoulli a, Bernoulli b)
		{
			return LogAverageFactor(or, a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool or, Bernoulli a, bool b)
		{
			return LogAverageFactor(or, a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool or, bool a, Bernoulli b)
		{
			return LogEvidenceRatio(or, b, a);
		}

		//-- VMP ---------------------------------------------------------------------------------------

		const string NotSupportedMessage = "Variational Message Passing does not support an Or factor with fixed output.";

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(or,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'or'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(or,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageLogarithm(Bernoulli A, Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return OrAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'or'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(or,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageLogarithm(bool A, Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return OrAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'or'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'or' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'or' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(or,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli OrAverageLogarithm(Bernoulli A, bool B)
		{
			// same as BP if you use John Winn's rule.
			return OrAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'or' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_or p(or) factor(or,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli or, Bernoulli B)
		{
			// when 'or' is marginalized, the factor is proportional to exp((A|B)*or.LogOdds)
			return Bernoulli.FromLogOdds(or.LogOdds * B.GetProbFalse());
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'or' integrated out.
		/// The formula is <c>sum_or p(or) factor(or,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli or, bool B)
		{
			return Bernoulli.FromLogOdds(B ? 0.0 : or.LogOdds);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(or,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(BooleanOrOp.NotSupportedMessage)]
		public static Bernoulli AAverageLogarithm(bool or, Bernoulli B)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli AAverageLogarithm(bool or, bool B)
		{
			return AAverageConditional(or, B);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'or' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_or p(or) factor(or,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli or, Bernoulli A)
		{
			return AAverageLogarithm(or, A);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="or">Incoming message from 'or'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'or' integrated out.
		/// The formula is <c>sum_or p(or) factor(or,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="or"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli or, bool A)
		{
			return AAverageLogarithm(or, A);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(or,a,b)))</c>.
		/// </para></remarks>
		public static Bernoulli BAverageLogarithm(bool or, Bernoulli A)
		{
			return AAverageLogarithm(or, A);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="or">Constant value for 'or'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli BAverageLogarithm(bool or, bool A)
		{
			return AAverageLogarithm(or, A);
		}
	}
}
