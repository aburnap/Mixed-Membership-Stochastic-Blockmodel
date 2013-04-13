// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Not"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Not")]
	[Quality(QualityBand.Stable)]
	public static class BooleanNotOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(not,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool not, bool b)
		{
			return (not == Factor.Not(b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(not,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool not, bool b) { return LogAverageFactor(not, b); }

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(not,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool not, bool b) { return LogAverageFactor(not, b); }

		/// <summary>
		/// EP message to 'not'.
		/// </summary>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'not' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'not'.
		/// The formula is <c>int f(not,x) q(x) dx</c> where <c>x = (b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Bernoulli NotAverageConditional([SkipIfUniform] Bernoulli b)
		{
			return Bernoulli.FromLogOdds(-b.LogOdds);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="not">Incoming message from 'not'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (not)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="not"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli not)
		{
			return Bernoulli.FromLogOdds(-not.LogOdds);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (not)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool not)
		{
			return Bernoulli.PointMass(!not);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Incoming message from 'not'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="to_not">Outgoing message to 'not'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(not,b) p(not,b) factor(not,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli not, Bernoulli b, [Fresh] Bernoulli to_not)
		{
			return to_not.GetLogAverageOf(not);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Incoming message from 'not'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(not) p(not) factor(not,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli not, bool b)
		{
			return LogAverageFactor(b, not);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (not,b)</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool not, Bernoulli b)
		{
			return not ? b.GetLogProbFalse() : b.GetLogProbTrue();
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Incoming message from 'not'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(not) p(not) factor(not,b) / sum_not p(not) messageTo(not))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli not) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(not,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool not, Bernoulli b)
		{
			return LogAverageFactor(not, b);
		}


		//- VMP ---------------------------------------------------------------------------

		/// <summary>
		/// VMP message to 'not'.
		/// </summary>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'not' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'not'.
		/// The formula is <c>int log(f(not,x)) q(x) dx</c> where <c>x = (b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Bernoulli NotAverageLogarithm([SkipIfUniform] Bernoulli b)
		{
			return NotAverageConditional(b);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="not">Incoming message from 'not'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (not)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="not"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli not)
		{
			return BAverageConditional(not);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="not">Constant value for 'not'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (not)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageLogarithm(bool not)
		{
			return BAverageConditional(not);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (not,b)</c>.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }
	}
}
