// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.And"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "And")]
	[Quality(QualityBand.Mature)]
	public static class BooleanAndOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_() p() factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool and, bool a, bool b)
		{
			return (and == Factor.And(a,b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(and,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool and, bool a, bool b)
		{
			return LogAverageFactor(and, a, b);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		public static double AverageLogFactor(bool and, bool a, bool b)
		{
			return LogAverageFactor(and, a, b);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Incoming message from 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(and) p(and) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli and, bool a, bool b)
		{
			return and.GetLogProb(Factor.And(a, b));
		}

		/// <summary>
		/// EP message to 'and'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'and'.
		/// The formula is <c>int f(and,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageConditional(Bernoulli A, Bernoulli B)
		{
			return Bernoulli.FromLogOdds(-Bernoulli.Or(-A.LogOdds, -B.LogOdds));
		}

		/// <summary>
		/// EP message to 'and'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'and'.
		/// The formula is <c>int f(and,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageConditional(bool A, Bernoulli B)
		{
			if (A) return B;
			else return Bernoulli.PointMass(false);
		}
		/// <summary>
		/// EP message to 'and'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'and'.
		/// The formula is <c>int f(and,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageConditional(Bernoulli A, bool B)
		{
			return AndAverageConditional(B, A);
		}

		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli and, Bernoulli B)
		{
			if (B.IsPointMass) return AAverageConditional(and, B.Point);
			return Bernoulli.FromLogOdds(-Bernoulli.Gate(-and.LogOdds, -B.LogOdds));
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageConditional(bool and, Bernoulli B)
		{
			if (B.IsPointMass) return AAverageConditional(and, B.Point);
			return AAverageConditional(Bernoulli.PointMass(and), B);
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli and, bool B)
		{
			if (and.IsPointMass) return AAverageConditional(and.Point, B);
			return Bernoulli.FromLogOdds(B ? and.LogOdds : 0);
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageConditional(bool and, bool B)
		{
			if (B) return Bernoulli.PointMass(and);
			else if (!and) return Bernoulli.Uniform();
			else throw new AllZeroException();
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli and, Bernoulli A)
		{
			return AAverageConditional(and, A);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool and, Bernoulli A)
		{
			return AAverageConditional(and, A);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli and, bool A)
		{
			return AAverageConditional(and, A);
		}

		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool and, bool A)
		{
			return AAverageConditional(and, A);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Incoming message from 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(and,a,b) p(and,a,b) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli and, Bernoulli a, Bernoulli b)
		{
			Bernoulli to_and = AndAverageConditional(a, b);
			return to_and.GetLogAverageOf(and);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Incoming message from 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(and,b) p(and,b) factor(and,a,b))</c>.
		/// </para></remarks>
#pragma warning disable 1591
		public static double LogAverageFactor(Bernoulli and, bool a, Bernoulli b)
		{
			Bernoulli to_and = AndAverageConditional(a, b);
			return to_and.GetLogAverageOf(and);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Incoming message from 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(and,a) p(and,a) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli and, Bernoulli a, bool b)
		{
			Bernoulli to_and = AndAverageConditional(a, b);
			return to_and.GetLogAverageOf(and);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool and, Bernoulli a, Bernoulli b)
		{
			Bernoulli to_and = AndAverageConditional(a, b);
			return to_and.GetLogProb(and);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool and, Bernoulli a, bool b)
		{
			Bernoulli to_and = AndAverageConditional(a, b);
			return to_and.GetLogProb(and);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(and,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool and, bool a, Bernoulli b)
		{
			return LogAverageFactor(and, b, a);
		}
#pragma warning restore 1591

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Incoming message from 'and'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(and,a,b) p(and,a,b) factor(and,a,b) / sum_and p(and) messageTo(and))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli and)		{			return 0.0;		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(and,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool and, Bernoulli a, Bernoulli b)
		{
			return LogAverageFactor(and, a, b);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(and,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool and, Bernoulli a, bool b)
		{
			return LogAverageFactor(and, a, b);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(and,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool and, bool a, Bernoulli b)
		{
			return LogEvidenceRatio(and, b, a);
		}

		//-- VMP ---------------------------------------------------------------------------------------------

		const string NotSupportedMessage = "Variational Message Passing does not support an And factor with fixed output.";

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'and'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'and' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(and,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageLogarithm(Bernoulli A, Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return AndAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'and'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'and'.
		/// The formula is <c>int log(f(and,x)) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageLogarithm(bool A, Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return AndAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'and'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'and' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'and'.
		/// The formula is <c>int log(f(and,x)) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AndAverageLogarithm(Bernoulli A, bool B)
		{
			// same as BP if you use John Winn's rule.
			return AndAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli and, Bernoulli B)
		{
			// when 'and' is marginalized, the factor is proportional to exp(A*B*and.LogOdds)
			return Bernoulli.FromLogOdds(and.LogOdds * B.GetProbTrue());
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli and, bool B)
		{
			return Bernoulli.FromLogOdds(B ? and.LogOdds : 0.0);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		[NotSupported(BooleanAndOp.NotSupportedMessage)]
		public static Bernoulli AAverageLogarithm(bool and, Bernoulli B)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (and,b)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageLogarithm(bool and, bool B)
		{
			return AAverageConditional(and, B);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli and, Bernoulli A)
		{
			return AAverageLogarithm(and, A);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="and">Incoming message from 'and'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="and"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli and, bool A)
		{
			return AAverageLogarithm(and, A);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		[NotSupported(BooleanAndOp.NotSupportedMessage)]
		public static Bernoulli BAverageLogarithm(bool and, Bernoulli A)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="and">Constant value for 'and'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (and,a)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageLogarithm(bool and, bool A)
		{
			return AAverageLogarithm(and, A);
		}
	}
}
