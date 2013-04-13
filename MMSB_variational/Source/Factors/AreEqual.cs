// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	// this factor is symmetric among all three arguments.
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.AreEqual(bool,bool)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "AreEqual", typeof(bool), typeof(bool))]
	[Quality(QualityBand.Stable)]
	public static class BooleanAreEqualOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_() p() factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, bool a, bool b)
		{
			return (areEqual == Factor.AreEqual(a, b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, bool a, bool b)
		{
			return LogAverageFactor(areEqual, a, b);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		public static double AverageLogFactor(bool areEqual, bool a, bool b)
		{
			return LogAverageFactor(areEqual, a, b);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli areEqual, bool a, bool b)
		{
			return areEqual.GetLogProb(Factor.AreEqual(a, b));
		}

		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageConditional([SkipIfUniform] Bernoulli A, [SkipIfUniform] Bernoulli B)
		{
			return Bernoulli.FromLogOdds(Bernoulli.LogitProbEqual(A.LogOdds, B.LogOdds));
		}
		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageConditional(bool A, [SkipIfUniform] Bernoulli B)
		{
			return Bernoulli.FromLogOdds(A ? B.LogOdds : -B.LogOdds);
		}
		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageConditional([SkipIfUniform] Bernoulli A, bool B)
		{
			return AreEqualAverageConditional(B, A);
		}
		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(bool A, bool B)
		{
			return Bernoulli.PointMass(Factor.AreEqual(A, B));
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli areEqual, [SkipIfUniform] Bernoulli B)
		{
			return AreEqualAverageConditional(areEqual, B);
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional(bool areEqual, [SkipIfUniform] Bernoulli B)
		{
			return AreEqualAverageConditional(areEqual, B);
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Bernoulli AAverageConditional([SkipIfUniform] Bernoulli areEqual, bool B)
		{
			return AreEqualAverageConditional(areEqual, B);
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageConditional(bool areEqual, bool B)
		{
			return AreEqualAverageConditional(areEqual, B);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli areEqual, [SkipIfUniform] Bernoulli A)
		{
			return AAverageConditional(areEqual, A);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional(bool areEqual, [SkipIfUniform] Bernoulli A)
		{
			return AAverageConditional(areEqual, A);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Bernoulli BAverageConditional([SkipIfUniform] Bernoulli areEqual, bool A)
		{
			return AAverageConditional(areEqual, A);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageConditional(bool areEqual, bool A)
		{
			return AAverageConditional(areEqual, A);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <param name="to_areEqual">Outgoing message to 'areEqual'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli areEqual, [Fresh] Bernoulli to_areEqual)
		{
			return to_areEqual.GetLogAverageOf(areEqual);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Bernoulli A, Bernoulli B, [Fresh] Bernoulli to_A)
		{
			//Bernoulli to_A = AAverageConditional(areEqual, B);
			return A.GetLogAverageOf(to_A);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Bernoulli A, bool B, [Fresh] Bernoulli to_A)
		{
			//Bernoulli toA = AAverageConditional(areEqual, B);
			return A.GetLogAverageOf(to_A);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_B">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, bool A, Bernoulli B, [Fresh] Bernoulli to_B)
		{
			return LogAverageFactor(areEqual, B, A, to_B);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual,a,b) p(areEqual,a,b) factor(areEqual,a,b) / sum_areEqual p(areEqual) messageTo(areEqual))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli areEqual)		{			return 0.0;		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, Bernoulli A, Bernoulli B, [Fresh] Bernoulli to_A)
		{
			return LogAverageFactor(areEqual, A, B, to_A);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_B">Outgoing message to 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, bool A, Bernoulli B, [Fresh] Bernoulli to_B)
		{
			return LogAverageFactor(areEqual, A, B, to_B);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, Bernoulli A, bool B, [Fresh] Bernoulli to_A)
		{
			return LogEvidenceRatio(areEqual, B, A, to_A);
		}

		//- VMP -------------------------------------------------------------------------------------------------------------

		const string NotSupportedMessage = "Variational Message Passing does not support an AreEqual factor with fixed output.";

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
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageLogarithm([SkipIfUniform] Bernoulli A, [SkipIfUniform] Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return AreEqualAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int log(f(areEqual,x)) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageLogarithm(bool A, [SkipIfUniform] Bernoulli B)
		{
			// same as BP if you use John Winn's rule.
			return AreEqualAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int log(f(areEqual,x)) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Bernoulli AreEqualAverageLogarithm([SkipIfUniform] Bernoulli A, bool B)
		{
			// same as BP if you use John Winn's rule.
			return AreEqualAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int log(f(areEqual,x)) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(bool A, bool B)
		{
			// same as BP if you use John Winn's rule.
			return AreEqualAverageConditional(A, B);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli areEqual, [SkipIfUniform] Bernoulli B)
		{
			if (areEqual.IsPointMass) return AAverageLogarithm(areEqual.Point, B);
			// when AreEqual is marginalized, the factor is proportional to exp((A==B)*areEqual.LogOdds)
			return Bernoulli.FromLogOdds(areEqual.LogOdds * (2 * B.GetProbTrue() - 1));
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Bernoulli AAverageLogarithm([SkipIfUniform] Bernoulli areEqual, bool B)
		{
			return Bernoulli.FromLogOdds(B ? areEqual.LogOdds : -areEqual.LogOdds);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(BooleanAreEqualOp.NotSupportedMessage)]
		public static Bernoulli AAverageLogarithm(bool areEqual, [SkipIfUniform] Bernoulli B)
		{
			if (B.IsPointMass) return AAverageLogarithm(areEqual, B.Point);
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int log(f(a,x)) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		public static Bernoulli AAverageLogarithm(bool areEqual, bool B)
		{
			return AAverageConditional(areEqual, B);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli areEqual, [SkipIfUniform] Bernoulli A)
		{
			return AAverageLogarithm(areEqual, A);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Bernoulli BAverageLogarithm([SkipIfUniform] Bernoulli areEqual, bool A)
		{
			return AAverageLogarithm(areEqual, A);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[NotSupported(BooleanAreEqualOp.NotSupportedMessage)]
		public static Bernoulli BAverageLogarithm(bool areEqual, [SkipIfUniform] Bernoulli A)
		{
			return AAverageLogarithm(areEqual, A);
		}
		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int log(f(b,x)) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		public static Bernoulli BAverageLogarithm(bool areEqual, bool A)
		{
			return AAverageLogarithm(areEqual, A);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.AreEqual(int, int)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "AreEqual", typeof(int), typeof(int))]
#if !SILVERLIGHT
	[FactorMethod(typeof(EnumSupport), "AreEqual<>")]
#endif
	[Quality(QualityBand.Preview)]
	public static class DiscreteAreEqualOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_() p() factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, int a, int b)
		{
			return (areEqual == Factor.AreEqual(a, b)) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, int a, int b)
		{
			return LogAverageFactor(areEqual, a, b);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		public static double AverageLogFactor(bool areEqual, int a, int b)
		{
			return LogAverageFactor(areEqual, a, b);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli areEqual, int a, int b)
		{
			return areEqual.GetLogProb(Factor.AreEqual(a, b));
		}

		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(Discrete A, Discrete B)
		{
			return Bernoulli.FromLogOdds(MMath.Logit(A.ProbEqual(B)));
		}
		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(int A, Discrete B)
		{
			return Bernoulli.FromLogOdds(MMath.Logit(B[A]));
		}
		/// <summary>
		/// EP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'areEqual'.
		/// The formula is <c>int f(areEqual,x) q(x) dx</c> where <c>x = (a,b)</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageConditional(Discrete A, int B)
		{
			return AreEqualAverageConditional(B, A);
		}

		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete AAverageConditional([SkipIfUniform] Bernoulli areEqual, Discrete B, Discrete result)
		{
			if (areEqual.IsPointMass) return AAverageConditional(areEqual.Point, B, result);
			if (result == default(Discrete)) result = Distributions.Discrete.Uniform(B.Dimension, B.Sparsity);
			double p = areEqual.GetProbTrue();
			Vector probs = result.GetWorkspace();
			probs = B.GetProbs(probs);
			probs.SetToProduct(probs, 2.0 * p - 1.0);
			probs.SetToSum(probs, 1.0 - p);
			result.SetProbs(probs);
			return result;
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete AAverageConditional([SkipIfUniform] Bernoulli areEqual, int B, Discrete result)
		{
			if (areEqual.IsPointMass) return AAverageConditional(areEqual.Point, B, result);
			Vector probs = result.GetWorkspace();
			double p = areEqual.GetProbTrue();
			probs.SetAllElementsTo(1 - p);
			probs[B] = p;
			result.SetProbs(probs);
			return result;
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		public static Discrete AAverageConditional(bool areEqual, Discrete B, Discrete result)
		{
			if (B.IsPointMass) return AAverageConditional(areEqual, B.Point, result);
			if (result == default(Discrete)) result = Distributions.Discrete.Uniform(B.Dimension, B.Sparsity);
			if (areEqual) result.SetTo(B);
			else {
				Vector probs = result.GetWorkspace();
				probs = B.GetProbs(probs);
				probs.SetToDifference(1.0, probs);
				result.SetProbs(probs);
			}
			return result;
		}
		/// <summary>
		/// EP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'a'.
		/// The formula is <c>int f(a,x) q(x) dx</c> where <c>x = (areEqual,b)</c>.
		/// </para></remarks>
		public static Discrete AAverageConditional(bool areEqual, int B, Discrete result)
		{
			if (areEqual) result.Point = B;
			else if (result.Dimension == 2) result.Point = 1 - B;
			else {
				Vector probs = result.GetWorkspace();
				probs.SetAllElementsTo(1);
				probs[B] = 0;
				result.SetProbs(probs);
			}
			return result;
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete BAverageConditional([SkipIfUniform] Bernoulli areEqual, Discrete A, Discrete result)
		{
			return AAverageConditional(areEqual, A, result);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete BAverageConditional([SkipIfUniform] Bernoulli areEqual, int A, Discrete result)
		{
			return AAverageConditional(areEqual, A, result);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(bool areEqual, Discrete A, Discrete result)
		{
			return AAverageConditional(areEqual, A, result);
		}
		/// <summary>
		/// EP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'b'.
		/// The formula is <c>int f(b,x) q(x) dx</c> where <c>x = (areEqual,a)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(bool areEqual, int A, Discrete result)
		{
			return AAverageConditional(areEqual, A, result);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <param name="to_areEqual">Outgoing message to 'areEqual'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli areEqual, [Fresh] Bernoulli to_areEqual)
		{
			return to_areEqual.GetLogAverageOf(areEqual);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Discrete A, Discrete B)
		{
			Bernoulli to_areEqual = AreEqualAverageConditional(A, B);
			return to_areEqual.GetLogProb(areEqual);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_A">Previous outgoing message to 'A'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Discrete A, Discrete B, [Fresh] Discrete to_A)
		{
			return A.GetLogAverageOf(to_A);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, int A, Discrete B)
		{
			Bernoulli to_areEqual = AreEqualAverageConditional(A, B);
			return to_areEqual.GetLogProb(areEqual);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_B">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, int A, Discrete B, [Fresh] Discrete to_B)
		{
			return B.GetLogAverageOf(to_B);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Discrete A, int B)
		{
			Bernoulli to_areEqual = AreEqualAverageConditional(A, B);
			return to_areEqual.GetLogProb(areEqual);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool areEqual, Discrete A, int B, [Fresh] Discrete to_A)
		{
			return A.GetLogAverageOf(to_A);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(areEqual) p(areEqual) factor(areEqual,a,b) / sum_areEqual p(areEqual) messageTo(areEqual))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli areEqual) { return 0.0; }
		//public static double LogEvidenceRatio(bool areEqual, Discrete A, Discrete B) { return LogAverageFactor(areEqual, A, B); }
		//public static double LogEvidenceRatio(bool areEqual, int A, Discrete B) { return LogAverageFactor(areEqual, A, B); }
		//public static double LogEvidenceRatio(bool areEqual, Discrete A, int B) { return LogAverageFactor(areEqual, A, B); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, Discrete A, Discrete B, [Fresh] Discrete to_A) { return LogAverageFactor(areEqual, A, B, to_A); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="to_B">Outgoing message to 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool areEqual, int A, Discrete B, [Fresh] Discrete to_B) { return LogAverageFactor(areEqual, A, B, to_B); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="to_A">Outgoing message to 'a'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(areEqual,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>	
		public static double LogEvidenceRatio(bool areEqual, Discrete A, int B, [Fresh] Discrete to_A) { return LogAverageFactor(areEqual, A, B, to_A); }

		//-- VMP ----------------------------------------------------------------------------------------

		const string NotSupportedMessage = "Variational Message Passing does not support an AreEqual factor with fixed output.";

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		// same as BP if you use John Winn's rule.
		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(Discrete A, Discrete B)
		{
			return AreEqualAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(int A, Discrete B)
		{
			return AreEqualAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'areEqual'.
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'areEqual' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'areEqual' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(areEqual,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli AreEqualAverageLogarithm(Discrete A, int B)
		{
			return AreEqualAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'areEqual' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_areEqual p(areEqual) factor(areEqual,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete AAverageLogarithm([SkipIfUniform] Bernoulli areEqual, Discrete B, Discrete result)
		{
			if (areEqual.IsPointMass) return AAverageLogarithm(areEqual.Point, B, result);
			if (result == default(Discrete)) result = Discrete.Uniform(B.Dimension, B.Sparsity);
			// when AreEqual is marginalized, the factor is proportional to exp((A==B)*areEqual.LogOdds)
			Vector probs = result.GetWorkspace();
			probs = B.GetProbs(probs);
			probs.SetToFunction(probs, x => Math.Exp(x * areEqual.LogOdds));
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'areEqual' integrated out.
		/// The formula is <c>sum_areEqual p(areEqual) factor(areEqual,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete AAverageLogarithm([SkipIfUniform] Bernoulli areEqual, int B, Discrete result)
		{
			return AAverageConditional(areEqual, B, result);
		}

		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(areEqual,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(DiscreteAreEqualOp.NotSupportedMessage)]
		public static Discrete AAverageLogarithm(bool areEqual, Discrete B, Discrete result)
		{
			if (B.IsPointMass) return AAverageLogarithm(areEqual, B.Point, result);
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'a'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Discrete AAverageLogarithm(bool areEqual, int B, Discrete result)
		{
			return AAverageConditional(areEqual, B, result);
		}

		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'areEqual' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_areEqual p(areEqual) factor(areEqual,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete BAverageLogarithm([SkipIfUniform] Bernoulli areEqual, Discrete A, Discrete result)
		{
			return AAverageLogarithm(areEqual, A, result);
		}

		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Incoming message from 'areEqual'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'areEqual' integrated out.
		/// The formula is <c>sum_areEqual p(areEqual) factor(areEqual,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="areEqual"/> is not a proper distribution</exception>
		public static Discrete BAverageLogarithm([SkipIfUniform] Bernoulli areEqual, int A, Discrete result)
		{
			return AAverageLogarithm(areEqual, A, result);
		}

		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(areEqual,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(DiscreteAreEqualOp.NotSupportedMessage)]
		public static Discrete BAverageLogarithm(bool areEqual, Discrete A, Discrete result)
		{
			return AAverageLogarithm(areEqual, A, result);
		}

		/// <summary>
		/// VMP message to 'b'.
		/// </summary>
		/// <param name="areEqual">Constant value for 'areEqual'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Discrete BAverageLogarithm(bool areEqual, int A, Discrete result)
		{
			return AAverageLogarithm(areEqual, A, result);
		}
	}
}
