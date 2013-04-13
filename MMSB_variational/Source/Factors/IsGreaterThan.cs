// (C) Copyright 2008-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.IsGreaterThan"/>, given random arguments to the function.
	/// </summary>
	/// <remarks>A and B need not have the same dimension.</remarks>
	[FactorMethod(typeof(Factor), "IsGreaterThan", typeof(int), typeof(int))]
	[Quality(QualityBand.Preview)]
	static public class IsGreaterThanOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool isGreaterThan, int a, int b)
		{
			return (isGreaterThan == Factor.IsGreaterThan(a, b)) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(isGreaterThan,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool isGreaterThan, int a, int b) { return LogAverageFactor(isGreaterThan, a, b); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(isGreaterThan,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool isGreaterThan, int a, int b) { return LogAverageFactor(isGreaterThan, a, b); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isGreaterThan) p(isGreaterThan) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli isGreaterThan, int a, int b)
		{
			return isGreaterThan.GetLogProb(Factor.IsGreaterThan(a, b));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="to_isGreaterThan">Outgoing message to 'isGreaterThan'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isGreaterThan,a) p(isGreaterThan,a) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli isGreaterThan, Discrete a, [Fresh] Bernoulli to_isGreaterThan)
		{
			return to_isGreaterThan.GetLogAverageOf(isGreaterThan);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="to_isGreaterThan">Outgoing message to 'isGreaterThan'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isGreaterThan,b) p(isGreaterThan,b) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Bernoulli isGreaterThan, int a, Discrete b, [Fresh] Bernoulli to_isGreaterThan)
		{
			return to_isGreaterThan.GetLogAverageOf(isGreaterThan);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(isGreaterThan) p(isGreaterThan) factor(isGreaterThan,a,b) / sum_isGreaterThan p(isGreaterThan) messageTo(isGreaterThan))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Bernoulli isGreaterThan) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool isGreaterThan, Discrete a, Discrete b)
		{
			if (isGreaterThan) { // a > b
				double sum = 0.0;
				for (int i = 1; i < a.Dimension; i++) {
					double bSum = 0.0;
					for (int j = 0; (j < b.Dimension) && (j < i); j++) {
						bSum += b[j];
					}
					sum += a[i]*bSum;
				}
				return Math.Log(sum);
			} else { // a <= b
				double sum = 0.0;
				for (int i = 0; i < a.Dimension; i++) {
					double bSum = 0.0;
					for (int j = i; j < b.Dimension; j++) {
						bSum += b[j];
					}
					sum += a[i]*bSum;
				}
				return Math.Log(sum);
			}
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool isGreaterThan, int a, Discrete b)
		{
			if (isGreaterThan) { // a > b
				double sum = 0.0;
				for (int j = 0; (j < b.Dimension) && (j < a); j++) {
					sum += b[j];
				}
				return Math.Log(sum);
			} else { // a <= b
				double sum = 0.0;
				for (int j = a; j < b.Dimension; j++) {
					sum += b[j];
				}
				return Math.Log(sum);
			}
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(isGreaterThan,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(bool isGreaterThan, Discrete a, int b)
		{
			if (isGreaterThan) { // a > b
				double sum = 0.0;
				for (int i = b+1; i < a.Dimension; i++) {
					sum += a[i];
				}
				return Math.Log(sum);
			} else { // a <= b
				double sum = 0.0;
				for (int i = 0; (i < a.Dimension) && (i <= b); i++) {
					sum += a[i];
				}
				return Math.Log(sum);
			}
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(isGreaterThan,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool isGreaterThan, Discrete a, Discrete b)
		{
			return LogAverageFactor(isGreaterThan, a, b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(isGreaterThan,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool isGreaterThan, int a, Discrete b)
		{
			return LogAverageFactor(isGreaterThan, a, b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(isGreaterThan,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(bool isGreaterThan, Discrete a, int b)
		{
			return LogAverageFactor(isGreaterThan, a, b);
		}

		/// <summary>
		/// EP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[p(isGreaterThan) sum_(b) p(b) factor(isGreaterThan,a,b)]/p(isGreaterThan)</c>.
		/// </para></remarks>
		public static Bernoulli IsGreaterThanAverageConditional(int a, Discrete b)
		{
			if (b.IsPointMass) return Bernoulli.PointMass(a > b.Point);
			double sum = 0.0;
			for (int i = 0; (i < a) && (i < b.Dimension); i++) {
				sum += b[i];
			}
			if (sum > 1) sum = 1; // this can happen due to round-off errors
			return new Bernoulli(sum);
		}

		/// <summary>
		/// EP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[p(isGreaterThan) sum_(a) p(a) factor(isGreaterThan,a,b)]/p(isGreaterThan)</c>.
		/// </para></remarks>
		public static Bernoulli IsGreaterThanAverageConditional(Discrete a, int b)
		{
			if (a.IsPointMass) return Bernoulli.PointMass(a.Point > b);
			double sum = 0.0;
			for (int i = b + 1; i < a.Dimension; i++) {
				sum += a[i];
			}
			return new Bernoulli(sum);
		}

		/// <summary>
		/// EP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[p(isGreaterThan) sum_(a,b) p(a,b) factor(isGreaterThan,a,b)]/p(isGreaterThan)</c>.
		/// </para></remarks>
		static public Bernoulli IsGreaterThanAverageConditional(Discrete a, Discrete b)
		{
			if (a.IsPointMass) return IsGreaterThanAverageConditional(a.Point, b);
			if (b.IsPointMass) return IsGreaterThanAverageConditional(a, b.Point);
			double sum = 0.0;
			for (int i = 1; i < a.Dimension; i++) {
				double xsum = 0.0;
				for (int j = 0; (j < i) && (j < b.Dimension); j++) {
					xsum += b[j];
				}
				sum += xsum * a[i];
			}
			return new Bernoulli(sum);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(isGreaterThan,b) p(isGreaterThan,b) factor(isGreaterThan,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		static public Discrete AAverageConditional([SkipIfUniform] Bernoulli isGreaterThan, Discrete b, Discrete result)
		{
			if (b.IsPointMass) return AAverageConditional(isGreaterThan, b.Point, result);
			Vector aProbs = result.GetWorkspace();
			double probTrue = isGreaterThan.GetProbTrue();
			double probFalse = 1 - probTrue;
			for (int i = 0; i < aProbs.Count; i++) {
				double sum1 = 0.0;
				int j = 0;
				for (; (j < i) && (j < b.Dimension); j++) {
					sum1 += b[j];
				}
				double sum0 = 0.0;
				for (; j < b.Dimension; j++) {
					sum0 += b[j];
				}
				aProbs[i] = probTrue*sum1 + probFalse*sum0;
			}
			result.SetProbs(aProbs);
			return result;
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(isGreaterThan) p(isGreaterThan) factor(isGreaterThan,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		static public Discrete AAverageConditional([SkipIfUniform] Bernoulli isGreaterThan, int b, Discrete result)
		{
			double probTrue = isGreaterThan.GetProbTrue();
			double probFalse = 1 - probTrue;
			Vector aProbs = result.GetWorkspace();
			for (int i = 0; i < aProbs.Count; i++) {
				aProbs[i] = (i > b) ? probTrue : probFalse;
			}
			result.SetProbs(aProbs);
			return result;
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(b) p(b) factor(isGreaterThan,a,b)]/p(a)</c>.
		/// </para></remarks>
		static public Discrete AAverageConditional(bool isGreaterThan, Discrete b, Discrete result)
		{
			return AAverageConditional(Bernoulli.PointMass(isGreaterThan), b, result);
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		static public Discrete AAverageConditional(bool isGreaterThan, int b, Discrete result)
		{
			return AAverageConditional(Bernoulli.PointMass(isGreaterThan), b, result);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(isGreaterThan,a) p(isGreaterThan,a) factor(isGreaterThan,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		static public Discrete BAverageConditional([SkipIfUniform] Bernoulli isGreaterThan, Discrete a, Discrete result)
		{
			if (a.IsPointMass) return BAverageConditional(isGreaterThan, a.Point, result);
			Vector bProbs = result.GetWorkspace();
			double probTrue = isGreaterThan.GetProbTrue();
			double probFalse = 1 - probTrue;
			for (int j = 0; j < bProbs.Count; j++) {
				double sum0 = 0.0;
				int i = 0;
				for (; (i <= j) && (i < a.Dimension); i++) {
					sum0 += a[i];
				}
				double sum1 = 0.0;
				for (; i < a.Dimension; i++) {
					sum1 += a[i];
				}
				bProbs[j] = probTrue*sum1 + probFalse*sum0;
			}
			result.SetProbs(bProbs);
			return result;
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(isGreaterThan) p(isGreaterThan) factor(isGreaterThan,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		static public Discrete BAverageConditional([SkipIfUniform] Bernoulli isGreaterThan, int a, Discrete result)
		{
			Vector bProbs = result.GetWorkspace();
			double probTrue = isGreaterThan.GetProbTrue();
			double probFalse = 1 - probTrue;
			for (int j = 0; j < bProbs.Count; j++) {
				bProbs[j] = (a > j) ? probTrue : probFalse;
			}
			result.SetProbs(bProbs);
			return result;
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(a) p(a) factor(isGreaterThan,a,b)]/p(b)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(bool isGreaterThan, Discrete a, Discrete result)
		{
			return BAverageConditional(Bernoulli.PointMass(isGreaterThan), a, result);
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Discrete BAverageConditional(bool isGreaterThan, int a, Discrete result)
		{
			return BAverageConditional(Bernoulli.PointMass(isGreaterThan), a, result);
		}

		// VMP /////////////////////////////////////////////////////////////////////////////////////////////////////////
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Bernoulli isGreaterThan, Discrete a, Discrete b) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Bernoulli isGreaterThan, int a, Discrete b) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Bernoulli isGreaterThan, Discrete a, int b) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Bernoulli isGreaterThan, int a, int b) { return 0.0; }

		/// <summary>
		/// VMP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(isGreaterThan,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli IsGreaterThanAverageLogarithm(int a, Discrete b) { return IsGreaterThanAverageConditional(a, b); }
		/// <summary>
		/// VMP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(isGreaterThan,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli IsGreaterThanAverageLogarithm(Discrete a, int b) { return IsGreaterThanAverageConditional(a, b); }
		/// <summary>
		/// VMP message to 'isGreaterThan'
		/// </summary>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'isGreaterThan' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'isGreaterThan' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(isGreaterThan,a,b)]</c>.
		/// </para></remarks>
		public static Bernoulli IsGreaterThanAverageLogarithm(Discrete a, Discrete b) { return IsGreaterThanAverageConditional(a, b); }
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'isGreaterThan' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_isGreaterThan p(isGreaterThan) factor(isGreaterThan,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		public static Discrete AAverageLogarithm([SkipIfUniform] Bernoulli isGreaterThan, Discrete b, Discrete result)
		{
			if (b.IsPointMass) return AAverageLogarithm(isGreaterThan, b.Point, result);
			if (isGreaterThan.IsPointMass) return AAverageLogarithm(isGreaterThan.Point, b, result);
			// f(a,b) = p(c=1) I(a > b) + p(c=0) I(a <= b)
			// message to a = exp(sum_b q(b) log f(a,b))
			Vector aProbs = result.GetWorkspace();
			double logProbTrue = isGreaterThan.GetLogProbTrue();
			double logProbFalse = isGreaterThan.GetLogProbFalse();
			for (int i = 0; i < aProbs.Count; i++) {
				double sum = 0.0;
				int j = 0;
				for (; (j < i) && (j < b.Dimension); j++) {
					sum += logProbTrue*b[j];
				}
				for (; j < b.Dimension; j++) {
					sum += logProbFalse*b[j];
				}
				aProbs[i] = Math.Exp(sum);
			}
			result.SetProbs(aProbs);
			return result;
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'isGreaterThan' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_isGreaterThan p(isGreaterThan) factor(isGreaterThan,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		public static Discrete BAverageLogarithm([SkipIfUniform] Bernoulli isGreaterThan, Discrete a, Discrete result)
		{
			if (a.IsPointMass) return BAverageLogarithm(isGreaterThan, a.Point, result);
			if (isGreaterThan.IsPointMass) return BAverageLogarithm(isGreaterThan.Point, a, result);
			// f(a,b) = p(c=1) I(a > b) + p(c=0) I(a <= b)
			// message to b = exp(sum_a q(a) log f(a,b))
			Vector bProbs = result.GetWorkspace();
			double logProbTrue = isGreaterThan.GetLogProbTrue();
			double logProbFalse = isGreaterThan.GetLogProbFalse();
			for (int j = 0; j < bProbs.Count; j++) {
				double sum = 0.0;
				int i = 0;
				for (; (i <= j) && (i < a.Dimension); i++) {
					sum += logProbFalse*a[i];
				}
				for (; i < a.Dimension; i++) {
					sum += logProbTrue*a[i];
				}
				bProbs[j] = Math.Exp(sum);
			}
			result.SetProbs(bProbs);
			return result;
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'isGreaterThan' integrated out.
		/// The formula is <c>sum_isGreaterThan p(isGreaterThan) factor(isGreaterThan,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		public static Discrete AAverageLogarithm([SkipIfUniform] Bernoulli isGreaterThan, int b, Discrete result) { return AAverageConditional(isGreaterThan, b, result); }
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Incoming message from 'isGreaterThan'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'isGreaterThan' integrated out.
		/// The formula is <c>sum_isGreaterThan p(isGreaterThan) factor(isGreaterThan,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="isGreaterThan"/> is not a proper distribution</exception>
		public static Discrete BAverageLogarithm([SkipIfUniform] Bernoulli isGreaterThan, int a, Discrete result) { return BAverageConditional(isGreaterThan, a, result); }

		const string NotSupportedMessage = "Variational Message Passing does not support an IsGreaterThan factor with fixed output.";
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static double AverageLogFactor(bool isGreaterThan, Discrete a, Discrete b)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool isGreaterThan, int a, Discrete b)
		{
			return LogAverageFactor(isGreaterThan, a, b);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(bool isGreaterThan, Discrete a, int b)
		{
			return LogAverageFactor(isGreaterThan, a, b);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(isGreaterThan,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static Discrete AAverageLogarithm(bool isGreaterThan, Discrete b, Discrete result) 
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(isGreaterThan,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static Discrete BAverageLogarithm(bool isGreaterThan, Discrete a, Discrete result) 
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Discrete AAverageLogarithm(bool isGreaterThan, int b, Discrete result) { return AAverageConditional(isGreaterThan, b, result); }
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="isGreaterThan">Constant value for 'isGreaterThan'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Discrete BAverageLogarithm(bool isGreaterThan, int a, Discrete result) { return BAverageConditional(isGreaterThan, a, result); }
	}
}
