// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Plus(int,int)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Plus", typeof(int), typeof(int))]
	[Quality(QualityBand.Stable)]
	public static class IntegerPlusOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum) p(Sum) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sum, int a, int b)
		{
			return sum.GetLogProb(Factor.Plus(a, b));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum,A) p(Sum,A) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sum, Discrete a, [Fresh] Discrete to_sum)
		{
			return to_sum.GetLogAverageOf(sum);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="to_sum">Outgoing message to 'sum'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum,B) p(Sum,B) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Discrete sum, int a, Discrete b, [Fresh] Discrete to_sum)
		{
			return to_sum.GetLogAverageOf(sum);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(Sum) p(Sum) factor(Sum,A,B) / sum_Sum p(Sum) messageTo(Sum))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Discrete sum) { return 0.0; }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A,B) p(A,B) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sum, Discrete a, Discrete b)
		{
			if (a.IsPointMass) return LogAverageFactor(sum, a.Point, b);
			double z = 0.0;
			for (int i = 0; (i < a.Dimension) && (sum-i < b.Dimension); i++)
			{
				z += a[i] * b[sum-i];
			}
			return Math.Log(z);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sum, int a, Discrete b)
		{
			if (b.IsPointMass) return LogAverageFactor(sum, a, b.Point);
			int j = sum-a;
			if (j < b.Dimension) return b.GetLogProb(j);
			else return Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sum, Discrete a, int b)
		{
			if (a.IsPointMass) return LogAverageFactor(sum, a.Point, b);
			int i = sum-b;
			if (i < a.Dimension) return a.GetLogProb(i);
			else return Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sum,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(int sum, int a, int b)
		{
			return (sum == a+b) ? 0.0 : Double.NegativeInfinity;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A,B) p(A,B) factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sum, Discrete a, Discrete b) { return LogAverageFactor(sum, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sum, int a, Discrete b) { return LogAverageFactor(sum, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sum, Discrete a, int b) { return LogAverageFactor(sum, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(Sum,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(int sum, int a, int b) { return LogAverageFactor(sum, a, b); }

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(A,B) p(A,B) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		public static Discrete SumAverageConditional(Discrete a, Discrete b, Discrete result)
		{
			if (a.IsPointMass) return SumAverageConditional(a.Point, b, result);
			Vector probs = result.GetWorkspace();
			probs.SetAllElementsTo(0.0);
			for (int i = 0; i < a.Dimension; i++)
			{
				double pa = a[i];
				for (int j = 0; j < b.Dimension; j++)
				{
					probs[i+j] += (pa * b[j]);
				}
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(B) p(B) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		public static Discrete SumAverageConditional(int a, Discrete b, Discrete result)
		{
			Vector probs = result.GetWorkspace();
			probs.SetAllElementsTo(0.0);
			for (int j = 0; j < b.Dimension; j++)
			{
				probs[a+j] = b[j];
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'Sum'
		/// </summary>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'Sum' as the random arguments are varied.
		/// The formula is <c>proj[p(Sum) sum_(A) p(A) factor(Sum,A,B)]/p(Sum)</c>.
		/// </para></remarks>
		public static Discrete SumAverageConditional(Discrete a, int b, Discrete result)
		{
			return SumAverageConditional(b, a, result);
		}

		[Skip]
		public static Discrete SumAverageConditionalInit([IgnoreDependency] Discrete a, [IgnoreDependency] Discrete b)
		{
			return Discrete.Uniform(a.Dimension+b.Dimension-1);
		}
		[Skip]
		public static Discrete SumAverageConditionalInit([IgnoreDependency] Discrete a, [IgnoreDependency] int b)
		{
			return Discrete.Uniform(a.Dimension+b);
		}
		[Skip]
		public static Discrete SumAverageConditionalInit([IgnoreDependency] int a, [IgnoreDependency] Discrete b)
		{
			return Discrete.Uniform(a+b.Dimension);
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(Sum,B) p(Sum,B) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		public static Discrete AAverageConditional(Discrete sum, Discrete b, Discrete result)
		{
			// message to a = sum_b p(sum = a+b) p(b)
			Vector probs = result.GetWorkspace();
			for (int i = 0; i < result.Dimension; i++)
			{
				double p = 0.0;
				int maxPlus1 = sum.Dimension - i;
				for (int j = 0; (j < b.Dimension) && (j < maxPlus1); j++)
				{
					p += b[j] * sum[i+j];
				}
				probs[i] = p;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(Sum) p(Sum) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		public static Discrete AAverageConditional(Discrete sum, int b, Discrete result)
		{
			// message to a = p(sum = a+b)
			Vector probs = result.GetWorkspace();
			int maxPlus1 = sum.Dimension - b;
			for (int i = 0; (i < result.Dimension) && (i < maxPlus1); i++)
			{
				probs[i] = sum[i+b];
			}
			for (int i = maxPlus1; i < result.Dimension; i++)
			{
				probs[i] = 0.0;
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(B) p(B) factor(Sum,A,B)]/p(A)</c>.
		/// </para></remarks>
		public static Discrete AAverageConditional(int sum, Discrete b, Discrete result)
		{
			// message to a = p(b = sum-a)
			Vector probs = result.GetWorkspace();
			int min = sum - (b.Dimension - 1);
			for (int i = 0; i < min; i++)
			{
				probs[i] = 0.0;
			}
			for (int i = min; (i < result.Dimension) && (i <= sum); i++)
			{
				probs[i] = b[sum-i];
			}
			result.SetProbs(probs);
			return result;
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		public static Discrete AAverageConditional(int sum, int b, Discrete result)
		{
			int a = sum-b;
			if (a < 0 || a > result.Dimension) throw new AllZeroException();
			result.Point = a;
			return result;
		}

		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(Sum,A) p(Sum,A) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(Discrete sum, Discrete a, Discrete result)
		{
			return AAverageConditional(sum, a, result);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="sum">Incoming message from 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(Sum) p(Sum) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(Discrete sum, int a, Discrete result)
		{
			return AAverageConditional(sum, a, result);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(A) p(A) factor(Sum,A,B)]/p(B)</c>.
		/// </para></remarks>
		public static Discrete BAverageConditional(int sum, Discrete a, Discrete result)
		{
			return AAverageConditional(sum, a, result);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="sum">Constant value for 'Sum'.</param>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		public static Discrete BAverageConditional(int sum, int a, Discrete result)
		{
			return AAverageConditional(sum, a, result);
		}

		// Poisson case -----------------------------------------------------------------------------------------

		public static Poisson SumAverageConditional(Poisson a, Poisson b)
		{
			if(a.Precision != 1 || b.Precision != 1) throw new NotImplementedException("precision != 1 not implemented");
			return new Poisson(a.Rate + b.Rate);
		}

		public static Poisson AAverageConditional(int sum, Poisson a, Poisson b)
		{
			// a has rate r
			// b has rate t
			// p(a) =propto r^a/a! t^(s-a)/(s-a)! =propto Binom(s, r/(r+t))
			// E[a] = s*r/(r+t)
			if (a.Precision != 1 || b.Precision != 1) throw new NotImplementedException("precision != 1 not implemented");
			double rate = a.Rate/(a.Rate + b.Rate);
			return new Poisson(sum*rate) / a;
		}

		public static Poisson BAverageConditional(int sum, Poisson a, Poisson b)
		{
			return AAverageConditional(sum, b, a);
		}

		public static Poisson AAverageConditional(Poisson sum, Poisson a, Poisson b)
		{
			if (sum.IsPointMass) return AAverageConditional(sum.Point, a, b);
			// p(a) = r^a/a! sum_(s >= a) q^s t^(s-a)/(s-a)!
			//      = r^a/a! q^a sum_(s >= a) q^(s-a) t^(s-a)/(s-a)!
			//      = r^a/a! q^a exp(q t)
			if (sum.Precision != 0) throw new NotImplementedException("sum.Precision != 0 not implemented");
			if (b.Precision != 1) throw new NotImplementedException("b.Precision != 1 not implemented");
			return new Poisson(sum.Rate, 0);
		}

		public static Poisson BAverageConditional(Poisson sum, Poisson a, Poisson b)
		{
			return AAverageConditional(sum, b, a);
		}
	}
}
