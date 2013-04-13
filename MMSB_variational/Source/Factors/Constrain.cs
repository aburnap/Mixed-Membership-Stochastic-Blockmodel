// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Exception which is thrown when a constraint is violated.  This
	/// occurs when an observation does not hold true or a weight is 0.
	/// </summary>
	public class ConstraintViolatedException : Exception
	{
		/// <summary>
		/// Constructs a constraint violation expception with a specified error message 
		/// </summary>
		/// <param name="s"></param>
		public ConstraintViolatedException(string s)
			: base(s)
		{
		}
	}

	/// <summary>
	/// A repository of commonly used constraint methods.
	/// </summary>
	public static class Constrain
	{
		/// <summary>
		/// Constrains a value to be equal to a sample from dist.
		/// </summary>
		/// <typeparam name="DomainType">Domain type</typeparam>
		/// <typeparam name="DistributionType">Distribution type</typeparam>
		/// <param name="value">Value</param>
		/// <param name="dist">Distribution instance</param>
		[Stochastic]
		public static void EqualRandom<DomainType, DistributionType>(DomainType value, DistributionType dist)
			where DistributionType : Sampleable<DomainType>
		{
			if (!value.Equals(dist.Sample())) throw new ConstraintViolatedException(value + " != " + dist);
		}

		/// <summary>
		/// Constrains a value to be equal to another value.
		/// </summary>
		/// <typeparam name="T">Value type</typeparam>
		/// <param name="A">First value</param>
		/// <param name="B">Second value</param>
		public static void Equal<T>(T A, T B)
		{
			if (!A.Equals(B)) throw new ConstraintViolatedException(A + " != " + B);
		}

		/// <summary>
		/// Constrains a set of integers to contain a particular integer.
		/// </summary>
		/// <param name="set">The set of integers, specified as a list</param>
		/// <param name="i">The integer which the set must contain</param>
		[Hidden]
		public static void Contain(IList<int> set, int i)
		{
			if (!set.Contains(i))
			{
				throw new ConstraintViolatedException(
					"Containment constraint violated (the supplied set did not contain the integer "+i+")");
			}
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Constrain.EqualRandom{DomainType,DistributionType}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Constrain), "EqualRandom<,>")]
	[Quality(QualityBand.Mature)]
	public static class ConstrainEqualRandomOp<DomainType>
	//		where DistributionType : Sampleable<DomainType>,  // from the factor definition
	//	CanGetLogAverageOf<DistributionType>, CanGetAverageLog<DistributionType>, CanGetLogProb<DomainType>
	{
		/// <summary>
		/// EP message to 'value'.
		/// </summary>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns>The outgoing EP message to the 'value' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int f(value,x) q(x) dx</c> where <c>x = (dist)</c>.
		/// </para></remarks>
		public static DistributionType ValueAverageConditional<DistributionType>([IsReturned] DistributionType dist)
		{
			return dist;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="value">Incoming message from 'value'.</param>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (value,dist)</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType>(DistributionType value, DistributionType dist)
			where DistributionType : CanGetLogAverageOf<DistributionType>
		{
			return value.GetLogAverageOf(dist);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="value">Constant value for 'value'.</param>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns><c>log(int f(x) qnotf(x) dx)</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(int f(x) qnotf(x) dx)</c>
		/// where <c>x = (value,dist)</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType>(DomainType value, [Proper] DistributionType dist)
			where DistributionType : CanGetLogProb<DomainType>
		{
			return dist.GetLogProb(value);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		public static double LogEvidenceRatio<DistributionType>(DistributionType value, DistributionType dist)
			where DistributionType : CanGetLogAverageOf<DistributionType>
		{
			return LogAverageFactor(value, dist);
		}
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		public static double LogEvidenceRatio<DistributionType>(DomainType value, DistributionType dist)
			where DistributionType : CanGetLogProb<DomainType>
		{
			return LogAverageFactor(value, dist);
		}

		//-- VMP --------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="value">Incoming message from 'value'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (value,dist)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="value"/> is not a proper distribution</exception>
		public static double AverageLogFactor<DistributionType>(DistributionType value, DistributionType dist)
			where DistributionType : CanGetAverageLog<DistributionType>
		{
			return value.GetAverageLog(dist);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="value">Constant value for 'value'.</param>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// where <c>x = (value,dist)</c>.
		/// </para></remarks>
		public static double AverageLogFactor<DistributionType>(DomainType value, [Proper] DistributionType dist)
			where DistributionType : CanGetLogProb<DomainType>
		{
			return dist.GetLogProb(value);
		}

		/// <summary>
		/// VMP message to 'value'.
		/// </summary>
		/// <param name="dist">Constant value for 'dist'.</param>
		/// <returns>The outgoing VMP message to the 'value' argument.</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'value'.
		/// The formula is <c>int log(f(value,x)) q(x) dx</c> where <c>x = (dist)</c>.
		/// </para></remarks>
		public static DistributionType ValueAverageLogarithm<DistributionType>([IsReturned] DistributionType dist)
		{
			return dist;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Constrain.Equal{T}"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Constrain), "Equal<>")]
	[Quality(QualityBand.Mature)]
	public static class ConstrainEqualOp<T>
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A,B) p(A,B) factor(A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor<DistributionType>(DistributionType a, DistributionType b)
			where DistributionType : CanGetLogAverageOf<DistributionType>
		{
			return a.GetLogAverageOf(b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(A,B))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogAverageFactor<DistributionType>(T a, DistributionType b)
			where DistributionType : CanGetLogProb<T>
		{
			return b.GetLogProb(a);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Incoming message from 'A'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(A,B))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogAverageFactor<DistributionType>(DistributionType a, T b)
			where DistributionType : CanGetLogProb<T>
		{
			return a.GetLogProb(b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(T a, T b)
		{
			IEqualityComparer<T> equalityComparer = Utils.Util.GetEqualityComparer<T>();
			return (equalityComparer.Equals(a,b) ? 0.0 : Double.NegativeInfinity);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A,B) p(A,B) factor(A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType>(DistributionType a, DistributionType b)
			where DistributionType : CanGetLogAverageOf<DistributionType>
		{
			return LogAverageFactor(a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType>(T a, DistributionType b)
			where DistributionType : CanGetLogProb<T>
		{
			return LogAverageFactor(a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Incoming message from 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio<DistributionType>(DistributionType a, T b)
			where DistributionType : CanGetLogProb<T>
		{
			return LogAverageFactor(a, b);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="a">Constant value for 'A'.</param>
		/// <param name="b">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(T a, T b)
		{
			return LogAverageFactor(a, b);
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(B) p(B) factor(A,B)]/p(A)</c>.
		/// </para></remarks>
		public static DistributionType AAverageConditional<DistributionType>([IsReturned] DistributionType B, DistributionType result)
	where DistributionType : SettableTo<DistributionType>
		{
			result.SetTo(B);
			return result;
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		public static DistributionType AAverageConditional<DistributionType>(T B, DistributionType result)
	where DistributionType : HasPoint<T>
		{
			result.Point = B;
			return result;
		}

		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(A) p(A) factor(A,B)]/p(B)</c>.
		/// </para></remarks>
		public static DistributionType BAverageConditional<DistributionType>([IsReturned] DistributionType A, DistributionType result)
	where DistributionType : SettableTo<DistributionType>
		{
			result.SetTo(A);
			return result;
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		public static DistributionType BAverageConditional<DistributionType>(T A, DistributionType result)
	where DistributionType : HasPoint<T>
		{
			return AAverageConditional<DistributionType>(A, result);
		}

		//-- VMP -----------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		const string NotSupportedMessage = "VMP does not support Constrain.Equal between random variables";

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'A'.
		/// The formula is <c>exp(sum_(B) p(B) log(factor(A,B)))</c>.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static DistributionType AAverageLogarithm<DistributionType>(DistributionType B, DistributionType result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'B'.
		/// The formula is <c>exp(sum_(A) p(A) log(factor(A,B)))</c>.
		/// </para></remarks>
		[NotSupported(NotSupportedMessage)]
		public static DistributionType BAverageLogarithm<DistributionType>(DistributionType A, DistributionType result)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		public static DistributionType AAverageLogarithm<DistributionType>(T B, DistributionType result)
			where DistributionType : HasPoint<T>
		{
			result.Point = B;
			return result;
		}

		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		public static DistributionType BAverageLogarithm<DistributionType>(T A, DistributionType result)
			where DistributionType : HasPoint<T>
		{
			return AAverageLogarithm<DistributionType>(A, result);
		}

		//-- Max product ----------------------------------------------------------------------
		/// <summary>
		/// Max product message to 'A' 
		/// </summary>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static DistributionType AMaxConditional<DistributionType>([IsReturned] DistributionType B, DistributionType result)
			where DistributionType : SettableTo<DistributionType>
		{
			result.SetTo(B);
			return result;
		}

		/// <summary>
		/// Max product message to 'B'
		/// </summary>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static DistributionType BMaxConditional<DistributionType>([IsReturned] DistributionType A, DistributionType result)
			where DistributionType : SettableTo<DistributionType>
		{
			result.SetTo(A);
			return result;
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Constrain.Contain"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Constrain), "Contain")]
	[Quality(QualityBand.Experimental)]
	public class ConstrainContainOp
	{
		/// <summary>
		/// VMP message to 'set'
		/// </summary>
		/// <param name="i">Constant value for 'i'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'set' conditioned on the given values.
		/// </para></remarks>
		public static BernoulliIntegerSubset SetAverageLogarithm(int i, BernoulliIntegerSubset result)
		{
			result.SetToUniform();
			result.LogOddsVector[i] = double.PositiveInfinity;
			return result;
		}
	}

}
