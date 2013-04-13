// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;

[assembly: MicrosoftResearch.Infer.Factors.HasMessageFunctions]

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Random{DomainType}"/>, given random arguments to the function.
	/// </summary>
	/// <typeparam name="DomainType">Domain type</typeparam>
	[FactorMethod(typeof(Factor), "Random<>")]
	[Quality(QualityBand.Mature)]
	public static class UnaryOp<DomainType>
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Incoming message from 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		public static double LogAverageFactor<T>(T random, T dist)
			where T : CanGetLogAverageOf<T>
		{
			return dist.GetLogAverageOf(random);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Constant value for 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		public static double LogAverageFactor<T>(DomainType random, T dist)
			where T : CanGetLogProb<DomainType>
		{
			return dist.GetLogProb(random);
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Incoming message from 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		[Skip]
		public static double LogEvidenceRatio<T>(T random, T dist) { return 0.0; }

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Constant value for 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		public static double LogEvidenceRatio<T>(DomainType random, T dist)
			where T : CanGetLogProb<DomainType>
		{
			return LogAverageFactor(random, dist);
		}

		/// <summary>
		/// EP message to 'random'.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="dist">Incoming message from 'dist'</param>
		public static T RandomAverageConditional<T>([IsReturned] T dist) { return dist; }

		//-- VMP ---------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Incoming message from 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		public static double AverageLogFactor<T>(T random, T dist)
			where T : CanGetAverageLog<T>
		{
			return random.GetAverageLog(dist);
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="random">Constant value for 'random'</param>
		/// <param name="dist">Incoming message from 'dist'</param>
		/// <returns></returns>
		public static double AverageLogFactor<T>(DomainType random, T dist)
	where T : CanGetLogProb<DomainType>
		{
			return dist.GetLogProb(random);
		}

		/// <summary>
		/// VMP message to 'random'.
		/// </summary>
		/// <typeparam name="T">Distribution type</typeparam>
		/// <param name="dist">Incoming message from 'dist'</param>
		public static T RandomAverageLogarithm<T>([IsReturned] T dist) { return dist; }

		//-- Max product ---------------------------------------------------------------------------------------

		public static Bernoulli RandomMaxConditional([IsReturned] Bernoulli dist)
		{
			return dist;
		}
		public static UnnormalizedDiscrete RandomMaxConditional([SkipIfUniform] Discrete dist)
		{
			return UnnormalizedDiscrete.FromDiscrete(dist);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Constant{DomainType}"/>, given random arguments to the function.
	/// </summary>
	//[FactorMethod(typeof(Factor), "Constant<>")]
	internal static class ConstantOp<DomainType>
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <returns></returns>
		[Skip]
		public static double LogEvidenceRatio() { return 0.0; }

		/// <summary>
		/// EP message to 'constant'
		/// </summary>
		/// <typeparam name="DistributionType">Distribution type</typeparam>
		/// <param name="value">Constant value</param>
		/// <param name="result">Where to put result</param>
		/// <returns></returns>
		public static DistributionType ConstantAverageConditional<DistributionType>(DomainType value, DistributionType result)
			 where DistributionType : HasPoint<DomainType>
		{
			result.Point = value;
			return result;
		}

		//-- VMP ---------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns></returns>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'constant'
		/// </summary>
		/// <typeparam name="DistributionType">Distribution type</typeparam>
		/// <param name="value">Constant value</param>
		/// <param name="result">Where to put result</param>
		/// <returns></returns>
		public static DistributionType ConstantAverageLogarithm<DistributionType>(DomainType value, DistributionType result)
	 where DistributionType : HasPoint<DomainType>
		{
			result.Point = value;
			return result;
		}
	}
}
