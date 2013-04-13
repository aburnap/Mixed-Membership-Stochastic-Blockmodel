// (C) Copyright 2009 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Internal factors relating to Gibbs message passing
	/// </summary>
	public static class GibbsFactor
	{
		/// <summary>
		/// Accumulate factor - for accumulating samples
		/// </summary>
		/// <typeparam name="TDist">The type of the distribution</typeparam>
		/// <typeparam name="T">The type of the domain</typeparam>
		/// <typeparam name="TEst">The type of the estimator</typeparam>
		/// <param name="item">The item to accumulate</param>
		/// <param name="estimator">The estimator</param>
		/// <returns></returns>
		[Stochastic]
		public static T Accumulate<TDist, T, TEst>(T item, TEst estimator)
			where TDist : IDistribution<T>, Sampleable<T>
			where TEst : Estimator<TDist>, Accumulator<T>
		{
			return item;
		}
	}

	/// <summary>
	/// Accumulate operator - used for accumulating samples for a Gibbs algorithm
	/// </summary>
	/// <typeparam name="TDist">Distribution type</typeparam>
	/// <typeparam name="T">Domian type</typeparam>
	/// <typeparam name="TEst">Estimator type</typeparam>
	[FactorMethod(new string[] { "Estimate", "Sample", "Estimator" }, typeof(GibbsFactor), "Accumulate<,,>")]
	public static class AccumulateOp<TDist, T, TEst>
		where TDist : IDistribution<T>, Sampleable<T>
		where TEst : Estimator<TDist>, Accumulator<T>
	{

		/// <summary>
		/// EP message to estimate. The incoming sample is added
		/// to the accumulator, and the resulting estimated distribution
		/// is returned
		/// </summary>
		/// <param name="sample">The sample</param>
		/// <param name="estimator">The estimator and accumulator</param>
		/// <param name="result">The result</param>
		/// <returns></returns>
		public static TDist EstimateAverageConditional(
			T sample,
			TEst estimator,
			TDist result)
		{
			estimator.Add(sample);
			return estimator.GetDistribution(result);
		}

		/// <summary>
		/// EP message to sample
		/// </summary>
		/// <param name="sample">The sample</param>
		/// <returns></returns>
		public static TDist SampleAverageConditional(
			TDist result)
		{
			// return estimator.GetDistribution(result);
			result.SetToUniform();
			return result;
		}

		/// <summary>
		/// EP message to estimate. The incoming sample is added
		/// to the accumulator, and the resulting estimated distribution
		/// is returned
		/// </summary>
		/// <param name="sample">The sample</param>
		/// <param name="estimator">The estimator and accumulator</param>
		/// <param name="result">The result</param>
		/// <returns></returns>
		public static TDist EstimateAverageLogarithm(
			T sample,
			TEst estimator,
			TDist result)
		{
			return EstimateAverageConditional(sample, estimator, result);
		}

		/// <summary>
		/// VMP message to sample
		/// </summary>
		/// <param name="result">The result</param>
		/// <returns></returns>
		public static TDist SampleAverageLogarithm(TDist result)
		{
			return SampleAverageConditional(result);
		}
	}
}
