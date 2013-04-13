// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// A mixture of distributions of the same type 
	/// </summary>
	/// <typeparam name="T">The distribution type</typeparam>
	public class Mixture<T> : Estimator<Mixture<T>>, Accumulator<T>
	{
		/// <summary>
		/// The components
		/// </summary>
		public List<T> Components;
		/// <summary>
		/// The mixing weight of each component.  Does not necessarily sum to 1.
		/// </summary>
		public List<double> Weights;

		/// <summary>
		/// Add a component to the mixture with a given weight
		/// </summary>
		/// <param name="item">The component to add</param>
		/// <param name="weight">The weight</param>
		public void Add(T item, double weight)
		{
			Components.Add(item);
			Weights.Add(weight);
		}

		/// <summary>
		/// Add a component to the mixture. A weight of 1 is assumed
		/// </summary>
		/// <param name="item">The component to add</param>
		public void Add(T item)
		{
			Add(item, 1.0);
		}

		/// <summary>
		/// The sum of the component weights
		/// </summary>
		/// <returns></returns>
		public double WeightSum()
		{
			double sum = 0.0;
			foreach (double w in Weights) {
				sum += w;
			}
			return sum;
		}

		/// <summary>
		/// Normalize the weights to add to 1
		/// </summary>
		public void Normalize()
		{
			double sum = WeightSum();
			if (sum > 0) {
				for (int i = 0; i < Weights.Count; i++) {
					Weights[i] /= sum;
				}
			}
		}

		/// <summary>
		/// The the resulting mixture
		/// </summary>
		/// <param name="result">Where to put the resulting mixture</param>
		/// <returns></returns>
		public Mixture<T> GetDistribution(Mixture<T> result)
		{
			Normalize();
			return this;
		}

		/// <summary>
		/// Create a mixture model
		/// </summary>
		public Mixture()
		{
			Components = new List<T>();
			Weights = new List<double>();
		}

		/// <summary>
		/// Clears the estimator
		/// </summary>
		public void Clear()
		{
			Components.Clear();
			Weights.Clear();
		}
	}
	public class Mixture<T, DomainType> : Mixture<T>, CanGetAverageLog<T>, CanGetLogProb<DomainType>, Sampleable<DomainType>
		where T : CanGetAverageLog<T>, CanGetLogProb<DomainType>, Sampleable<DomainType>
	{
		public double GetAverageLog(T that)
		{
			double sum = 0;
			double weightSum = 0;
			for (int i = 0; i < Weights.Count; i++) {
				weightSum += Weights[i];
				sum += Weights[i]*Components[i].GetAverageLog(that);
			}
			return sum/weightSum;
		}

		public double GetLogProb(DomainType value)
		{
			double logProb = double.NegativeInfinity;
			double weightSum = 0;
			for (int i = 0; i < Weights.Count; i++) {
				weightSum += Weights[i];
				logProb = MMath.LogSumExp(logProb, Math.Log(Weights[i]) + Components[i].GetLogProb(value));
			}
			return logProb - Math.Log(weightSum);
		}

		public DomainType Sample()
		{
			int i = Rand.Sample(Weights, WeightSum());
			return Components[i].Sample();
		}

		public DomainType Sample(DomainType result)
		{
			int i = Rand.Sample(Weights, WeightSum());
			return Components[i].Sample(result);
		}
	}
}
