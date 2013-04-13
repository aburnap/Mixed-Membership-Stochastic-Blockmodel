// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Factors;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// A discrete distribution over the values of an enum.
	/// </summary>
	/// <typeparam name="TEnum"></typeparam>
	public class DiscreteEnum<TEnum> : Discrete, IDistribution<TEnum>, SettableTo<DiscreteEnum<TEnum>>, SettableToProduct<DiscreteEnum<TEnum>>,
		SettableToRatio<DiscreteEnum<TEnum>>, SettableToPower<DiscreteEnum<TEnum>>, SettableToWeightedSum<DiscreteEnum<TEnum>>, 
		CanGetLogAverageOf<DiscreteEnum<TEnum>>, CanGetLogAverageOfPower<DiscreteEnum<TEnum>>,
		CanGetAverageLog<DiscreteEnum<TEnum>>, Sampleable<TEnum>
	{
		/// <summary>
		/// Creates a uniform distribution over the enum values.
		/// </summary>
		public DiscreteEnum()
			: base(Enum.GetValues(typeof(TEnum)).Length)
		{
		}

		/// <summary>
		/// Creates a distribution over the enum values using the probabilities from the given discrete distribution.
		/// </summary>
		/// <param name="d"></param>
		public DiscreteEnum(Discrete d)
			: this(d.GetProbs())
		{
		}

		/// <summary>
		/// Creates a discrete distribution over enum values from the given probabilities.
		/// </summary>
		/// <param name="probs"></param>
		public DiscreteEnum(params double[] probs)
			: base(probs)
		{
			int len = Enum.GetValues(typeof(TEnum)).Length;
			if (probs.Length!=len) throw new ArgumentException("Invalid length of probability vector, was "+probs.Length+", should be "+len);
		}

		/// <summary>
		/// Creates a discrete distribution over enum values from the given vector of probabilities.
		/// </summary>
		/// <param name="probs"></param>
		[Construction("GetProbs")]
		public DiscreteEnum(Vector probs)
			: base(probs)
		{
			int len = Enum.GetValues(typeof(TEnum)).Length;
			if (probs.Count!=len) throw new ArgumentException("Invalid length of probability vector, was "+probs.Count+", should be "+len);
		}

		/// <summary>
		/// Creates a uniform distribution over the enum values.
		/// </summary>
		[Skip]
		public static DiscreteEnum<TEnum> Uniform()
		{
			return new DiscreteEnum<TEnum>();
		}

		/// <summary>
		/// Creates a Discrete distribution which allows only one enum value.
		/// </summary>
		/// <param name="value">The allowed value.</param>
		/// <returns></returns>
		public static DiscreteEnum<TEnum> PointMass(TEnum value)
		{
			var d = DiscreteEnum<TEnum>.Uniform();
			((HasPoint<TEnum>)d).Point = value;
			return d;
		}

		object ICloneable.Clone()
		{
			DiscreteEnum<TEnum> clone = new DiscreteEnum<TEnum>();
			clone.SetProbs(GetProbs());
			return clone;
		}


		TEnum HasPoint<TEnum>.Point
		{
			get
			{
				Array values = Enum.GetValues(typeof(TEnum));
				return (TEnum)values.GetValue(base.Point);
			}
			set
			{
				base.Point = (int)(object)value;
			}
		}




		double CanGetLogProb<TEnum>.GetLogProb(TEnum value)
		{
			return base.GetLogProb((int)(object)value);
		}

		/// <summary>
		/// Sets the parameters of this instance to the parameters of that instance
		/// </summary>
		/// <param name="value">That instance</param>
		public void SetTo(DiscreteEnum<TEnum> value) { base.SetTo(value); }
		/// <summary>
		/// Sets the parameters to represent the product of two discrete enum distributions.
		/// </summary>
		/// <param name="a">The first discrete enum distribution</param>
		/// <param name="b">The second discrete enum distribution</param>
		public void SetToProduct(DiscreteEnum<TEnum> a, DiscreteEnum<TEnum> b) { base.SetToProduct(a, b); }
		/// <summary>
		/// Sets the parameters to represent the ratio of two discrete enum distributions.
		/// </summary>
		/// <param name="numerator">The numerator discrete enum distribution</param>
		/// <param name="denominator">The denominator discrete enum distribution</param>
		public void SetToRatio(DiscreteEnum<TEnum> numerator, DiscreteEnum<TEnum> denominator) { base.SetToRatio(numerator, denominator); }

		/// <summary>
		/// Sets the parameters to represent the power of a discrete enum distributions.
		/// </summary>
		/// <param name="value">The discrete enum distribution</param>
		/// <param name="exponent">The exponent</param>
		public void SetToPower(DiscreteEnum<TEnum> value, double exponent) { base.SetToPower(value, exponent); }

		/// <summary>
		/// Sets the parameters to represent the weighted sum of two discrete enum distributions.
		/// </summary>
		/// <param name="value1">The first discrete enum distribution</param>
		/// <param name="weight1">The first weight</param>
		/// <param name="value2">The second discrete enum distribution</param>
		/// <param name="weight2">The second weight</param>
		public void SetToSum(double weight1, DiscreteEnum<TEnum> value1, double weight2, DiscreteEnum<TEnum> value2) { base.SetToSum(weight1, value1, weight2, value2); }

		/// <summary>
		/// The log of the integral of the product of this discrete enum and that discrete enum
		/// </summary>
		/// <param name="that">That discrete enum distribution</param>
		/// <returns>The log inner product</returns>
		public double GetLogAverageOf(DiscreteEnum<TEnum> that) { return base.GetLogAverageOf(that); }

		/// <summary>
		/// Get the integral of this distribution times another distribution raised to a power.
		/// </summary>
		/// <param name="that"></param>
		/// <param name="power"></param>
		/// <returns></returns>
		public double GetLogAverageOfPower(DiscreteEnum<TEnum> that, double power) { return base.GetLogAverageOfPower(that, power); }

		/// <summary>
		/// The expected logarithm of that distribution under this distribution.
		/// </summary>
		/// <param name="that">The distribution to take the logarithm of.</param>
		/// <returns><c>sum_x this.Evaluate(x)*Math.Log(that.Evaluate(x))</c></returns>
		/// <remarks>This is also known as the cross entropy.</remarks>
		public double GetAverageLog(DiscreteEnum<TEnum> that) { return base.GetAverageLog(that); }

		/// <summary>
		/// Returns a string representation of this discrete enum distribution.
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			StringBuilder sb = new StringBuilder();
			string[] names = Enum.GetNames(typeof(TEnum));
			int maxLen = 0; foreach(string s in names) maxLen = Math.Max(maxLen, s.Length);
			for(int i=0;i<names.Length;i++) 
			{
				if (sb.Length>0) sb.AppendLine();
				sb.Append(" P("+names[i]+")"); sb.Append(' ', maxLen-names[i].Length);
				sb.Append(" = "+this[i]);
			}
			return sb.ToString();
		}

		/// <summary>
		/// Returns a sample from the distribution
		/// </summary>
		/// <returns>The sample value</returns>
		[Stochastic]
		public new TEnum Sample()
		{
			Array values = Enum.GetValues(typeof(TEnum));
			return (TEnum)values.GetValue(base.Sample());
		}

		/// <summary>
		/// Returns a sample from the distribution
		/// </summary>
		/// <param name="result">Not used</param>
		/// <returns>The sample value</returns>
		[Stochastic]
		public TEnum Sample(TEnum result)
		{
			return Sample();
		}
	}

}
