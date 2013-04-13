// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// Represents a list of Beta distributions, optimised for the case where many share 
	/// the same pseudocount values.
	/// </summary>
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public class SparseBetaList : IDistribution<SparseVector>, Sampleable<SparseVector>,
		SettableToProduct<SparseBetaList>,
				SettableTo<SparseBetaList>, SettableToPower<SparseBetaList>, SettableToRatio<SparseBetaList>,
				SettableToWeightedSum<SparseBetaList>, CanGetLogAverageOf<SparseBetaList>, CanGetLogAverageOfPower<SparseBetaList>,
								CanGetAverageLog<SparseBetaList>
	{
		#region Properties
		/// <summary>
		/// The sparse vector holding the true counts of each Beta in the list
		/// </summary>
		[System.Xml.Serialization.XmlElement(Type=typeof(SerializableVector))]
		public ApproximateSparseVector TrueCounts { get; set; }

		/// <summary>
		/// The sparse vector holding the false counts of each Beta in the list
		/// </summary>
		[System.Xml.Serialization.XmlElement(Type=typeof(SerializableVector))]
		public ApproximateSparseVector FalseCounts { get; set; }

		/// <summary>
		/// The number of elements in the list
		/// </summary>
		public int Count { get { return TrueCounts.Count; } }
		#endregion

		#region Constructors and factory methods

		/// <summary>
		/// Parameterless constructor required for serialization 
		/// </summary>
		private SparseBetaList() { }

		/// <summary>
		/// Creates a list of Betas of the specified size, each set to uniform.
		/// </summary>
		/// <param name="size">The size of the list</param>
		public SparseBetaList(int size)
			: this(size, 1.0, 1.0)
		{
		}

		/// <summary>
		/// The default sparsity settings.
		/// </summary>
		public static Sparsity DefaultSparsity = Sparsity.ApproximateWithTolerance(0.01);

		/// <summary>
		/// Creates a list of Betas of the specified size, each set to have the specified true and false counts.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="trueCount">The true count for each Beta distribution</param>
		/// <param name="falseCount">The false count for each Beta distribution</param>
		public SparseBetaList(int size, double trueCount, double falseCount)
		{
			TrueCounts = ApproximateSparseVector.Constant(size, trueCount, DefaultSparsity);
			FalseCounts = ApproximateSparseVector.Constant(size, falseCount, DefaultSparsity);
		}

		/// <summary>
		/// Creates a list of Betas using the supplied vectors of true and false counts.
		/// </summary>
		/// <param name="trueCounts">The vector of true counts</param>
		/// <param name="falseCounts">The vector of false counts</param>
		public SparseBetaList(Vector trueCounts, Vector falseCounts)
		{
			if (trueCounts.Count != falseCounts.Count) {
				throw new ArgumentException("True and false count vectors must be the same size " +
                    trueCounts.Count + "!=" + falseCounts.Count);
			}
			TrueCounts = ApproximateSparseVector.Copy(trueCounts);
			FalseCounts = ApproximateSparseVector.Copy(falseCounts);
		}

		/// <summary>
		/// Creates a point mass SparseBetaList
		/// </summary>
		/// <param name="probsTrue"></param>
		/// <returns></returns>
		public static SparseBetaList PointMass(SparseVector probsTrue)
		{
			var falseCounts = SparseVector.Constant(probsTrue.Count, double.PositiveInfinity);
			SparseBetaList sbl = new SparseBetaList(probsTrue, falseCounts);
			return sbl;
		}

		#endregion

		#region Element set/get

		/// <summary>Gets or sets a Beta element.</summary>
		public Beta this[int index]
		{
			get
			{
				return new Beta(TrueCounts[index], FalseCounts[index]);
			}
			set
			{
				TrueCounts[index] = value.TrueCount;
				FalseCounts[index] = value.FalseCount;
			}
		}

		#endregion

		#region Expected value, expected logs
		/// <summary>
		/// The expected value E[p].
		/// </summary>
		/// <returns>Sparse vector of the means of the elements</returns>
		/// <remarks>Elements of the result must be between 0 and 1.</remarks>
		public SparseVector GetMean()
		{
			var mean = SparseVector.Zero(Count);
			mean.SetToFunction(TrueCounts, FalseCounts, ComputeMean);
			return mean;
		}

		static double ComputeMean(double trueCount, double falseCount)
		{
			if (double.IsPositiveInfinity(falseCount)) return trueCount;
			if ((trueCount == 0.0) && (falseCount == 0.0)) return 0.5;
			return trueCount / (trueCount + falseCount);
		}

		/// <summary>
		/// The expected value E[log(p)].
		/// </summary>
		public SparseVector GetMeanLog()
		{
			var meanLog = SparseVector.Zero(Count);
			meanLog.SetToFunction(TrueCounts, FalseCounts, ComputeMeanLog);
			return meanLog;
		}

		static internal double ComputeMeanLog(double trueCount, double falseCount)
		{
			if (double.IsPositiveInfinity(falseCount)) return Math.Log(trueCount);
			if ((trueCount == 0.0) && (falseCount == 0.0)) throw new ImproperDistributionException(new Beta(trueCount, falseCount));
			return MMath.Digamma(trueCount) - MMath.Digamma(trueCount + falseCount);
		}

		SparseVector meanLogOneMinus;
		/// <summary>
		/// The expected value E[log(1-p)].
		/// </summary>
		public SparseVector GetMeanLogOneMinus()
		{
			if (meanLogOneMinus==null) meanLogOneMinus = SparseVector.Zero(Count);
			meanLogOneMinus.SetToFunction(TrueCounts, FalseCounts, ComputeMeanLogOneMinus);
			return meanLogOneMinus;
		}

		static internal double ComputeMeanLogOneMinus(double trueCount, double falseCount)
		{
			if (double.IsPositiveInfinity(falseCount)) return Math.Log(1 - trueCount);
			if ((trueCount == 0.0) && (falseCount == 0.0)) throw new ImproperDistributionException(new Beta(trueCount, falseCount));
			return MMath.Digamma(falseCount) - MMath.Digamma(trueCount + falseCount);
		}


		#endregion

		#region IDistribution<SparseVector> Members
		/// <summary>
		/// Gets/sets the distribution as a point value list of doubles
		/// </summary>
		[NonSerializedProperty, System.Xml.Serialization.XmlIgnore]
		public SparseVector Point
		{
			get
			{
				return TrueCounts;
			}
			set
			{
				foreach (double d in value) {
					if (d < 0.0 || d > 1.0) throw new ArgumentException("Supplied point value has at least one element outside [0,1]", "value");
				}
				FalseCounts.SetAllElementsTo(Double.PositiveInfinity);
				TrueCounts.SetTo(value);
			}
		}

		/// <summary>
		/// Whether all the Beta elements are point masses
		/// </summary>
		public bool IsPointMass
		{
			get { return FalseCounts.All(x => Double.IsPositiveInfinity(x)); }
		}

		/// <summary>
		/// The maximum 'difference' between corresponding Beta elements.
		/// </summary>
		public double MaxDiff(object thatd)
		{
			SparseBetaList that = (SparseBetaList)thatd;
			// differences in true counts
			SparseVector trueDiff =SparseVector.Zero(Count);
			trueDiff.SetToFunction(TrueCounts, that.TrueCounts, (a, b) => MMath.AbsDiff(a, b));
			// differences in false counts
			SparseVector falseDiff = SparseVector.Zero(Count);
			falseDiff.SetToFunction(FalseCounts, that.FalseCounts, (a, b) => MMath.AbsDiff(a, b));
			// max overall difference
			return Math.Max(trueDiff.Max(), falseDiff.Max());
		}

		/// <summary>
		/// Sets all Beta elements to uniform.
		/// </summary>
		public void SetToUniform()
		{
			TrueCounts.SetAllElementsTo(1);
			FalseCounts.SetAllElementsTo(1);
		}

		/// <summary>
		/// Tests if all Bernoulli elements are uniform.
		/// </summary>
		/// <returns></returns>
		public bool IsUniform()
		{
			return TrueCounts.EqualsAll(1) && FalseCounts.EqualsAll(1);
		}

		/// <summary>
		/// Evaluates the logarithm of the density function
		/// </summary>
		/// <param name="value"></param>
		/// <returns></returns>
		public double GetLogProb(SparseVector value)
		{
			if (value.Count != Count) throw new ArgumentException("Point value had incorrect size of " + value.Count + " != " + Count);
			// TODO: replace by an efficient sparse implementation
			double logprob = 0;
			for (int i = 0; i < Count; i++) {
				logprob += this[i].GetLogProb(value[i]);
			}
			return logprob;
		}

		#endregion

		#region Product, power

		/// <summary>
		/// Sets to the product of two sparse lists of Betas.
		/// </summary>
		/// <param name="a">The first sparse list of Betas</param>
		/// <param name="b">The second sparse list of Betas</param>
		/// <remarks>
		/// The result may not be proper, i.e. its parameters may be negative.
		/// For example, if you multiply Beta(0.1,0.1) by itself you get Beta(-0.8, -0.8).
		/// No error is thrown in this case.
		/// </remarks>
		public void SetToProduct(SparseBetaList a, SparseBetaList b)
		{
			// TODO: handle point masses
			TrueCounts.SetToFunction(a.TrueCounts, b.TrueCounts, (x, y) => x + y - 1);
			FalseCounts.SetToFunction(a.FalseCounts, b.FalseCounts, (x, y) => x + y - 1);
		}

		/// <summary>
		/// Sets the parameters to represent the list of sparse Betas raised to some power.
		/// </summary>
		/// <param name="betaList">The list of Betas</param>
		/// <param name="exponent">The exponent</param>
		public void SetToPower(SparseBetaList betaList, double exponent)
		{
			// TODO: handle point masses
			TrueCounts.SetToFunction(betaList.TrueCounts, x => exponent * (x - 1) + 1);
			FalseCounts.SetToFunction(betaList.FalseCounts, x => exponent * (x - 1) + 1);
		}
		#endregion

		#region SetTo, Clone

		/// <summary>
		/// Sets this Beta list to a copy of another Beta list.
		/// </summary>
		/// <param name="value">The list to copy</param>
		public void SetTo(SparseBetaList value)
		{
			TrueCounts.SetTo(value.TrueCounts);
			FalseCounts.SetTo(value.FalseCounts);
		}

		/// <summary>
		/// Clones this sparse Beta list.
		/// </summary>
		public object Clone()
		{
			return new SparseBetaList((SparseVector)TrueCounts.Clone(), (SparseVector)FalseCounts.Clone());
		}

		#endregion

		#region Sampleable<SparseVector> Members

		/// <summary>
		/// Samples a vector from this distribution.
		/// </summary>
		/// <returns></returns>
		/// <remarks>The samples themselves will not be sparse</remarks>
		public SparseVector Sample()
		{
			// TODO: change to return a dense vector
			// Should not be used for sampling since the samples themselves will not be sparse 
			SparseVector sv = SparseVector.Zero(Count);
			return Sample(sv);
		}

		/// <summary>
		/// Samples a vector from this distribution.
		/// </summary>
		/// <param name="sv">Where to put the result</param>
		/// <returns></returns>
		/// <remarks>The samples themselves will not be sparse</remarks>
		public SparseVector Sample(SparseVector sv)
		{
			// Should not be used for sampling since the samples themselves will not be sparse 
			for (int i=0; i<Count; i++) sv[i] = this[i].Sample();
			return sv;
		}

		#endregion

		#region SettableToWeightedSum<SparseBetaList> Members

		/// <summary>
		/// Sets this sparse Bernoulli list to the weighted sum of two others (which must be of the same size)
		/// </summary>
		/// <param name="weight1"></param>
		/// <param name="value1"></param>
		/// <param name="weight2"></param>
		/// <param name="value2"></param>
		public void SetToSum(double weight1, SparseBetaList value1, double weight2, SparseBetaList value2)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region SettableToRatio<SparseBetaList,SparseBetaList> Members

		/// <summary>
		/// Sets this sparse Bernoulli list to the ratio of two others (which must be of the same size)
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(SparseBetaList numerator, SparseBetaList denominator)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region CanGetLogAverageOf<SparseBetaList> Members

		/// <summary>
		/// Gets the integral of the product of two sparse Beta list distributions.
		/// </summary>
		/// <remarks>Not yet implemented</remarks>
		public double GetLogAverageOf(SparseBetaList that)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region CanGetAverageLog<SparseBetaList> Members

		/// <summary>
		/// Gets the expected logarithm of that distribution under this distribution.
		/// </summary>
		/// <param name="that">The distribution to take the logarithm of.</param>
		/// <returns><c>sum_x this.Evaluate(x)*Math.Log(that.Evaluate(x))</c></returns>
		/// <remarks>Not yet implemented.</remarks>
		public double GetAverageLog(SparseBetaList that)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region Conversions: ToString, ToArray
		/// <summary>
		/// ToString override
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			string s = "BetaList(TrueCounts:" + TrueCounts + ", FalseCounts:" + FalseCounts + ")";
			return s;
		}

		/// <summary>
		/// Converts this sparse array to a dense .NET array.
		/// </summary>
		/// <returns></returns>
		public Beta[] ToArray()
		{
			Beta[] b = new Beta[Count];
			for (int i=0; i<b.Length; i++) b[i]=this[i];
			return b;
		}

		#endregion

		public double GetLogAverageOfPower(SparseBetaList that, double power)
		{
			throw new NotImplementedException();
		}
	}
}
