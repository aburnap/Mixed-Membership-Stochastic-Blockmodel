// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// Represents a list of Bernoulli distributions, optimised for the case where many share 
	/// the same probability of being true.
	/// This class can be used as a distribution over a fixed-sized list of booleans or sparsely 
	/// as a distribution over a variable-sized list of integers, which are the indices of elements 
	/// in the boolean list with value 'true'.
	/// </summary>
	/// <remarks>
	/// One use for this class is as a distribution over subsets of a set of objects.
	/// </remarks>
	/// TODO: operators
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public class SparseBernoulliList : SparseBernoulliListBase, IDistribution<IList<bool>>, Sampleable<IList<bool>>,
				 SettableTo<SparseBernoulliList>, SettableToProduct<SparseBernoulliList>,
		SettableToPower<SparseBernoulliList>, SettableToRatio<SparseBernoulliList>,
		SettableToWeightedSum<SparseBernoulliList>, CanGetLogAverageOf<SparseBernoulliList>, CanGetLogAverageOfPower<SparseBernoulliList>,
				CanGetAverageLog<SparseBernoulliList>
	{
		/// <summary>
		/// Default sparsity specification for this distribution
		/// </summary>
		public static Sparsity DefaultSparsity = Sparsity.ApproximateWithTolerance(0.01);

		#region Constructors
		/// <summary>
		/// Parameterless constructor required for serialization 
		/// </summary>
		private SparseBernoulliList() { }

		/// <summary>
		/// Creates a list of Bernoullis of the specified size, each set to uniform.
		/// </summary>
		/// <param name="size">The size of the list</param>
		public SparseBernoulliList(int size)
		{
			LogOddsVector = ApproximateSparseVector.Zero(size, DefaultSparsity);
		}

		/// <summary>
		/// Creates a list of Bernoullis of the specified size, each set to have the specified probability of being true.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="probTrue">The probability of each Bernoulli being true</param>
		public SparseBernoulliList(int size, double probTrue)
		{
			LogOddsVector = ApproximateSparseVector.Constant(size, MMath.Logit(probTrue), DefaultSparsity);
		}

		/// <summary>
		/// Creates a list of Bernoullis using the supplied vector of log odds (which will *not* be cloned).
		/// </summary>
		/// <param name="logOddsVector">The log odds vector to use internally</param>
		[Construction("LogOddsVector")]
		public SparseBernoulliList(Vector logOddsVector)
			: this(logOddsVector.Count)
		{
			LogOddsVector.SetTo(logOddsVector);
		}
		#endregion

		# region Conversions e.g. ToString
		/// <summary>
		/// Returns a human readable string form of this distribution.
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			return "SparseBernoulliList(ProbsTrue:"+GetProbTrueVector()+")";
		}
		#endregion

		/// <summary>
		/// Converts this sparse array to an ordinary dense .NET array.
		/// </summary>
		/// <returns>A .NET array of Bernoulli values</returns>
		public Bernoulli[] ToArray()
		{
			Bernoulli[] array = new Bernoulli[Count];
			for (int i=0; i<array.Length; i++) array[i] = this[i];
			return array;
		}

		#region IDistribution<IList<bool>> Members

		/// <summary>
		/// Clones this sparse Bernoulli list.
		/// </summary>
		/// <returns></returns>
		public object Clone()
		{
			SparseBernoulliList clone = new SparseBernoulliList((SparseVector)LogOddsVector.Clone());
			return clone;
		}

		/// <summary>
		/// Gets/sets the distribution as a point value boolean array
		/// </summary>
		[NonSerializedProperty, System.Xml.Serialization.XmlIgnore]
		public IList<bool> Point
		{
			get
			{
				IList<bool> b = new bool[Count];
				for (int i = 0; i < Count; i++) b[i] = Double.IsPositiveInfinity(LogOddsVector[i]);
				return b;
			}
			set
			{
				if (value.Count != Count) throw new ArgumentException("Point value had incorrect size of "+value.Count+" != "+Count);
				// Set all to false
				LogOddsVector.SetAllElementsTo(Double.NegativeInfinity);
				for (int i = 0; i < Count; i++) if (value[i]) LogOddsVector[i] = Double.PositiveInfinity;
			}
		}

		/// <summary>
		/// Evaluates the logarithm of the density function
		/// </summary>
		/// <param name="value">a list of boolean values</param>
		/// <returns>Log of the probability density for the given event</returns>
		public double GetLogProb(IList<bool> value)
		{
			if (value.Count != Count) throw new ArgumentException("Point value had incorrect size of " + value.Count + " != " + Count);
			SparseVector tv = SparseVector.Zero(LogOddsVector.Count);
			tv.SetToFunction(LogOddsVector, x => MMath.LogisticLn(-x));
			double logprob = tv.Sum();
			for (int i = 0; i < Count; i++) {
				if (!value[i]) continue;
				double logOdds = LogOddsVector[i];
				logprob -= MMath.LogisticLn(-logOdds);
				logprob += MMath.LogisticLn(logOdds);
			}
			return logprob;
		}

		#endregion

		#region Sampleable<IList<bool>> Members

		/// <summary>
		/// Samples a list of booleans from this distribution.
		/// </summary>
		/// <returns></returns>
		public IList<bool> Sample()
		{
			return Sample(new bool[Count]);
		}

		/// <summary>
		/// Samples a list of booleans from this distribution, using the supplied list for storage.
		/// </summary>
		/// <returns></returns>
		public IList<bool> Sample(IList<bool> result)
		{
			int i=0;
			foreach (double d in LogOddsVector) result[i++]=Bernoulli.Sample(MMath.Logistic(d));
			return result;
		}

		/// <summary>
		/// Samples from a list of Bernoulli distributions with the specified vector of P(true) values
		/// </summary>
		/// <param name="probsTrue">the vector of P(true) values</param>
		/// <returns>The sample</returns>
		[Stochastic]
		public static IList<bool> Sample(SparseVector probsTrue)
		{
			var result = new bool[probsTrue.Count];
			int i = 0;
			foreach (double d in probsTrue) result[i++]=Bernoulli.Sample(d);
			return result;
		}

		#endregion

		#region SettableToProduct<SparseBernoulliList,SparseBernoulliList> Members

		/// <summary>
		/// Sets this list to the elementwise product of the two supplied lists.
		/// </summary>
		/// <param name="a">The first list</param>
		/// <param name="b">The second list</param>
		/// <returns>The resulting list of Bernoulli distributions</returns>
		public void SetToProduct(SparseBernoulliList a, SparseBernoulliList b)
		{
			LogOddsVector.SetToSum(a.LogOddsVector, b.LogOddsVector);
		}

		#endregion

		#region SettableTo<SparseBernoulliList> Members

		/// <summary>
		/// Sets this lists to copy the state of another list (which must be the same size).
		/// </summary>
		/// <param name="otherList">The list to copy</param>
		public void SetTo(SparseBernoulliList otherList)
		{
			LogOddsVector.SetTo(otherList.LogOddsVector);
		}

		#endregion

		#region SettableToWeightedSum<SparseBernoulliList> Members

		/// <summary>
		/// Sets this sparse Bernoulli list to the weighted sum of two others (which must be of the same size)
		/// </summary>
		/// <param name="weight1"></param>
		/// <param name="value1"></param>
		/// <param name="weight2"></param>
		/// <param name="value2"></param>
		public void SetToSum(double weight1, SparseBernoulliList value1, double weight2, SparseBernoulliList value2)
		{
			if (weight1 + weight2 == 0) SetToUniform();
			else {
				LogOddsVector.SetToFunction(value1.GetProbTrueVector(), value2.GetProbTrueVector(), (x, y) => MMath.Logit((weight1 * x + weight2 *y)/ (weight1+weight2)));
			}
		}

		#endregion

		#region SettableToRatio<SparseBernoulliList,SparseBernoulliList> Members

		/// <summary>
		/// Sets this sparse Bernoulli list to the ratio of two others (which must be of the same size)
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(SparseBernoulliList numerator, SparseBernoulliList denominator)
		{
			// TODO: point masses
			LogOddsVector.SetToDifference(numerator.LogOddsVector, denominator.LogOddsVector);
		}

		#endregion

		#region SettableToPower<SparseBernoulliList> Members

		/// <summary>
		/// Sets this sparse Bernoulli list to the power of another (which must be of the same size)
		/// </summary>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		/// <remarks>Not yet implemented</remarks>
		public void SetToPower(SparseBernoulliList value, double exponent)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region CanGetLogAverageOf<SparseBernoulliList> Members

		/// <summary>
		/// Gets the integral of the product of two sparse Bernoulli list distributions.
		/// </summary>
		/// <remarks>Not yet implemented</remarks>
		public double GetLogAverageOf(SparseBernoulliList that)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region CanGetAverageLog<SparseBernoulliList> Members

		/// <summary>
		/// Gets the expected logarithm of that distribution under this distribution.
		/// </summary>
		/// <param name="that">The distribution to take the logarithm of.</param>
		/// <returns><c>sum_x this.Evaluate(x)*Math.Log(that.Evaluate(x))</c></returns>
		/// <remarks>This is also known as the cross entropy.</remarks>
		public double GetAverageLog(SparseBernoulliList that)
		{
			var res = SparseVector.Zero(LogOddsVector.Count);
			var res2 = SparseVector.Zero(LogOddsVector.Count);
			res.SetToFunction(LogOddsVector, that.GetLogProbTrueVector(), (logpr, logProbTrue) => double.IsNegativeInfinity(logpr) ?0 :MMath.Logistic(logpr)*logProbTrue);
			res2.SetToFunction(LogOddsVector, that.GetLogProbFalseVector(), (logpr, logProbFalse) => double.IsPositiveInfinity(logpr) ?0 :MMath.Logistic(-logpr)*logProbFalse);
			return res.Sum() + res2.Sum();
		}

		#endregion

		public double GetLogAverageOfPower(SparseBernoulliList that, double power)
		{
			throw new NotImplementedException();
		}
	}

	/// <summary>
	/// Represents a list of Bernoulli distributions considered as a distribution over a variable-sized list of
	/// integers, which are the indices of elements in the boolean list with value 'true'
	/// </summary>
	[Quality(QualityBand.Experimental)]
	public class BernoulliIntegerSubset : SparseBernoulliListBase, IDistribution<IList<int>>, Sampleable<IList<int>>,
				SettableToProduct<BernoulliIntegerSubset>,
		SettableTo<BernoulliIntegerSubset>, SettableToPower<BernoulliIntegerSubset>, SettableToRatio<BernoulliIntegerSubset>,
		SettableToWeightedSum<BernoulliIntegerSubset>, CanGetLogAverageOf<BernoulliIntegerSubset>,
				CanGetAverageLog<BernoulliIntegerSubset>
	{

		#region Constructors
		/// <summary>
		/// Creates a list of Bernoullis of the specified size, each set to uniform.
		/// </summary>
		/// <param name="size">The size of the list</param>
		public BernoulliIntegerSubset(int size)
		{
			LogOddsVector = ApproximateSparseVector.Zero(size, SparseBernoulliList.DefaultSparsity);
		}

		/// <summary>
		/// Creates a list of Bernoullis of the specified size, each set to have the specified probability of being true.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="probTrue">The probability of each Bernoulli being true</param>
		public BernoulliIntegerSubset(int size, double probTrue)
		{
			LogOddsVector = ApproximateSparseVector.Constant(size, MMath.Logit(probTrue), SparseBernoulliList.DefaultSparsity);
		}

		/// <summary>
		/// Creates a list of Bernoullis using the supplied vector of log odds (which will *not* be cloned).
		/// </summary>
		/// <param name="logOddsVector">The log odds vector to use internally</param>
		[Construction("LogOddsVector")]
		public BernoulliIntegerSubset(Vector logOddsVector)
			: this(logOddsVector.Count)
		{
			LogOddsVector.SetTo(logOddsVector);
		}
		#endregion

		#region IDistribution<IList<int>> Members
		/// <summary>
		/// Gets/sets the distribution as a point value int array. This is a sparse method - the ints specify
		/// the indices where elements are true - all other elements are false.
		/// </summary>
		[NonSerializedProperty, System.Xml.Serialization.XmlIgnore]
		public IList<int> Point
		{
			get
			{
				var trueIndices = LogOddsVector.IndexOfAll(Double.IsPositiveInfinity);
				return new List<int>(trueIndices);
			}
			set
			{
				// Set all to false
				LogOddsVector.SetAllElementsTo(Double.NegativeInfinity);
				// Set specified elements to true 
				for (int i = 0; i < value.Count; i++) LogOddsVector[value[i]] = Double.PositiveInfinity;
			}
		}

		/// <summary>
		/// Gets the log probability of the given value under this distribution
		/// </summary>
		/// <param name="value"></param>
		/// <returns></returns>
		public double GetLogProb(IList<int> value)
		{
			SparseVector tv = SparseVector.Zero(LogOddsVector.Count);
			tv.SetToFunction(LogOddsVector, x => MMath.LogisticLn(-x));
			double logprob = tv.Sum();
			foreach (int index in value) {
				double logOdds = LogOddsVector[index];
				logprob -= MMath.LogisticLn(-logOdds);
				logprob += MMath.LogisticLn(logOdds);
			}
			return logprob;
		}

		/// <summary>
		/// Clones this object.
		/// </summary>
		/// <returns></returns>
		public object Clone()
		{
			BernoulliIntegerSubset clone = new BernoulliIntegerSubset((SparseVector)LogOddsVector.Clone());
			return clone;
		}
		#endregion

		#region Sampleable<IList<int>> Members

		/// <summary>
		/// Samples a list of ints from this distribution.
		/// </summary>
		/// <returns></returns>
		public IList<int> Sample()
		{
			return Sample(new List<int>());
		}

		/// <summary>
		/// Samples a list of ints from this distribution
		/// </summary>
		/// <param name="result">Where to put the resulting sample</param>
		/// <returns></returns>
		public IList<int> Sample(IList<int> result)
		{
			// TODO: efficient sparse implementation
			result.Clear();
			int i=0;
			foreach (double d in LogOddsVector) {
				if (Bernoulli.Sample(MMath.Logistic(d))) result.Add(i);
				i++;
			}
			return result;
		}

		/// <summary>
		/// Samples from a list of Bernoulli distributions with the specified vector of P(true) values
		/// </summary>
		/// <param name="probsTrue">the vector of P(true) values</param>
		/// <returns>The sample</returns>
		[Stochastic]
		public static IList<int> Sample(SparseVector probsTrue)
		{
			var result = new bool[probsTrue.Count];
			var list = new List<int>();
			int i=0;
			foreach (double d in probsTrue) {
				if (Bernoulli.Sample(d)) list.Add(i);
				i++;
			}
			return list;
		}

		#endregion

		#region SettableToProduct<BernoulliIntegerSubset,BernoulliIntegerSubset> Members

		/// <summary>
		/// Sets this BernoulliIntegerSubset distribution to a product of two other such distributions
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public void SetToProduct(BernoulliIntegerSubset a, BernoulliIntegerSubset b)
		{
			LogOddsVector.SetToSum(a.LogOddsVector, b.LogOddsVector);
		}

		#endregion

		#region SettableTo<BernoulliIntegerSubset> Members

		/// <summary>
		/// Sets this BernoulliIntegerSubset distribution to another such distribution
		/// </summary>
		/// <param name="value"></param>
		public void SetTo(BernoulliIntegerSubset value)
		{
			LogOddsVector.SetTo(value.LogOddsVector);
		}

		#endregion

		#region SettableToPower<BernoulliIntegerSubset> Members

		/// <summary>
		/// Sets this BernoulliIntegerSubset distribution to the power another such distribution
		/// </summary>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		public void SetToPower(BernoulliIntegerSubset value, double exponent)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region SettableToRatio<BernoulliIntegerSubset,BernoulliIntegerSubset> Members

		/// <summary>
		/// Sets this BernoulliIntegerSubset distribution to the ratio of two other such distributions
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(BernoulliIntegerSubset numerator, BernoulliIntegerSubset denominator)
		{
			LogOddsVector.SetToDifference(numerator.LogOddsVector, denominator.LogOddsVector);
		}

		#endregion

		#region SettableToWeightedSum<BernoulliIntegerSubset> Members

		/// <summary>
		/// Sets this BernoulliIntegerSubset distribution to the weighted sum of two other such distributions
		/// </summary>
		/// <param name="weight1"></param>
		/// <param name="value1"></param>
		/// <param name="weight2"></param>
		/// <param name="value2"></param>
		/// <remarks>Not yet implemented</remarks>
		public void SetToSum(double weight1, BernoulliIntegerSubset value1, double weight2, BernoulliIntegerSubset value2)
		{
			if (weight1 + weight2 == 0) SetToUniform();
			else {
				LogOddsVector.SetToFunction(value1.GetProbTrueVector(), value2.GetProbTrueVector(), (x, y) => MMath.Logit((weight1 * x + weight2 *y)/ (weight1+weight2)));
			}
		}

		#endregion

		#region CanGetLogAverageOf<BernoulliIntegerSubset> Members

		/// <summary>
		/// Gets the log of the integral of the product of this BernoulliIntegerSubset distribution and another such distribution
		/// </summary>
		/// <param name="that"></param>
		/// <returns></returns>
		/// <remarks>Not yet implemented</remarks>
		public double GetLogAverageOf(BernoulliIntegerSubset that)
		{
			throw new NotImplementedException();
		}

		#endregion

		#region CanGetAverageLog<BernoulliIntegerSubset> Members

		/// <summary>
		/// The expected logarithm of that distribution under this distribution.
		/// </summary>
		/// <param name="that">The distribution to take the logarithm of.</param>
		/// <returns><c>sum_x this.Evaluate(x)*Math.Log(that.Evaluate(x))</c></returns>
		/// <remarks>This is also known as the cross entropy.</remarks>
		public double GetAverageLog(BernoulliIntegerSubset that)
		{
			var res = SparseVector.Zero(LogOddsVector.Count);
			var res2 = SparseVector.Zero(LogOddsVector.Count);
			res.SetToFunction(LogOddsVector, that.GetLogProbTrueVector(), (logpr, logProbTrue) => double.IsNegativeInfinity(logpr) ?0 :MMath.Logistic(logpr)*logProbTrue);
			res2.SetToFunction(LogOddsVector, that.GetLogProbFalseVector(), (logpr, logProbFalse) => double.IsPositiveInfinity(logpr) ?0 :MMath.Logistic(-logpr)*logProbFalse);
			return res.Sum() + res2.Sum();
		}

		#endregion
	}

	/// <summary>
	/// Base class for <see cref="BernoulliIntegerSubset"/> and <see cref="SparseBernoulliList"/>.
	/// </summary>
	[Serializable]
	public class SparseBernoulliListBase : SettableToUniform
	{
		#region Properties
		/// <summary>
		/// The sparse vector holding the log odds of each Bernoulli in the list
		/// </summary>
		[System.Xml.Serialization.XmlElement(Type=typeof(SerializableVector))]
		public ApproximateSparseVector LogOddsVector { get; set; }

		/// <summary>
		/// The number of elements in the list
		/// </summary>
		public int Count { get { return LogOddsVector.Count; } }
		#endregion

		#region Element set/get

		/// <summary>Gets or sets an element.</summary>
		public Bernoulli this[int index]
		{
			get
			{
				return Bernoulli.FromLogOdds(LogOddsVector[index]);
			}
			set
			{
				LogOddsVector[index] = value.LogOdds;
			}
		}

		/// <summary>
		/// Gets the probability of the binary variable being true
		/// </summary>
		/// <returns>p(x=true)</returns>
		public double GetProbTrue(int index)
		{
			return MMath.Logistic(LogOddsVector[index]);
		}
		/// <summary>
		/// Sets the probability of the binary variable being true
		/// </summary>
		public void SetProbTrue(int index, double probTrue)
		{
			if (probTrue < 0 || probTrue > 1) throw new ArgumentException(String.Format("probTrue = {0} is not in [0,1]", probTrue), "probTrue");
			LogOddsVector[index] = MMath.Logit(probTrue);
		}

		/// <summary>
		/// Gets the probability of the binary variable being false
		/// </summary>
		/// <returns>p(x=false)</returns>
		public double GetProbFalse(int index)
		{
			return MMath.Logistic(-LogOddsVector[index]);
		}
		/// <summary>
		/// Sets the probability of the binary variable being false
		/// </summary>
		public void SetProbFalse(int index, double probFalse)
		{
			if (probFalse < 0 || probFalse > 1) throw new ArgumentException(String.Format("probFalse = {0} is not in [0,1]", probFalse), "probTrue");
			LogOddsVector[index] = -MMath.Logit(probFalse);
		}
		#endregion

		/// <summary>
		/// Sets all Bernoulli elements to uniform.
		/// </summary>
		public void SetToUniform()
		{
			LogOddsVector.SetAllElementsTo(0.0);
		}

		/// <summary>
		/// Tests if all Bernoulli elements are uniform.
		/// </summary>
		/// <returns></returns>
		public bool IsUniform()
		{
			return LogOddsVector.EqualsAll(0);
		}


		/// <summary>
		/// Whether all the Bernoulli elements are point masses
		/// </summary>
		public bool IsPointMass
		{
			get
			{
				return LogOddsVector.All(Double.IsInfinity);
			}
		}

		/// <summary>
		/// The maximum 'difference' between this instance and that instance.
		/// This returns the maximum absolute difference between the Log-odds of any element
		/// </summary>
		/// <param name="that">The other distribution</param>
		/// <returns>The resulting maximum difference</returns>
		/// <remarks><c>a.MaxDiff(b) == b.MaxDiff(a)</c></remarks>
		public double MaxDiff(object that)
		{
			SparseBernoulliListBase sbl = (SparseBernoulliListBase)that;
			// differences in log odds
			SparseVector absDiff = SparseVector.Zero(LogOddsVector.Count);
			absDiff.SetToFunction(LogOddsVector, sbl.LogOddsVector, (a, b) => MMath.AbsDiff(a, b));
			return absDiff.Max();
		}

		#region Vectors of prob true, false etc.
		/// <summary>
		/// Gets a vector of P(true) values.
		/// </summary>
		/// <returns></returns>
		public SparseVector GetProbTrueVector()
		{
			SparseVector sv = SparseVector.Zero(Count);
			sv.SetToFunction(LogOddsVector, MMath.Logistic);
			return sv;
		}

		/// <summary>
		/// Gets a vector of P(true) values.
		/// </summary>
		/// <returns></returns>
		public void SetProbTrueVector(SparseVector probTrueVector)
		{
			LogOddsVector.SetToFunction(probTrueVector, MMath.Logit);
		}

		/// <summary>
		/// Gets a vector of log P(true) values.
		/// </summary>
		/// <returns></returns>
		public SparseVector GetLogProbTrueVector()
		{
			SparseVector sv = SparseVector.Zero(Count);
			sv.SetToFunction(LogOddsVector, MMath.LogisticLn);
			return sv;
		}

		/// <summary>
		/// Gets a vector of log P(false) values.
		/// </summary>
		/// <returns></returns>
		public SparseVector GetLogProbFalseVector()
		{
			SparseVector sv = SparseVector.Zero(Count);
			sv.SetToFunction(LogOddsVector, x => MMath.LogisticLn(-x));
			return sv;
		}

		#endregion
	}
}
