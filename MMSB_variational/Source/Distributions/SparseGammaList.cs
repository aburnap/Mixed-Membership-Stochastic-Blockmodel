using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// Represents a sparse list of Gamma distributions, optimised for the case where many share 
	/// the same parameterisation.
	/// </summary>
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public class SparseGammaList : ApproximateSparseList<Gamma>,
		IDistribution<SparseVector>, Sampleable<SparseVector>,
		SettableTo<SparseGammaList>, SettableToProduct<SparseGammaList>,
		SettableToPower<SparseGammaList>, SettableToRatio<SparseGammaList>,
		SettableToWeightedSum<SparseGammaList>, CanGetLogAverageOf<SparseGammaList>, CanGetLogAverageOfPower<SparseGammaList>,
		CanGetAverageLog<SparseGammaList>
	{
		/// <summary>
		/// Default tolerance for sparsity approximation
		/// </summary>
		public new static double DefaultTolerance = 0.01;
		/// <summary>
		/// Default constructor
		/// </summary>
		public SparseGammaList() : this(DefaultTolerance) { }

		/// <summary>
		/// The dimension of the SparseGammaList domain
		/// </summary>
		public int Dimension
		{
			get { return Count; }
		}
	
		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size">The size of the list</param>
		protected SparseGammaList(int size) : this(size, DefaultTolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma list of the given size, and assigns all
		/// elements to the specified common value
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		protected SparseGammaList(int size, Gamma commonValue) : this(size, commonValue, DefaultTolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		protected SparseGammaList(int size, Gamma commonValue, List<SparseElement<Gamma>> sortedSparseValues) :
			this(size, commonValue, sortedSparseValues, DefaultTolerance) { }

		/// <summary>
		/// Copy constructor
		/// </summary>
		/// <param name="that"></param>
		protected SparseGammaList(SparseGammaList that) : base(that) { }

		/// <summary>
		/// Constructs a sparse Gamma list with the given tolerance
		/// </summary>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public SparseGammaList(double tolerance) : base(tolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		protected SparseGammaList(int size, double tolerance) : base(size, Gamma.Uniform(), tolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma list of the given size, and assigns all
		/// elements to the specified common value
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		protected SparseGammaList(int size, Gamma commonValue, double tolerance) : base(size, commonValue, tolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		protected SparseGammaList(int size, Gamma commonValue, List<SparseElement<Gamma>> sortedSparseValues, double tolerance) :
			base(size, commonValue, sortedSparseValues, tolerance) { }

		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size"></param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		/// <returns></returns>
		public new static SparseGammaList FromSize(int size, double tolerance)
		{
			return new SparseGammaList(size, tolerance);
		}

		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to the specified Gamma
		/// </summary>
		/// <param name="size">Size</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		/// <returns></returns>
		public static SparseGammaList FromGamma(int size, Gamma commonValue, double tolerance)
		{
			return new SparseGammaList(size, commonValue, tolerance);
		}

		/// <summary>
		/// Returns a sparse Gamma list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		[Construction("Count", "CommonValue", "SparseValues", "Tolerance")]
		public new static SparseGammaList FromSparseValues(int size, Gamma commonValue,
			List<SparseElement<Gamma>> sortedSparseValues, double tolerance)
		{
			return new SparseGammaList(size, commonValue, sortedSparseValues, tolerance);
		}

		/// <summary>
		///Creates a sparse Gamma list of a given size, with each
		/// element having a given mean and variance.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="variance">The desired variance.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromMeanAndVariance(
			int size, double mean, double variance, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromMeanAndVariance(mean, variance), tolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given mean and mean log.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">Desired expected value.</param>
		/// <param name="meanLog">Desired expected logarithm.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromMeanAndMeanLog(
			int size, double mean, double meanLog, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromMeanAndMeanLog(mean, meanLog), tolerance);
		}


		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given log mean and mean log.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="logMean">Log of desired expected value.</param>
		/// <param name="meanLog">Desired expected logarithm.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromLogMeanAndMeanLog(
			int size, double logMean, double meanLog, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromLogMeanAndMeanLog(logMean, meanLog), tolerance);
		}
	
		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given shape and rate.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shape">shape</param>
		/// <param name="rate">rate = 1/scale</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromShapeAndRate(
			int size, double shape, double rate, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromShapeAndRate(shape, rate), tolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size with given shape and rate Vectors.
		/// </summary>
		/// <param name="shape">shape</param>
		/// <param name="rate">rate = 1/scale</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromShapeAndRate(
			Vector shape, Vector rate, double tolerance)
		{
			var result = new SparseGammaList(shape.Count, tolerance);
			result.SetToFunction(shape, rate, (s, r) => Gamma.FromShapeAndRate(s, r));
			return result;
		}
	
		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given shape and scale.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shape">shape</param>
		/// <param name="scale">scale</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromShapeAndScale(
			int size, double shape, double scale, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromShapeAndScale(shape, scale), tolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size with given shape and scale Vectors.
		/// </summary>
		/// <param name="shape">shape</param>
		/// <param name="scale">scale</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromShapeAndScale(
			SparseVector shape, SparseVector scale, double tolerance)
		{
			var result = new SparseGammaList(shape.Count, tolerance);
			result.SetToFunction(shape, scale, (s, sc) => Gamma.FromShapeAndScale(s, sc));
			return result;
		}
	
		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having the specified natural parameters
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shapeMinus1">shape - 1</param>
		/// <param name="rate">rate = 1/scale</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGammaList FromNatural(int size, double shapeMinus1, double rate, double tolerance)
		{
			return new SparseGammaList(size, Gamma.FromNatural(shapeMinus1, rate), tolerance);
		}

		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size"></param>
		/// <returns></returns>
		public new static SparseGammaList FromSize(int size)
		{
			return FromSize(size, DefaultTolerance);
		}

		/// <summary>
		/// Constructs a sparse Gamma with the specified number of elements
		/// all of which are set to the specified Gamma
		/// </summary>
		/// <param name="size">Size</param>
		/// <param name="commonValue">The common value</param>
		/// <returns></returns>
		public static SparseGammaList FromGamma(int size, Gamma commonValue)
		{
			return FromGamma(size, commonValue, DefaultTolerance);
		}

		/// <summary>
		/// Returns a sparse Gamma list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		public new static SparseGammaList FromSparseValues(int size, Gamma commonValue,
			List<SparseElement<Gamma>> sortedSparseValues)
		{
			return FromSparseValues(size, commonValue, sortedSparseValues, DefaultTolerance);
		}

		/// <summary>
		///Creates a sparse Gamma list of a given size, with each
		/// element having a given mean and variance.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="variance">The desired variance.</param>
		public static SparseGammaList FromMeanAndVariance(
			int size, double mean, double variance)
		{
			return FromMeanAndVariance(size, mean, variance, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given mean and mean log.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">Desired expected value.</param>
		/// <param name="meanLog">Desired expected logarithm.</param>
		public static SparseGammaList FromMeanAndMeanLog(
			int size, double mean, double meanLog)
		{
			return FromMeanAndMeanLog(size, mean, meanLog, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given log mean and mean log.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="logMean">Log of desired expected value.</param>
		/// <param name="meanLog">Desired expected logarithm.</param>
		public static SparseGammaList FromLogMeanAndMeanLog(
			int size, double logMean, double meanLog)
		{
			return FromLogMeanAndMeanLog(size, logMean, meanLog, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given shape and rate.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shape">shape</param>
		/// <param name="rate">rate = 1/scale</param>
		public static SparseGammaList FromShapeAndRate(
			int size, double shape, double rate)
		{
			return FromShapeAndRate(size, shape, rate, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size with given shape and rate Vectors.
		/// </summary>
		/// <param name="shape">shape</param>
		/// <param name="rate">rate = 1/scale</param>
		public static SparseGammaList FromShapeAndRate(
			Vector shape, Vector rate)
		{
			return FromShapeAndRate(shape, rate, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having a given shape and scale.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shape">shape</param>
		/// <param name="scale">scale</param>
		public static SparseGammaList FromShapeAndScale(
			int size, double shape, double scale)
		{
			return FromShapeAndScale(size, shape, scale, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gamma list of a given size with given shape and scale Vectors.
		/// </summary>
		/// <param name="shape">shape</param>
		/// <param name="scale">scale</param>
		public static SparseGammaList FromShapeAndScale(
			SparseVector shape, SparseVector scale)
		{
			return FromShapeAndScale(shape, scale, DefaultTolerance);
		}
	
		/// <summary>
		/// Creates a sparse Gamma list of a given size, with each
		/// element having the specified natural parameters
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="shapeMinus1">shape - 1</param>
		/// <param name="rate">rate = 1/scale</param>
		public static SparseGammaList FromNatural(int size, double shapeMinus1, double rate)
		{
			return FromNatural(size, shapeMinus1, rate, DefaultTolerance);
		}

		/// <summary>
		/// Sets/gets the instance as a point mass
		/// </summary>
		[NonSerializedProperty, System.Xml.Serialization.XmlIgnore]
		public SparseVector Point
		{
			get
			{
				return SparseVector.FromSparseValues(
					Count, CommonValue.Point,
					SparseValues.Select(x => new SparseElement<double>(x.Index, x.Value.Point)).ToList());
			}
			set
			{
				CommonValue = Gamma.PointMass(value.CommonValue);
				SparseValues = value.SparseValues.Select(x => new SparseElement<Gamma>(x.Index, Gamma.PointMass(x.Value))).ToList();
			}
		}
		/// <summary>
		/// Asks whether the instance is a point mass
		/// </summary>
		public bool IsPointMass
		{
			get { return All(x => IsPointMass); }
		}
		/// <summary>
		/// Returns the maximum difference between the parameters of this sparse Gamma list
		/// and another
		/// </summary>
		/// <param name="thatd">The other sparse Gamma list</param>
		/// <returns>The maximum difference</returns>
		/// <remarks><c>a.MaxDiff(b) == b.MaxDiff(a)</c></remarks>
		public double MaxDiff(object thatd)
		{
			SparseGammaList that = (SparseGammaList)thatd;
			return Reduce<double, Gamma>(
				double.NegativeInfinity, that, (x, y, z) => Math.Max(x, y.MaxDiff(z)));
		}

		/// <summary>
		/// Sets this sparse Gamma list to be a uniform distribution
		/// </summary>
		public void SetToUniform()
		{
			CommonValue.SetToUniform();
			SparseValues.Clear();
		}

		/// <summary>
		/// Asks whether this sparse Gamma list is uniform
		/// </summary>
		public bool IsUniform()
		{
			return All(x => IsUniform());
		}

		/// <summary>
		/// Evaluates the log of the density function
		/// </summary>
		/// <param name="value">The point at which to evaluate the density</param>
		/// <returns></returns>
		public double GetLogProb(SparseVector value)
		{
			if (value.Count != Count) throw new ArgumentException("Point value had incorrect size of " + value.Count + " != " + Count);
			return Reduce<double, double>(0.0, value, (x, y, z) => x + y.GetLogProb(z));
		}

		/// <summary>
		/// Samples from this sparse Gamma list
		/// </summary>
		/// <returns></returns>
		/// <remarks>This method is inefficient in that the result will be dense even though the return type is sparse.</remarks>
		public SparseVector Sample()
		{
			SparseVector sv = SparseVector.Zero(Count);
			return Sample(sv);
		}

		/// <summary>
		/// Samples from this sparse Gamma list
		/// </summary>
		/// <param name="result">Where to put the result</param>
		/// <returns></returns>
		/// <remarks>This method is inefficient in that the result will be dense even though the return type is sparse.</remarks>
		public SparseVector Sample(SparseVector result)
		{
			if (result.Count != Count) throw new ArgumentException("Result vector has incorrect size of " + result.Count + " != " + Count);
			IEnumerator<Gamma> e = GetEnumerator();
			int i = 0;
			while (e.MoveNext()) result[i++] = e.Current.Sample();
			return result;
		}

		/// <summary>
		/// Samples from a list of Gamma distributions with the specified vectors
		/// of shapes and rates
		/// </summary>
		/// <param name="shapes">Vector of shapes</param>
		/// <param name="rates">Vector of rates</param>
		/// <returns>The sample</returns>
		[Stochastic]
		public static SparseVector Sample(SparseVector shapes, SparseVector rates)
		{
			SparseVector sample = SparseVector.Zero(shapes.Count);
			sample.SetToFunction(shapes, rates, (s, r) => Gamma.Sample(s, r));
			return sample;
		}
	
		/// <summary>
		/// Sets this sparse Gamma list to another sparse Gamma list
		/// </summary>
		/// <param name="value"></param>
		public void SetTo(SparseGammaList value)
		{
			base.SetTo(value);
		}

		/// <summary>
		/// Clones this sparse Gamma list.
		/// </summary>
		public new object Clone()
		{
			return new SparseGammaList(this);
		}

		/// <summary>
		/// Sets this sparse Gamma list to the product of two other sparse Gamma lists
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public void SetToProduct(SparseGammaList a, SparseGammaList b)
		{
			SetToFunction(a, b, (ga, gb) => ga * gb);
		}
		/// <summary>
		/// Sets this sparse Gamma list to the power of another sparse Gamma list
		/// </summary>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		public void SetToPower(SparseGammaList value, double exponent)
		{
			SetToFunction(value, g => g ^ exponent);
		}

		/// <summary>
		/// Sets this sparse Gamma list to the ratio of two other sparse Gamma lists
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(SparseGammaList numerator, SparseGammaList denominator)
		{
			SetToFunction(numerator, denominator, (gn, gd) => gn / gd);
		}

		/// <summary>
		/// Creates a sparse Gamma list whose elements match the means and variances
		/// of the weighted sums of the elements of two other sparse Gamma lists.
		/// </summary>
		/// <param name="weight1">The first weight</param>
		/// <param name="value1">The first sparse Gamma list</param>
		/// <param name="weight2">The second weight</param>
		/// <param name="value2">The second sparse Gamma list</param>
		public void SetToSum(double weight1, SparseGammaList value1, double weight2, SparseGammaList value2)
		{
			SetToFunction(value1, value2, (g1, g2) => { var g = new Gamma(); g.SetToSum(weight1, g1, weight2, g2); return g; });
		}

		/// <summary>
		/// Returns the log of the integral of the product of this sparse Gamma list and another sparse Gamma list
		/// </summary>
		/// <param name="that"></param>
		public double GetLogAverageOf(SparseGammaList that)
		{
			return Reduce<double, Gamma>(0.0, that, (x, y, z) => x + y.GetLogAverageOf(z));
		}

		/// <summary>
		/// Returns the log of the integral of the product of this sparse Gamma list and another sparse Gamma list raised to a power
		/// </summary>
		/// <param name="that"></param>
		public double GetLogAverageOfPower(SparseGammaList that, double power)
		{
			return Reduce<double, Gamma>(0.0, that, (x, y, z) => x + y.GetLogAverageOfPower(z, power));
		}

		/// <summary>
		/// The expected logarithm of that sparse Gamma list under this sparse Gamma list.
		/// </summary>
		/// <param name="that"></param>
		/// <returns></returns>
		public double GetAverageLog(SparseGammaList that)
		{
			return Reduce<double, Gamma>(0.0, that, (x, y, z) => x + y.GetAverageLog(z));
		}
	}
}
