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
	/// Represents a sparse list of Gaussian distributions, optimised for the case where many share 
	/// the same parameterisation.
	/// </summary>
	[Serializable]
	[Quality(QualityBand.Experimental)]
	public class SparseGaussianList : ApproximateSparseList<Gaussian>,
		IDistribution<IList<double>>, Sampleable<IList<double>>,
		SettableTo<SparseGaussianList>, SettableToProduct<SparseGaussianList>,
		SettableToPower<SparseGaussianList>, SettableToRatio<SparseGaussianList>,
		SettableToWeightedSum<SparseGaussianList>, CanGetLogAverageOf<SparseGaussianList>, CanGetLogAverageOfPower<SparseGaussianList>,
		CanGetAverageLog<SparseGaussianList>
	{
		/// <summary>
		/// Default tolerance for sparsity approximation
		/// </summary>
		public new static double DefaultTolerance = 0.01;
		/// <summary>
		/// Default constructor
		/// </summary>
		public SparseGaussianList() : this(DefaultTolerance) { }

		/// <summary>
		/// The dimension of the SparseGaussianList domain
		/// </summary>
		public int Dimension
		{
			get { return Count; }
		}
	
		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size">The size of the list</param>
		protected SparseGaussianList(int size) : this(size, DefaultTolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian list of the given size, and assigns all
		/// elements to the specified common value
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		protected SparseGaussianList(int size, Gaussian commonValue) : this(size, commonValue, DefaultTolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		protected SparseGaussianList(int size, Gaussian commonValue, List<SparseElement<Gaussian>> sortedSparseValues) :
			this(size, commonValue, sortedSparseValues, DefaultTolerance) { }

		/// <summary>
		/// Copy constructor
		/// </summary>
		/// <param name="that"></param>
		protected SparseGaussianList(SparseGaussianList that) : base(that) { }

		/// <summary>
		/// Constructs a sparse Gaussian list with the given tolerance
		/// </summary>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public SparseGaussianList(double tolerance) : base(tolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		protected SparseGaussianList(int size, double tolerance) : base(size, Gaussian.Uniform(), tolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian list of the given size, and assigns all
		/// elements to the specified common value
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		protected SparseGaussianList(int size, Gaussian commonValue, double tolerance) : base(size, commonValue, tolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		protected SparseGaussianList(int size, Gaussian commonValue, List<SparseElement<Gaussian>> sortedSparseValues, double tolerance) :
			base(size, commonValue, sortedSparseValues, tolerance) { }

		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size"></param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		/// <returns></returns>
		public new static SparseGaussianList FromSize(int size, double tolerance)
		{
			return new SparseGaussianList(size, tolerance);
		}

		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to the specified Gaussian
		/// </summary>
		/// <param name="size">Size</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		/// <returns></returns>
		public static SparseGaussianList FromGaussian(int size, Gaussian commonValue, double tolerance)
		{
			return new SparseGaussianList(size, commonValue, tolerance);
		}

		/// <summary>
		/// Returns a sparse Gaussian list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		[Construction("Count", "CommonValue", "SparseValues", "Tolerance")]
		public new static SparseGaussianList FromSparseValues(int size, Gaussian commonValue,
			List<SparseElement<Gaussian>> sortedSparseValues, double tolerance)
		{
			return new SparseGaussianList(size, commonValue, sortedSparseValues, tolerance);
		}

		/// <summary>
		///Creates a sparse Gaussian list of a given size, with each
		/// element having a given mean and variance.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="variance">The desired variance.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGaussianList FromMeanAndVariance(
			int size, double mean, double variance, double tolerance)
		{
			return new SparseGaussianList(size, Gaussian.FromMeanAndVariance(mean, variance), tolerance);
		}

		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with each
		/// element having a given mean and precision.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="precision">The desired precision.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGaussianList FromMeanAndPrecision(
			int size, double mean, double precision, double tolerance)
		{
			return new SparseGaussianList(size, Gaussian.FromMeanAndPrecision(mean, precision), tolerance);
		}

		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with given mean and precision vectors.
		/// </summary>
		/// <param name="mean">The desired mean.</param>
		/// <param name="precision">The desired precision.</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGaussianList FromMeanAndPrecision(
			IList<double> mean, IList<double> precision, double tolerance)
		{
			var result = new SparseGaussianList(mean.Count, tolerance);
			result.SetToFunction(mean, precision, (m, p) => Gaussian.FromMeanAndPrecision(m, p));
			return result;
		}
	
		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with each
		/// element having the specified natural parameters
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="meanTimesPrecision">Mean time precision</param>
		/// <param name="precision">Precision</param>
		/// <param name="tolerance">The tolerance for the approximation</param>
		public static SparseGaussianList FromNatural(int size, double meanTimesPrecision, double precision, double tolerance)
		{
			return new SparseGaussianList(size, Gaussian.FromNatural(meanTimesPrecision, precision), tolerance);
		}

		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to uniform
		/// </summary>
		/// <param name="size"></param>
		/// <returns></returns>
		public new static SparseGaussianList FromSize(int size)
		{
			return FromSize(size, DefaultTolerance);
		}

		/// <summary>
		/// Constructs a sparse Gaussian with the specified number of elements
		/// all of which are set to the specified Gaussian
		/// </summary>
		/// <param name="size">Size</param>
		/// <param name="commonValue">The common value</param>
		/// <returns></returns>
		public static SparseGaussianList FromGaussian(int size, Gaussian commonValue)
		{
			return FromGaussian(size, commonValue, DefaultTolerance);
		}

		/// <summary>
		/// Returns a sparse Gaussian list of a given length and assigns all elements the given value,
		/// except for the specified list of sparse values. This list is stored internally as is
		/// so MUST be sorted by index and must not be modified externally after being passed in.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="commonValue">The common value</param>
		/// <param name="sortedSparseValues">The sorted list of non-common values</param>
		public new static SparseGaussianList FromSparseValues(int size, Gaussian commonValue,
			List<SparseElement<Gaussian>> sortedSparseValues)
		{
			return FromSparseValues(size, commonValue, sortedSparseValues, DefaultTolerance);
		}

		/// <summary>
		///Creates a sparse Gaussian list of a given size, with each
		/// element having a given mean and variance.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="variance">The desired variance.</param>
		public static SparseGaussianList FromMeanAndVariance(
			int size, double mean, double variance)
		{
			return FromMeanAndVariance(size, mean, variance, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with each
		/// element having a given mean and precision.
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="mean">The desired mean.</param>
		/// <param name="precision">The desired precision.</param>
		public static SparseGaussianList FromMeanAndPrecision(
			int size, double mean, double precision)
		{
			return FromMeanAndPrecision(size, mean, precision, DefaultTolerance);
		}

		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with given mean and precision vectors.
		/// </summary>
		/// <param name="mean">The desired mean.</param>
		/// <param name="precision">The desired precision.</param>
		public static SparseGaussianList FromMeanAndPrecision(
			IList<double> mean, IList<double> precision)
		{
			return FromMeanAndPrecision(mean, precision, DefaultTolerance);
		}
	
		/// <summary>
		/// Creates a sparse Gaussian list of a given size, with each
		/// element having the specified natural parameters
		/// </summary>
		/// <param name="size">The size of the list</param>
		/// <param name="meanTimesPrecision">Mean time precision</param>
		/// <param name="precision">Precision</param>
		public static SparseGaussianList FromNatural(int size, double meanTimesPrecision, double precision)
		{
			return FromNatural(size, meanTimesPrecision, precision, DefaultTolerance);
		}

		/// <summary>
		/// Sets/gets the instance as a point mass
		/// </summary>
		[NonSerializedProperty, System.Xml.Serialization.XmlIgnore]
		public IList<double> Point
		{
			get
			{
				return SparseVector.FromSparseValues(
					Count, CommonValue.Point, 
					SparseValues.Select(x => new SparseElement<double>(x.Index, x.Value.Point)).ToList());
			}
			set
			{
				ISparseList<double> p;
				if (!(value is ISparseList<double>))
					p = SparseList<double>.Copy(value);
				else p = (ISparseList<double>)value;
				var sen = p.GetSparseEnumerator();
				CommonValue = Gaussian.PointMass(sen.CommonValue);
				SparseValues.Clear();
				while (sen.MoveNext())
					SparseValues.Add(new SparseElement<Gaussian>(sen.CurrentIndex, Gaussian.PointMass(sen.Current)));
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
		/// Returns the maximum difference between the parameters of this sparse Gaussian list
		/// and another
		/// </summary>
		/// <param name="thatd">The other sparse Gaussian list</param>
		/// <returns>The maximum difference</returns>
		/// <remarks><c>a.MaxDiff(b) == b.MaxDiff(a)</c></remarks>
		public double MaxDiff(object thatd)
		{
			SparseGaussianList that = (SparseGaussianList)thatd;
			return Reduce<double, Gaussian>(
				double.NegativeInfinity, that, (x, y, z) => Math.Max(x, y.MaxDiff(z)));
		}

		/// <summary>
		/// Sets this sparse Gaussian list to be a uniform distribution
		/// </summary>
		public void SetToUniform()
		{
			CommonValue.SetToUniform();
			SparseValues.Clear();
		}

		/// <summary>
		/// Asks whether this sparse Gaussian list is uniform
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
		public double GetLogProb(IList<double> value)
		{
			if (value.Count != Count) throw new ArgumentException("Point value had incorrect size of " + value.Count + " != " + Count);
			return Reduce<double, double>(0.0, value, (x, y, z) => x + y.GetLogProb(z));
		}

		/// <summary>
		/// Samples from this sparse Gaussian list
		/// </summary>
		/// <returns></returns>
		/// <remarks>This method is inefficient in that the result will be dense even though the return type is sparse.</remarks>
		public IList<double> Sample()
		{
			SparseList<double> sv = SparseList<double>.FromSize(Count);
			return Sample(sv);
		}

		/// <summary>
		/// Samples from this sparse Gaussian list
		/// </summary>
		/// <param name="result">Where to put the result</param>
		/// <returns></returns>
		/// <remarks>This method is inefficient in that the result will be dense even though the return type is sparse.</remarks>
		public IList<double> Sample(IList<double> result)
		{
			if (result.Count != Count) throw new ArgumentException("Result list has incorrect size of " + result.Count + " != " + Count);
			IEnumerator<Gaussian> e = GetEnumerator();
			int i = 0;
			while (e.MoveNext()) result[i++] = e.Current.Sample();
			return result;
		}

		/// <summary>
		/// Samples from a list of Gaussian distributions with the specified vectors
		/// of means and precisions
		/// </summary>
		/// <param name="means">Vector of means</param>
		/// <param name="precs">Vector of precisions</param>
		/// <returns>The sample</returns>
		[Stochastic]
		public static IList<double> Sample(IList<double> means, IList<double> precs)
		{
			SparseList<double> sample = SparseList<double>.FromSize(means.Count);
			sample.SetToFunction(means, precs, (m, p) => Gaussian.Sample(m, p));
			return sample;
		}
	
		/// <summary>
		/// Sets this sparse Gaussian list to another sparse Gaussian list
		/// </summary>
		/// <param name="value"></param>
		public void SetTo(SparseGaussianList value)
		{
			base.SetTo(value);
		}

		/// <summary>
		/// Clones this sparse Gaussian list.
		/// </summary>
		public new object Clone()
		{
			return new SparseGaussianList(this);
		}

		/// <summary>
		/// Sets this sparse Gaussian list to the product of two other sparse Gaussian lists
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public void SetToProduct(SparseGaussianList a, SparseGaussianList b)
		{
			SetToFunction(a, b, (ga, gb) => ga * gb);
		}
		/// <summary>
		/// Sets this sparse Gaussian list to the power of another sparse Gaussian list
		/// </summary>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		public void SetToPower(SparseGaussianList value, double exponent)
		{
			SetToFunction(value, g => g ^ exponent);
		}

		/// <summary>
		/// Sets this sparse Gaussian list to the ratio of two other sparse Gaussian lists
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(SparseGaussianList numerator, SparseGaussianList denominator)
		{
			SetToFunction(numerator, denominator, (gn, gd) => gn / gd);
		}

		/// <summary>
		/// Creates a sparse Gaussian list whose elements match the means and variances
		/// of the weighted sums of the elements of two other sparse Gaussian lists.
		/// </summary>
		/// <param name="weight1">The first weight</param>
		/// <param name="value1">The first sparse Gaussian list</param>
		/// <param name="weight2">The second weight</param>
		/// <param name="value2">The second sparse Gaussian list</param>
		public void SetToSum(double weight1, SparseGaussianList value1, double weight2, SparseGaussianList value2)
		{
			SetToFunction(value1, value2, (g1, g2) => { var g = new Gaussian(); g.SetToSum(weight1, g1, weight2, g2); return g; });
		}

		/// <summary>
		/// Returns the log of the integral of the product of this sparse Gaussian list and another sparse Gaussian list
		/// </summary>
		/// <param name="that"></param>
		public double GetLogAverageOf(SparseGaussianList that)
		{
			return Reduce<double, Gaussian>(0.0, that, (x, y, z) => x + y.GetLogAverageOf(z));
		}

		/// <summary>
		/// Returns the log of the integral of the product of this sparse Gaussian list and another sparse Gaussian list raised to a power
		/// </summary>
		/// <param name="that"></param>
		public double GetLogAverageOfPower(SparseGaussianList that, double power)
		{
			return Reduce<double, Gaussian>(0.0, that, (x, y, z) => x + y.GetLogAverageOfPower(z, power));
		}

		/// <summary>
		/// The expected logarithm of that sparse Gaussian list under this sparse Gaussian list.
		/// </summary>
		/// <param name="that"></param>
		/// <returns></returns>
		public double GetAverageLog(SparseGaussianList that)
		{
			return Reduce<double, Gaussian>(0.0, that, (x, y, z) => x + y.GetAverageLog(z));
		}
	}
}
