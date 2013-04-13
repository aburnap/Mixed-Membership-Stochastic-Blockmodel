// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Distributions
{
	/// <summary>
	/// A point mass, which is the 'distribution' you get for an observed variable.
	/// All the probability mass is at the point given by observed value.
	/// </summary>
	[Serializable]
	[Quality(QualityBand.Stable)]
	public class PointMass<T> : IDistribution<T>, Sampleable<T>, CanGetMean<T>, SettableTo<PointMass<T>>, SettableToProduct<PointMass<T>>,
			SettableToRatio<PointMass<T>>, SettableToPower<PointMass<T>>, SettableToWeightedSum<PointMass<T>>, CanGetLogAverageOf<PointMass<T>>,
			CanGetAverageLog<PointMass<T>>
	{
		/// <summary>
		/// Parameterless constructor required for serialization 
		/// </summary>
		private PointMass() { }

		/// <summary>
		/// Creates a point mass at the specified location.
		/// </summary>
		/// <param name="point">The location of the point mass.</param>
		public PointMass(T point)
		{
			Point = point;
		}

		/// <summary>
		/// Creates a point mass at the specified location.
		/// </summary>
		/// <param name="point">The location of the point mass.</param>
		public static PointMass<T> Create(T point)
		{
			return new PointMass<T>(point);
		}

		/// <summary>
		/// The location of the point mass.
		/// </summary>
		public T Point { get; set; }

		/// <summary>
		/// Creates a copy of the point mass.
		/// </summary>
		/// <returns>The new PointMass object</returns>
		public object Clone()
		{
			return new PointMass<T>(Point);
		}

		/// <summary>
		/// Always returns true.
		/// </summary>
		public bool IsPointMass
		{
			get { return true; }
		}

		#region Diffable Members
		/// <summary>
		/// Returns 0 if the two distributions are the same, positive infinity otherwise.
		/// </summary>
		/// <param name="that"></param>
		/// <returns></returns>
		public double MaxDiff(object that)
		{
			if(!(that is PointMass<T>)) return Double.PositiveInfinity;
			PointMass<T> thatd = (PointMass<T>)that;
			return (this.Point.Equals(thatd.Point)) ? 0.0 : Double.PositiveInfinity;
		}
		#endregion

		#region SettableToUniform Members

		/// <summary>
		/// Always throws an exception, since a PointMass cannot be set to uniform.
		/// </summary>
		public void SetToUniform()
		{
			throw new InvalidOperationException("A PointMass cannot be set to uniform.");
		}

		/// <summary>
		/// Always returns false
		/// </summary>
		/// <returns>false</returns>
		public bool IsUniform()
		{
			return false;
		}

		#endregion

		#region CanGetLogProb<T> Members

		/// <summary>
		/// Returns 0 if the value is at the point mass and negative infinity elsewhere.
		/// </summary>
		/// <param name="value">The value at which to compute the log probability.</param>
		/// <returns>0 or negative infinity</returns>
		public double GetLogProb(T value)
		{
			if (value.Equals(Point)) return 0;
			return Double.NegativeInfinity;
		}

		#endregion

		/// <summary>
		/// 
		/// </summary>
		/// <returns></returns>
		public override string ToString()
		{
			return "PointMass("+Point+")";
		}

		/// <summary>
		/// Returns the location of the point mass
		/// </summary>
		/// <returns></returns>
		public T Sample()
		{
			return Point;
		}

		/// <summary>
		/// Returns the location of the point mass
		/// </summary>
		/// <param name="result">Where to put the result</param>
		/// <returns></returns>
		public T Sample(T result)
		{
			return Sample();
		}

		/// <summary>
		/// Returns the location of the point mass
		/// </summary>
		/// <returns></returns>
		public T GetMean()
		{
			return Point;
		}

		/// <summary>
		/// Sets this point mass to that point mass
		/// </summary>
		/// <param name="value"></param>
		public void SetTo(PointMass<T> value)
		{
			Point = value.Point;
		}

		/// <summary>
		/// Throws an exception unless the two point masses are equal 
		/// </summary>
		/// <param name="a"></param>
		/// <param name="b"></param>
		public void SetToProduct(PointMass<T> a, PointMass<T> b)
		{
			if(!a.Point.Equals(b.Point)) throw new AllZeroException();
			SetTo(a);
		}

		/// <summary>
		/// Not supported
		/// </summary>
		/// <param name="numerator"></param>
		/// <param name="denominator"></param>
		public void SetToRatio(PointMass<T> numerator, PointMass<T> denominator)
		{
			throw new NotSupportedException();
		}

		/// <summary>
		/// Sets to the value of the given point mass
		/// </summary>
		/// <param name="value"></param>
		/// <param name="exponent"></param>
		public void SetToPower(PointMass<T> value, double exponent)
		{
			if(exponent <= 0) throw new NotSupportedException();
			SetTo(value);
		}

		/// <summary>
		/// Throws an exception unless the point masses are equal
		/// </summary>
		/// <param name="weight1"></param>
		/// <param name="value1"></param>
		/// <param name="weight2"></param>
		/// <param name="value2"></param>
		public void SetToSum(double weight1, PointMass<T> value1, double weight2, PointMass<T> value2)
		{
			if(!value1.Equals(value2)) throw new NotSupportedException();
			SetTo(value1);
		}

		/// <summary>
		/// Returns 0 if the this and that point mass are equal, negative infinity otherwise
		/// </summary>
		/// <param name="that"></param>
		/// <returns>Not implemented</returns>
		public double GetLogAverageOf(PointMass<T> that)
		{
			if(Point.Equals(that.Point)) return 0.0;
			else return Double.NegativeInfinity;
		}

		/// <summary>
		/// Returns 0 if the this and that point mass are equal, negative infinity otherwise
		/// </summary>
		/// <param name="that"></param>
		/// <returns>Not implemented</returns>
		public double GetAverageLog(PointMass<T> that)
		{
			return GetLogAverageOf(that);
		}
	}
}
