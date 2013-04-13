// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Vector.FromArray(double[])"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "vector", "array" }, typeof(Vector), "FromArray", typeof(double[]))]
	[Quality(QualityBand.Preview)]
	public static class VectorFromArrayOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="vector">Constant value for 'fromArray'.</param>
		/// <param name="array">Constant value for 'data'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(fromArray,data))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector vector, double[] array)
		{
			for (int i = 0; i < array.Length; i++) {
				if (vector[i] != array[i]) return Double.NegativeInfinity;
			}
			return 0.0;
		}
		public static double LogEvidenceRatio(Vector vector, double[] array) { return LogAverageFactor(vector, array); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="vector">Constant value for 'fromArray'.</param>
		/// <param name="array">Incoming message from 'data'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(data) p(data) factor(fromArray,data))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector vector, IList<Gaussian> array)
		{
			double sum = 0.0;
			for (int i = 0; i < array.Count; i++) {
				sum += array[i].GetLogProb(vector[i]);
			}
			return sum;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="vector">Constant value for 'fromArray'.</param>
		/// <param name="array">Incoming message from 'data'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(data) p(data) factor(fromArray,data))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector vector, IList<Gaussian> array)
		{
			return LogAverageFactor(vector, array);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="vector">Incoming message from 'fromArray'.</param>
		/// <param name="to_vector">Outgoing message to 'vector'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(fromArray) p(fromArray) factor(fromArray,data))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(VectorGaussian vector, [Fresh] VectorGaussian to_vector)
		{
			return to_vector.GetLogAverageOf(vector);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="vector">Incoming message from 'fromArray'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(fromArray) p(fromArray) factor(fromArray,data) / sum_fromArray p(fromArray) messageTo(fromArray))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(VectorGaussian vector) { return 0.0; }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(fromArray,data))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// EP message to 'fromArray'
		/// </summary>
		/// <param name="array">Incoming message from 'data'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'fromArray' as the random arguments are varied.
		/// The formula is <c>proj[p(fromArray) sum_(data) p(data) factor(fromArray,data)]/p(fromArray)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		// TM: SkipIfAllUniform would be more accurate but leads to half-uniform distributions
		public static VectorGaussian VectorAverageConditional([SkipIfAnyUniform] IList<Gaussian> array, VectorGaussian result)
		{
			return ArrayFromVectorOp.VectorAverageConditional(array, result);
		}
		[Skip]
		public static VectorGaussian VectorAverageConditionalInit([IgnoreDependency] IList<Gaussian> array)
		{
			return new VectorGaussian(array.Count);
		}
		[Skip]
		public static VectorGaussian VectorAverageLogarithmInit([IgnoreDependency] IList<Gaussian> array)
		{
			return new VectorGaussian(array.Count);
		}

		/// <summary>
		/// VMP message to 'fromArray'
		/// </summary>
		/// <param name="array">Incoming message from 'data'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'fromArray' as the random arguments are varied.
		/// The formula is <c>proj[sum_(data) p(data) factor(fromArray,data)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		// TM: SkipIfAllUniform would be more accurate but leads to half-uniform distributions
		public static VectorGaussian VectorAverageLogarithm([SkipIfAnyUniform] IList<Gaussian> array, VectorGaussian result)
		{
			return VectorAverageConditional(array, result);
		}

		/// <summary>
		/// EP message to 'data'
		/// </summary>
		/// <param name="vector">Incoming message from 'fromArray'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'data'.</param>
		/// <param name="to_vector">Outgoing message to 'vector'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'data' as the random arguments are varied.
		/// The formula is <c>proj[p(data) sum_(fromArray) p(fromArray) factor(fromArray,data)]/p(data)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="vector"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageConditional<GaussianList>([SkipIfUniform] VectorGaussian vector, IList<Gaussian> array, [Fresh] VectorGaussian to_vector, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			return ArrayFromVectorOp.ArrayAverageConditional<GaussianList>(array, vector, to_vector, result);
		}
		/// <summary>
		/// VMP message to 'data'
		/// </summary>
		/// <param name="vector">Incoming message from 'fromArray'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="array">Incoming message from 'data'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'data' with 'fromArray' integrated out.
		/// The formula is <c>sum_fromArray p(fromArray) factor(fromArray,data)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="vector"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageLogarithm<GaussianList>([SkipIfUniform] VectorGaussian vector, [Proper] IList<Gaussian> array, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			// prec[i] = vector[i].Prec
			// meanTimesPrec = vector.MeanTimesPrec - vector.Prec[:,noti]*array[noti].Mean
			//               = vector.MeanTimesPrec - vector.Prec*array.Mean + diag(invdiag(vector.Prec))*array.Mean
			if(result.Count != vector.Dimension) throw new ArgumentException("vector.Dimension ("+vector.Dimension+") != result.Count ("+result.Count+")");
			if(result.Count != array.Count) throw new ArgumentException("array.Count ("+array.Count+") != result.Count ("+result.Count+")");
			if(vector.IsPointMass) return ArrayAverageLogarithm(vector.Point, result);
			int length = result.Count;
			Vector mean = Vector.Zero(length);
			for (int i = 0; i < length; i++) {
				Gaussian item = array[i];			
				mean[i] = item.GetMean();				
			}
			Vector meanTimesPrecision = vector.Precision * mean;
			for (int i = 0; i < length; i++) {
				double prec = vector.Precision[i,i];
				if(Double.IsPositiveInfinity(prec)) throw new NotSupportedException("Singular VectorGaussians not supported");
				double mprec = vector.MeanTimesPrecision[i] - meanTimesPrecision[i] + prec * mean[i];
				result[i] = Gaussian.FromNatural(mprec, prec);
			}
			return result;
		}
		/// <summary>
		/// EP message to 'data'
		/// </summary>
		/// <param name="vector">Constant value for 'fromArray'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'data' conditioned on the given values.
		/// </para></remarks>
		public static GaussianList ArrayAverageConditional<GaussianList>(Vector vector, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			return ArrayAverageLogarithm<GaussianList>(vector, result);
		}
		/// <summary>
		/// VMP message to 'data'
		/// </summary>
		/// <param name="vector">Constant value for 'fromArray'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'data' conditioned on the given values.
		/// </para></remarks>
		public static GaussianList ArrayAverageLogarithm<GaussianList>(Vector vector, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			int length = result.Count;
			for (int i = 0; i < length; i++) {
				result[i] = Gaussian.PointMass(vector[i]);
			}
			return result;
		}
	}
}
