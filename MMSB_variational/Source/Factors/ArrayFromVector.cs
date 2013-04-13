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
	/// Provides outgoing messages for <see cref="Factor.ArrayFromVector"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "ArrayFromVector")]
	[Quality(QualityBand.Preview)]
	public static class ArrayFromVectorOp
	{
		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Constant value for 'vector'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_() p() factor(array,vector))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double[] array, Vector vector)
		{
			if (array.Length != vector.Count) return Double.NegativeInfinity;
			for (int i = 0; i < array.Length; i++) {
				if (array[i] != vector[i]) return Double.NegativeInfinity;
			}
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Constant value for 'vector'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(array,vector))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double[] array, Vector vector) { return LogAverageFactor(array, vector); }

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Constant value for 'vector'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		public static double AverageLogFactor(double[] array, Vector vector) { return LogAverageFactor(array, vector); }

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(vector) p(vector) factor(array,vector))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double[] array, VectorGaussian vector)
		{
			return vector.GetLogProb(Vector.FromArray(array));
		}

		/// <summary>
		/// Evidence message for EP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence.</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(vector) p(vector) factor(array,vector))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double[] array, VectorGaussian vector)
		{
			return LogAverageFactor(array, vector);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <param name="to_vector">Outgoing message to 'vector'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array,vector) p(array,vector) factor(array,vector))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(IList<Gaussian> array, VectorGaussian vector, [Fresh] VectorGaussian to_vector)
		{
			return vector.GetLogAverageOf(to_vector);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <param name="to_vector">Outgoing message to 'vector'.</param>
		/// <param name="to_array">Outgoing message to 'array'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array,vector) p(array,vector) factor(array,vector) / sum_array p(array) messageTo(array))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(IList<Gaussian> array, VectorGaussian vector, [Fresh] VectorGaussian to_vector, [Fresh] IList<Gaussian> to_array)
		{
			double result = LogAverageFactor(array, vector, to_vector);
			int length = array.Count;
			for (int i = 0; i < length; i++) {
				result -= to_array[i].GetLogAverageOf(array[i]);
			}
			return result;
		}

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(double[] array, VectorGaussian vector) { return 0.0; }

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>	
		[Skip]
		public static double AverageLogFactor(IList<Gaussian> array, VectorGaussian vector) { return 0.0; }

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="vector">Constant value for 'vector'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// </para></remarks>	
		[Skip]
		public static double AverageLogFactor(IList<Gaussian> array, Vector vector) { return 0.0; }

		/// <summary>
		/// EP message to 'vector'.
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'vector' as the random arguments are varied.
		/// The formula is <c>proj[p(vector) sum_(array) p(array) factor(array,vector)]/p(vector)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static VectorGaussian VectorAverageConditional([SkipIfAllUniform] IList<Gaussian> array, VectorGaussian result)
		{
			return VectorAverageLogarithm(array, result);
		}

		/// <summary>
		/// VMP message to 'vector'.
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'vector' with 'array' integrated out.
		/// The formula is <c>sum_array p(array) factor(array,vector)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static VectorGaussian VectorAverageLogarithm([SkipIfAllUniform] IList<Gaussian> array, VectorGaussian result)
		{
			if (result.Dimension != array.Count) throw new ArgumentException("array.Count ("+array.Count+") != result.Dimension ("+result.Dimension+")");
			result.Precision.SetAllElementsTo(0.0);
			int length = array.Count;
			for (int i = 0; i < length; i++) {
				Gaussian item = array[i];
				result.Precision[i, i] = item.Precision;
				result.MeanTimesPrecision[i] = item.MeanTimesPrecision;
			}
			return result;
		}

		/// <summary>
		/// EP message to 'vector'.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'vector' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian VectorAverageConditional(double[] array, VectorGaussian result)
		{
			return VectorAverageLogarithm(array, result);
		}

		/// <summary>
		/// VMP message to 'vector'.
		/// </summary>
		/// <param name="array">Constant value for 'array'.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'vector' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian VectorAverageLogarithm(double[] array, VectorGaussian result)
		{
			result.Point = result.Point;
			result.Point.SetTo(array);
			return result;
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <param name="vector">Incoming message from 'vector'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_vector">Outgoing message to 'vector'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(vector) p(vector) factor(array,vector)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="vector"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageConditional<GaussianList>(IList<Gaussian> array, [SkipIfUniform] VectorGaussian vector, [Fresh] VectorGaussian to_vector, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			if (result.Count != vector.Dimension) throw new ArgumentException("vector.Dimension ("+vector.Dimension+") != result.Count ("+result.Count+")");
			int length = result.Count;
			VectorGaussian vectorTimesArray = new VectorGaussian(vector.Dimension);
			vectorTimesArray.SetToProduct(vector, to_vector);
			Vector mean = Vector.Zero(length);
			PositiveDefiniteMatrix variance = new PositiveDefiniteMatrix(length, length);
			vectorTimesArray.GetMeanAndVariance(mean, variance);
			for (int i = 0; i < length; i++) {
				Gaussian marginal = Gaussian.FromMeanAndVariance(mean[i], variance[i, i]);
				marginal.SetToRatio(marginal, array[i]);
				result[i] = marginal;
			}
			return result;
		}
		[Skip]
		public static DistributionStructArray<Gaussian, double> ArrayAverageConditionalInit([IgnoreDependency] VectorGaussian vector)
		{
			return new DistributionStructArray<Gaussian, double>(vector.Dimension);
		}

		/// <summary>
		/// VMP message to 'array'.
		/// </summary>
		/// <param name="vector">Incoming message from 'vector'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message.</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[sum_(vector) p(vector) factor(array,vector)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="vector"/> is not a proper distribution</exception>
		public static GaussianList ArrayAverageLogarithm<GaussianList>([SkipIfUniform] VectorGaussian vector, GaussianList result)
			where GaussianList : IList<Gaussian>
		{
			if (result.Count != vector.Dimension) throw new ArgumentException("vector.Dimension ("+vector.Dimension+") != result.Count ("+result.Count+")");
			int length = result.Count;
			Vector mean = Vector.Zero(length);
			PositiveDefiniteMatrix variance = new PositiveDefiniteMatrix(length, length);
			vector.GetMeanAndVariance(mean, variance);
			for (int i = 0; i < length; i++) {
				result[i] = Gaussian.FromMeanAndVariance(mean[i], variance[i, i]);
			}
			return result;
		}
	}
}
