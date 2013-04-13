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
	/// Provides outgoing messages for <see cref="Vector.Subvector(Vector, int, int)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Vector), "Subvector", typeof(Vector), typeof(int), typeof(int))]
	[Buffers("SourceMean", "SourceVariance")]
	[Quality(QualityBand.Preview)]
	public static class SubvectorOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="source">Constant value for 'source'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector subvector, Vector source, int startIndex)
		{
			for (int i = 0; i < subvector.Count; i++) {
				if (subvector[i] != source[startIndex+i]) return Double.NegativeInfinity;
			}
			return 0.0;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="source">Constant value for 'source'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector subvector, Vector source, int startIndex) { return LogAverageFactor(subvector, source, startIndex); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="source">Constant value for 'source'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Vector subvector, Vector source, int startIndex) { return LogAverageFactor(subvector, source, startIndex); }

		/// <summary>
		/// Initialise the buffer 'SourceVariance'
		/// </summary>
		/// <param name="Source">Incoming message from 'source'.</param>
		/// <returns>Initial value of buffer 'SourceVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix SourceVarianceInit([IgnoreDependency] VectorGaussian Source)
		{
			return new PositiveDefiniteMatrix(Source.Dimension, Source.Dimension);
		}
		/// <summary>
		/// Update the buffer 'SourceVariance'
		/// </summary>
		/// <param name="Source">Incoming message from 'source'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Source"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix SourceVariance([Proper] VectorGaussian Source, PositiveDefiniteMatrix result)
		{
			return Source.GetVariance(result);
		}

		/// <summary>
		/// Initialise the buffer 'SourceMean'
		/// </summary>
		/// <param name="Source">Incoming message from 'source'.</param>
		/// <returns>Initial value of buffer 'SourceMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector SourceMeanInit([IgnoreDependency] VectorGaussian Source)
		{
			return Vector.Zero(Source.Dimension);
		}
		/// <summary>
		/// Update the buffer 'SourceMean'
		/// </summary>
		/// <param name="Source">Incoming message from 'source'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SourceVariance">Buffer 'SourceVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Source"/> is not a proper distribution</exception>
		public static Vector SourceMean([Proper] VectorGaussian Source, [Fresh] PositiveDefiniteMatrix SourceVariance, Vector result)
		{
			return Source.GetMean(result, SourceVariance);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="SourceMean">Buffer 'SourceMean'.</param>
		/// <param name="SourceVariance">Buffer 'SourceVariance'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector subvector, [Fresh] Vector SourceMean, [Fresh] PositiveDefiniteMatrix SourceVariance, int startIndex)
		{
			double sum = 0.0;
			for (int i = startIndex; i < SourceMean.Count; i++) {
				sum += Gaussian.GetLogProb(subvector[i], SourceMean[i], SourceVariance[i, i]);
			}
			return sum;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="SourceMean">Buffer 'SourceMean'.</param>
		/// <param name="SourceVariance">Buffer 'SourceVariance'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector subvector, [Fresh] Vector SourceMean, [Fresh] PositiveDefiniteMatrix SourceVariance, int startIndex)
		{
			return LogAverageFactor(subvector, SourceMean, SourceVariance, startIndex);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Incoming message from 'subvector'.</param>
		/// <param name="to_subvector">Outgoing message to 'subvector'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(subvector) p(subvector) factor(subvector,source,startIndex,count))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(VectorGaussian subvector, [Fresh] VectorGaussian to_subvector)
		{
			return to_subvector.GetLogAverageOf(subvector);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="subvector">Incoming message from 'subvector'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(subvector) p(subvector) factor(subvector,source,startIndex,count) / sum_subvector p(subvector) messageTo(subvector))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(VectorGaussian subvector) { return 0.0; }

		[Skip]
		public static VectorGaussian SubvectorAverageConditionalInit(int count)
		{
			return new VectorGaussian(count);
		}
		[Skip]
		public static VectorGaussian SubvectorAverageLogarithmInit(int count)
		{
			return new VectorGaussian(count);
		}

		/// <summary>
		/// EP message to 'subvector'
		/// </summary>
		/// <param name="SourceMean">Buffer 'SourceMean'.</param>
		/// <param name="SourceVariance">Buffer 'SourceVariance'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'subvector' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SubvectorAverageConditional([Fresh] Vector SourceMean, [Fresh] PositiveDefiniteMatrix SourceVariance, int startIndex, VectorGaussian result)
		{
			PositiveDefiniteMatrix subVariance = new PositiveDefiniteMatrix(result.Dimension, result.Dimension);
			subVariance.SetToSubmatrix(SourceVariance, startIndex, startIndex);
			Vector subMean = Vector.Zero(result.Dimension);
			subMean.SetToSubvector(SourceMean, startIndex);
			result.SetMeanAndVariance(subMean, subVariance);
			return result;
		}
		/// <summary>
		/// EP message to 'source'
		/// </summary>
		/// <param name="subvector">Incoming message from 'subvector'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'source' as the random arguments are varied.
		/// The formula is <c>proj[p(source) sum_(subvector) p(subvector) factor(subvector,source,startIndex,count)]/p(source)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="subvector"/> is not a proper distribution</exception>
		public static VectorGaussian SourceAverageConditional([SkipIfUniform] VectorGaussian subvector, int startIndex, VectorGaussian result)
		{
			result.MeanTimesPrecision.SetAllElementsTo(0.0);
			result.MeanTimesPrecision.SetSubvector(startIndex, subvector.MeanTimesPrecision);
			result.Precision.SetAllElementsTo(0.0);
			result.Precision.SetSubmatrix(startIndex, startIndex, subvector.Precision);
			return result;
		}
		/// <summary>
		/// EP message to 'source'
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'source' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SourceAverageConditional(Vector subvector, int startIndex, VectorGaussian result)
		{
			result.MeanTimesPrecision.SetAllElementsTo(0.0);
			result.MeanTimesPrecision.SetSubvector(startIndex, subvector);
			result.Precision.SetAllElementsTo(0.0);
			int dim = result.Dimension;
			for (int i = startIndex; i < dim; i++) {
				result.Precision[i, i] = Double.PositiveInfinity;
			}
			return result;
		}

		//-- VMP ---------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(subvector,source,startIndex,count))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'subvector'
		/// </summary>
		/// <param name="SourceMean">Buffer 'SourceMean'.</param>
		/// <param name="SourceVariance">Buffer 'SourceVariance'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'subvector' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SubvectorAverageLogarithm([Fresh] Vector SourceMean, [Fresh] PositiveDefiniteMatrix SourceVariance, int startIndex, VectorGaussian result)
		{
			return SubvectorAverageConditional(SourceMean, SourceVariance, startIndex, result);
		}
		/// <summary>
		/// VMP message to 'source'
		/// </summary>
		/// <param name="subvector">Incoming message from 'subvector'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'source' with 'subvector' integrated out.
		/// The formula is <c>sum_subvector p(subvector) factor(subvector,source,startIndex,count)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="subvector"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="subvector"/> is not a proper distribution</exception>
		public static VectorGaussian SourceAverageLogarithm([SkipIfUniform] VectorGaussian subvector, int startIndex, VectorGaussian result)
		{
			return SourceAverageConditional(subvector, startIndex, result);
		}
		/// <summary>
		/// VMP message to 'source'
		/// </summary>
		/// <param name="subvector">Constant value for 'subvector'.</param>
		/// <param name="startIndex">Constant value for 'startIndex'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'source' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SourceAverageLogarithm(Vector subvector, int startIndex, VectorGaussian result)
		{
			return SourceAverageConditional(subvector, startIndex, result);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.GetItem{double}(Vector,int)/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "GetItem<>", typeof(double), typeof(Vector), typeof(int))]
	[Buffers("ArrayMean", "ArrayVariance")]
	[Quality(QualityBand.Preview)]
	public static class VectorElementOp
	{
		/// <summary>
		/// Initialise the buffer 'ArrayVariance'
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>Initial value of buffer 'ArrayVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix ArrayVarianceInit([IgnoreDependency] VectorGaussian array)
		{
			return new PositiveDefiniteMatrix(array.Dimension, array.Dimension);
		}
		/// <summary>
		/// Update the buffer 'ArrayVariance'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix ArrayVariance([Proper] VectorGaussian array, PositiveDefiniteMatrix result)
		{
			return array.GetVariance(result);
		}

		/// <summary>
		/// Initialise the buffer 'ArrayMean'
		/// </summary>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>Initial value of buffer 'ArrayMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector ArrayMeanInit([IgnoreDependency] VectorGaussian array)
		{
			return Vector.Zero(array.Dimension);
		}
		/// <summary>
		/// Update the buffer 'ArrayMean'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Vector ArrayMean([Proper] VectorGaussian array, [Fresh] PositiveDefiniteMatrix ArrayVariance, Vector result)
		{
			return array.GetMean(result, ArrayVariance);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayMean">Buffer 'ArrayMean'.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double item, [SkipIfUniform] VectorGaussian array, [Fresh] Vector ArrayMean, [Fresh] PositiveDefiniteMatrix ArrayVariance, int index)
		{
			Gaussian to_item = ItemAverageConditional(array, ArrayMean, ArrayVariance, index);
			return to_item.GetLogProb(item);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayMean">Buffer 'ArrayMean'.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(array) p(array) factor(item,array,index))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(double item, [SkipIfUniform] VectorGaussian array, [Fresh] Vector ArrayMean, [Fresh] PositiveDefiniteMatrix ArrayVariance, int index)
		{
			return LogAverageFactor(item, array, ArrayMean, ArrayVariance, index);
		}

		// this method must have 'array' as a parameter, even though it is not used, to disambiguate from the other overload of GetItem
		/// <summary>
		/// EP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayMean">Buffer 'ArrayMean'.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>The outgoing EP message to the 'item' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[p(item) sum_(array) p(array) factor(item,array,index)]/p(item)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Gaussian ItemAverageConditional([SkipIfUniform] VectorGaussian array, [Fresh] Vector ArrayMean, [Fresh] PositiveDefiniteMatrix ArrayVariance, int index)
		{
			return new Gaussian(ArrayMean[index], ArrayVariance[index,index]);
		}

		[Skip]
		public static Gaussian ItemAverageConditionalInit()
		{
			return Gaussian.Uniform();
		}
		[Skip]
		public static Gaussian ItemAverageLogarithmInit()
		{
			return Gaussian.Uniform();
		}

		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'array' as the random arguments are varied.
		/// The formula is <c>proj[p(array) sum_(item) p(item) factor(item,array,index)]/p(array)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static VectorGaussian ArrayAverageConditional([SkipIfUniform] Gaussian item, int index, VectorGaussian result)
		{
			result.MeanTimesPrecision.SetAllElementsTo(0.0);
			result.MeanTimesPrecision[index] = item.MeanTimesPrecision;
			result.Precision.SetAllElementsTo(0.0);
			result.Precision[index, index] = item.Precision;
			return result;
		}
		/// <summary>
		/// EP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian ArrayAverageConditional(double item, int index, VectorGaussian result)
		{
			result.MeanTimesPrecision.SetAllElementsTo(0.0);
			result.MeanTimesPrecision[index] = item;
			result.Precision.SetAllElementsTo(0.0);
			result.Precision[index, index] = Double.PositiveInfinity;
			return result;
		}

		//-- VMP ------------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayMean">Buffer 'ArrayMean'.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static double AverageLogFactor(double item, [SkipIfUniform] VectorGaussian array, [Fresh] Vector ArrayMean, [Fresh] PositiveDefiniteMatrix ArrayVariance, int index)
		{
			return LogAverageFactor(item, array, ArrayMean, ArrayVariance, index);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="item">Incoming message from 'item'.</param>
		/// <param name="array">Incoming message from 'array'.</param>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// In Variational Message Passing, the evidence contribution of a deterministic factor is zero.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Gaussian item, VectorGaussian array) { return 0.0; }

		/// <summary>
		/// VMP message to 'item'
		/// </summary>
		/// <param name="array">Incoming message from 'array'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="ArrayMean">Buffer 'ArrayMean'.</param>
		/// <param name="ArrayVariance">Buffer 'ArrayVariance'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <returns>The outgoing VMP message to the 'item' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'item' as the random arguments are varied.
		/// The formula is <c>proj[sum_(array) p(array) factor(item,array,index)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="array"/> is not a proper distribution</exception>
		public static Gaussian ItemAverageLogarithm([SkipIfUniform] VectorGaussian array, [Fresh] Vector ArrayMean, [Fresh] PositiveDefiniteMatrix ArrayVariance, int index)
		{
			return ItemAverageConditional(array, ArrayMean, ArrayVariance, index);
		}
		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Incoming message from 'item'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' with 'item' integrated out.
		/// The formula is <c>sum_item p(item) factor(item,array,index)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="item"/> is not a proper distribution</exception>
		public static VectorGaussian ArrayAverageLogarithm([SkipIfUniform] Gaussian item, int index, VectorGaussian result)
		{
			return ArrayAverageConditional(item, index, result);
		}
		/// <summary>
		/// VMP message to 'array'
		/// </summary>
		/// <param name="item">Constant value for 'item'.</param>
		/// <param name="index">Constant value for 'index'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'array' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian ArrayAverageLogarithm(double item, int index, VectorGaussian result)
		{
			return ArrayAverageConditional(item, index, result);
		}
	}
}
