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
	/// Provides outgoing messages for <see cref="Vector.Concat"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Vector), "Concat", typeof(Vector), typeof(Vector))]
	[Quality(QualityBand.Preview)]
	public static class ConcatOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector concat, Vector first, Vector second)
		{
			for (int i = 0; i < first.Count; i++) {
				if (concat[i] != first[i]) return Double.NegativeInfinity;
			}
			int dim1 = first.Count;
			for (int i = 0; i < second.Count; i++) {
				if (concat[i + dim1] != second[i]) return Double.NegativeInfinity;
			}
			return 0.0;
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(concat) p(concat) factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(VectorGaussian concat, Vector first, Vector second)
		{
			return concat.GetLogProb(Vector.Concat(first, second));
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(first,second) p(first,second) factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector concat, VectorGaussian first, VectorGaussian second)
		{
			Vector concat1 = Vector.Subvector(concat, 0, first.Dimension);
			Vector concat2 = Vector.Subvector(concat, first.Dimension, second.Dimension);
			return first.GetLogProb(concat1) + second.GetLogProb(concat2);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(second) p(second) factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector concat, Vector first, VectorGaussian second)
		{
			for (int i = 0; i < first.Count; i++) {
				if (concat[i] != first[i]) return Double.NegativeInfinity;
			}
			Vector concat2 = Vector.Subvector(concat, first.Count, second.Dimension);
			return second.GetLogProb(concat2);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(first) p(first) factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector concat, VectorGaussian first, Vector second)
		{
			int dim1 = first.Dimension;
			for (int i = 0; i < second.Count; i++) {
				if (concat[i + dim1] != second[i]) return Double.NegativeInfinity;
			}
			Vector concat1 = Vector.Subvector(concat, 0, first.Dimension);
			return first.GetLogProb(concat1);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(first,second) p(first,second) factor(concat,first,second))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector concat, VectorGaussian first, VectorGaussian second)
		{
			return LogAverageFactor(concat, first, second);
		}

		public static double LogEvidenceRatio(Vector concat, Vector first, Vector second)
		{
			return LogAverageFactor(concat, first, second);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(second) p(second) factor(concat,first,second))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector concat, Vector first, VectorGaussian second)
		{
			return LogAverageFactor(concat, first, second);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(first) p(first) factor(concat,first,second))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector concat, VectorGaussian first, Vector second)
		{
			return LogAverageFactor(concat, first, second);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'.</param>
		/// <param name="to_concat">Outgoing message to 'concat'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(concat) p(concat) factor(concat,first,second))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(VectorGaussian concat, [Fresh] VectorGaussian to_concat)
		{
			return to_concat.GetLogAverageOf(concat);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(concat) p(concat) factor(concat,first,second) / sum_concat p(concat) messageTo(concat))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(VectorGaussian concat) { return 0.0; }

		[Skip]
		public static VectorGaussian ConcatAverageConditionalInit([IgnoreDependency] VectorGaussian first, [IgnoreDependency] VectorGaussian second)
		{
			return new VectorGaussian(first.Dimension+second.Dimension);
		}
		[Skip]
		public static VectorGaussian ConcatAverageConditionalInit([IgnoreDependency] Vector first, [IgnoreDependency] VectorGaussian second)
		{
			return new VectorGaussian(first.Count+second.Dimension);
		}
		[Skip]
		public static VectorGaussian ConcatAverageConditionalInit([IgnoreDependency] VectorGaussian first, [IgnoreDependency] Vector second)
		{
			return new VectorGaussian(first.Dimension+second.Count);
		}

		/// <summary>
		/// EP message to 'concat'
		/// </summary>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[p(concat) sum_(first,second) p(first,second) factor(concat,first,second)]/p(concat)</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static VectorGaussian ConcatAverageConditional(VectorGaussian first, VectorGaussian second, VectorGaussian result)
		{
			int dim1 = first.Dimension;
			int dim2 = second.Dimension;
			if (result.Dimension != dim1 + dim2) throw new ArgumentException("concat.Dimension ("+result.Dimension+") != first.Dimension ("+first.Dimension+") + second.Dimension ("+second.Dimension+")");
			// assume result.Precision was initialized to 0.0?
			result.Precision.SetAllElementsTo(0.0);
			result.Precision.SetSubmatrix(0, 0, first.Precision);
			result.Precision.SetSubmatrix(dim1, dim1, second.Precision);
			result.MeanTimesPrecision.SetSubvector(0, first.MeanTimesPrecision);
			result.MeanTimesPrecision.SetSubvector(dim1, second.MeanTimesPrecision);
			return result;
		}

		/// <summary>
		/// EP message to 'concat'
		/// </summary>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[p(concat) sum_(second) p(second) factor(concat,first,second)]/p(concat)</c>.
		/// </para></remarks>
		public static VectorGaussian ConcatAverageConditional(Vector first, VectorGaussian second, VectorGaussian result)
		{
			int dim1 = first.Count;
			result.Precision.SetAllElementsTo(0.0);
			for (int i = 0; i < dim1; i++) {
				result.Precision[i, i] = Double.PositiveInfinity;
				result.MeanTimesPrecision[i] = first[i];
			}
			result.Precision.SetSubmatrix(dim1, dim1, second.Precision);
			result.MeanTimesPrecision.SetSubvector(dim1, second.MeanTimesPrecision);
			return result;
		}

		/// <summary>
		/// EP message to 'concat'
		/// </summary>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[p(concat) sum_(first) p(first) factor(concat,first,second)]/p(concat)</c>.
		/// </para></remarks>
		public static VectorGaussian ConcatAverageConditional(VectorGaussian first, Vector second, VectorGaussian result)
		{
			int dim1 = first.Dimension;
			int dim2 = second.Count;
			result.Precision.SetAllElementsTo(0.0);
			result.Precision.SetSubmatrix(0, 0, first.Precision);
			result.MeanTimesPrecision.SetSubvector(0, first.MeanTimesPrecision);
			for (int i = 0; i < dim2; i++) {
				int j = i+dim1;
				result.Precision[j, j] = Double.PositiveInfinity;
				result.MeanTimesPrecision[j] = second[i];
			}
			return result;
		}

		/// <summary>
		/// EP message to 'first'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'first' as the random arguments are varied.
		/// The formula is <c>proj[p(first) sum_(concat) p(concat) factor(concat,first,second)]/p(first)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian FirstAverageConditional([SkipIfUniform] VectorGaussian concat, Vector second, VectorGaussian result)
		{
			// joint distribution is proportional to: exp(-0.5 [first-mean1; second-mean2]' [prec11 prec12; prec21 prec22] [first-mean1; second-mean2])
			// posterior for first is proportional to: exp(-0.5 ((first-mean1)' prec11 (first-mean1) + 2 (first-mean1)' prec12 (second-mean2)))
			// = exp(-0.5 (first' prec11 first - 2 first' prec11 mean1 + 2 first' prec12 (second-mean2)))
			// first.precision = prec11
			// first.meanTimesPrecision = prec11 mean1 - prec12 (second-mean2) = [prec11; prec12] mean - prec12 second
			int dim1 = result.Dimension;
			int dim2 = second.Count;
			if (concat.Dimension != dim1 + dim2) throw new ArgumentException("concat.Dimension ("+concat.Dimension+") != first.Dimension ("+dim1+") + second.Dimension ("+dim2+")");
			result.Precision.SetToSubmatrix(concat.Precision, 0, 0);
			Matrix prec12 = new Matrix(dim1, dim2);
			prec12.SetToSubmatrix(concat.Precision, 0, dim1);
			Vector prec12second = Vector.Zero(dim1);
			prec12second.SetToProduct(prec12, second);
			result.MeanTimesPrecision.SetToSubvector(concat.MeanTimesPrecision, 0);
			result.MeanTimesPrecision.SetToDifference(result.MeanTimesPrecision, prec12second);
			return result;
		}

		/// <summary>
		/// EP message to 'first'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'first' as the random arguments are varied.
		/// The formula is <c>proj[p(first) sum_(concat,second) p(concat,second) factor(concat,first,second)]/p(first)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian FirstAverageConditional([SkipIfUniform] VectorGaussian concat, VectorGaussian second, VectorGaussian result)
		{
			if (second.IsPointMass) return FirstAverageConditional(concat, second.Point, result);
			int dim1 = result.Dimension;
			VectorGaussian concatTimesSecond = new VectorGaussian(concat.Dimension);
			concatTimesSecond.MeanTimesPrecision.SetSubvector(dim1, second.MeanTimesPrecision);
			concatTimesSecond.Precision.SetSubmatrix(dim1, dim1, second.Precision);
			concatTimesSecond.SetToProduct(concatTimesSecond, concat);
			concatTimesSecond.GetMarginal(0, result);
			return result;
		}

		/// <summary>
		/// EP message to 'second'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'second' as the random arguments are varied.
		/// The formula is <c>proj[p(second) sum_(concat) p(concat) factor(concat,first,second)]/p(second)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian SecondAverageConditional([SkipIfUniform] VectorGaussian concat, Vector first, VectorGaussian result)
		{
			// prec = concat.Precision[dim2,dim2]
			// meanTimesPrec = concat.MeanTimesPrecision[dim2] - concat.Precision[dim2,dim1]*first
			int dim1 = first.Count;
			int dim2 = result.Dimension;
			if (concat.Dimension != dim1 + dim2) throw new ArgumentException("concat.Dimension ("+concat.Dimension+") != first.Dimension ("+dim1+") + second.Dimension ("+dim2+")");
			result.Precision.SetToSubmatrix(concat.Precision, dim1, dim1);
			Matrix prec21 = new Matrix(dim2, dim1);
			prec21.SetToSubmatrix(concat.Precision, dim1, 0);
			Vector prec21first = Vector.Zero(dim2);
			prec21first.SetToProduct(prec21, first);
			result.MeanTimesPrecision.SetToSubvector(concat.MeanTimesPrecision, dim1);
			result.MeanTimesPrecision.SetToDifference(result.MeanTimesPrecision, prec21first);
			return result;
		}

		/// <summary>
		/// EP message to 'second'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'second' as the random arguments are varied.
		/// The formula is <c>proj[p(second) sum_(concat,first) p(concat,first) factor(concat,first,second)]/p(second)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian SecondAverageConditional([SkipIfUniform] VectorGaussian concat, VectorGaussian first, VectorGaussian result)
		{
			if (first.IsPointMass) return SecondAverageConditional(concat, first.Point, result);
			int dim1 = first.Dimension;
			VectorGaussian concatTimesFirst = new VectorGaussian(concat.Dimension);
			concatTimesFirst.MeanTimesPrecision.SetSubvector(0, first.MeanTimesPrecision);
			concatTimesFirst.Precision.SetSubmatrix(0, 0, first.Precision);
			concatTimesFirst.SetToProduct(concatTimesFirst, concat);
			concatTimesFirst.GetMarginal(dim1, result);
			return result;
		}

		/// <summary>
		/// EP message to 'first'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'first' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian FirstAverageConditional(Vector concat, VectorGaussian result)
		{
			result.Precision.SetAllElementsTo(0.0);
			for (int i = 0; i < result.Dimension; i++) {
				result.MeanTimesPrecision[i] = concat[i];
				result.Precision[i, i] = Double.PositiveInfinity;
			}
			return result;
		}

		/// <summary>
		/// EP message to 'second'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'second' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SecondAverageConditional(Vector concat, Vector first, VectorGaussian result)
		{
			int dim1 = first.Count;
			result.Precision.SetAllElementsTo(0.0);
			for (int i = 0; i < result.Dimension; i++) {
				result.MeanTimesPrecision[i] = concat[i + dim1];
				result.Precision[i, i] = Double.PositiveInfinity;
			}
			return result;
		}

		/// <summary>
		/// EP message to 'second'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'second' as the random arguments are varied.
		/// The formula is <c>proj[p(second) sum_(first) p(first) factor(concat,first,second)]/p(second)</c>.
		/// </para></remarks>
		public static VectorGaussian SecondAverageConditional(Vector concat, VectorGaussian first, VectorGaussian result)
		{
			int dim1 = first.Dimension;
			result.Precision.SetAllElementsTo(0.0);
			for (int i = 0; i < result.Dimension; i++) {
				result.MeanTimesPrecision[i] = concat[i + dim1];
				result.Precision[i, i] = Double.PositiveInfinity;
			}
			return result;
		}

		// VMP //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(concat,first,second))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		[Skip]
		public static VectorGaussian ConcatAverageLogarithmInit([IgnoreDependency] VectorGaussian first, [IgnoreDependency] VectorGaussian second)
		{
			return new VectorGaussian(first.Dimension+second.Dimension);
		}
		[Skip]
		public static VectorGaussian ConcatAverageLogarithmInit([IgnoreDependency] Vector first, [IgnoreDependency] VectorGaussian second)
		{
			return new VectorGaussian(first.Count+second.Dimension);
		}
		[Skip]
		public static VectorGaussian ConcatAverageLogarithmInit([IgnoreDependency] VectorGaussian first, [IgnoreDependency] Vector second)
		{
			return new VectorGaussian(first.Dimension+second.Count);
		}

		/// <summary>
		/// VMP message to 'concat'
		/// </summary>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[sum_(first,second) p(first,second) factor(concat,first,second)]</c>.
		/// </para></remarks>
		[SkipIfAllUniform]
		public static VectorGaussian ConcatAverageLogarithm(VectorGaussian first, VectorGaussian second, VectorGaussian result)
		{
			return ConcatAverageConditional(first, second, result);
		}

		/// <summary>
		/// VMP message to 'concat'
		/// </summary>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[sum_(second) p(second) factor(concat,first,second)]</c>.
		/// </para></remarks>
		public static VectorGaussian ConcatAverageLogarithm(Vector first, VectorGaussian second, VectorGaussian result)
		{
			return ConcatAverageConditional(first, second, result);
		}

		/// <summary>
		/// VMP message to 'concat'
		/// </summary>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'concat' as the random arguments are varied.
		/// The formula is <c>proj[sum_(first) p(first) factor(concat,first,second)]</c>.
		/// </para></remarks>
		public static VectorGaussian ConcatAverageLogarithm(VectorGaussian first, Vector second, VectorGaussian result)
		{
			return ConcatAverageConditional(first, second, result);
		}

		/// <summary>
		/// VMP message to 'first'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="second">Constant value for 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'first' with 'concat' integrated out.
		/// The formula is <c>sum_concat p(concat) factor(concat,first,second)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian FirstAverageLogarithm([SkipIfUniform] VectorGaussian concat, Vector second, VectorGaussian result)
		{
			return FirstAverageConditional(concat, second, result);
		}

		/// <summary>
		/// VMP message to 'first'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="second">Incoming message from 'second'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'first'.
		/// Because the factor is deterministic, 'concat' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(second) p(second) log(sum_concat p(concat) factor(concat,first,second)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian FirstAverageLogarithm([SkipIfUniform] VectorGaussian concat, VectorGaussian second, VectorGaussian result)
		{
			// prec = concat.Precision[dim1,dim1]
			// meanTimesPrec = concat.MeanTimesPrecision[dim1] - concat.Precision[dim1,dim2]*second.Mean
			Vector mSecond = second.GetMean();
			return FirstAverageConditional(concat, mSecond, result);
		}

		/// <summary>
		/// VMP message to 'second'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'second' with 'concat' integrated out.
		/// The formula is <c>sum_concat p(concat) factor(concat,first,second)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian SecondAverageLogarithm([SkipIfUniform] VectorGaussian concat, Vector first, VectorGaussian result)
		{
			return SecondAverageConditional(concat, first, result);
		}

		/// <summary>
		/// VMP message to 'second'
		/// </summary>
		/// <param name="concat">Incoming message from 'concat'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'second'.
		/// Because the factor is deterministic, 'concat' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(first) p(first) log(sum_concat p(concat) factor(concat,first,second)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="concat"/> is not a proper distribution</exception>
		public static VectorGaussian SecondAverageLogarithm([SkipIfUniform] VectorGaussian concat, VectorGaussian first, VectorGaussian result)
		{
			Vector mFirst = first.GetMean();
			return SecondAverageConditional(concat, mFirst, result);
		}

		/// <summary>
		/// VMP message to 'first'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'first' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian FirstAverageLogarithm(Vector concat, VectorGaussian result)
		{
			return FirstAverageConditional(concat, result);
		}

		/// <summary>
		/// VMP message to 'second'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Constant value for 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'second' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SecondAverageLogarithm(Vector concat, Vector first, VectorGaussian result)
		{
			return SecondAverageConditional(concat, first, result);
		}

		/// <summary>
		/// VMP message to 'second'
		/// </summary>
		/// <param name="concat">Constant value for 'concat'.</param>
		/// <param name="first">Incoming message from 'first'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'second'.
		/// The formula is <c>exp(sum_(first) p(first) log(factor(concat,first,second)))</c>.
		/// </para></remarks>
		public static VectorGaussian SecondAverageLogarithm(Vector concat, VectorGaussian first, VectorGaussian result)
		{
			return SecondAverageConditional(concat, first, result);
		}
	}
}
