// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.VectorGaussian"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(VectorGaussian), "Sample", typeof(Vector), typeof(PositiveDefiniteMatrix))]
	[FactorMethod(new string[] { "sample", "mean", "precision" }, typeof(Factor), "VectorGaussian")]
	[Buffers("SampleMean", "SampleVariance", "MeanMean", "MeanVariance", "PrecisionMean", "PrecisionMeanLogDet")]
	[Quality(QualityBand.Stable)]
	public static class VectorGaussianOp
	{
		/// <summary>
		/// Initialise the buffer 'SampleVariance'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'.</param>
		/// <returns>Initial value of buffer 'SampleVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix SampleVarianceInit([IgnoreDependency] VectorGaussian Sample)
		{
			return new PositiveDefiniteMatrix(Sample.Dimension, Sample.Dimension);
		}
		/// <summary>
		/// Update the buffer 'SampleVariance'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix SampleVariance([Proper] VectorGaussian Sample, PositiveDefiniteMatrix result)
		{
			return Sample.GetVariance(result);
		}
		/// <summary>
		/// Initialise the buffer 'SampleMean'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'.</param>
		/// <returns>Initial value of buffer 'SampleMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector SampleMeanInit([IgnoreDependency] VectorGaussian Sample)
		{
			return Vector.Zero(Sample.Dimension);
		}
		/// <summary>
		/// Update the buffer 'SampleMean'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		public static Vector SampleMean([Proper] VectorGaussian Sample, [Fresh] PositiveDefiniteMatrix SampleVariance, Vector result)
		{
			return Sample.GetMean(result, SampleVariance);
		}

		/// <summary>
		/// Initialise the buffer 'MeanVariance'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'.</param>
		/// <returns>Initial value of buffer 'MeanVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix MeanVarianceInit([IgnoreDependency] VectorGaussian Mean)
		{
			return new PositiveDefiniteMatrix(Mean.Dimension, Mean.Dimension);
		}
		/// <summary>
		/// Update the buffer 'MeanVariance'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix MeanVariance([Proper] VectorGaussian Mean, PositiveDefiniteMatrix result)
		{
			return Mean.GetVariance(result);
		}
		/// <summary>
		/// Initialise the buffer 'MeanMean'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'.</param>
		/// <returns>Initial value of buffer 'MeanMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector MeanMeanInit([IgnoreDependency] VectorGaussian Mean)
		{
			return Vector.Zero(Mean.Dimension);
		}
		/// <summary>
		/// Update the buffer 'MeanMean'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		public static Vector MeanMean([Proper] VectorGaussian Mean, [Fresh] PositiveDefiniteMatrix MeanVariance, Vector result)
		{
			return Mean.GetMean(result, MeanVariance);
		}

		/// <summary>
		/// Initialise the buffer 'PrecisionMean'
		/// </summary>
		/// <param name="Precision">Incoming message from 'precision'.</param>
		/// <returns>Initial value of buffer 'PrecisionMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix PrecisionMeanInit([IgnoreDependency] Wishart Precision)
		{
			return new PositiveDefiniteMatrix(Precision.Dimension, Precision.Dimension);
		}
		/// <summary>
		/// Update the buffer 'PrecisionMean'
		/// </summary>
		/// <param name="Precision">Incoming message from 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static PositiveDefiniteMatrix PrecisionMean([Proper] Wishart Precision, PositiveDefiniteMatrix result)
		{
			return Precision.GetMean(result);
		}
		/// <summary>
		/// Update the buffer 'PrecisionMeanLogDet'
		/// </summary>
		/// <param name="Precision">Incoming message from 'precision'.</param>
		/// <returns>New value of buffer 'PrecisionMeanLogDet'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double PrecisionMeanLogDet([Proper] Wishart Precision)
		{
			return Precision.GetMeanLogDeterminant();
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector sample, Vector mean, PositiveDefiniteMatrix precision)
		{
			int Dimension = sample.Count;
			LowerTriangularMatrix precL = new LowerTriangularMatrix(Dimension, Dimension);
			Vector iLb = Vector.Zero(Dimension);
			Vector precisionTimesMean = precision * mean;
			return VectorGaussian.GetLogProb(sample, precisionTimesMean, precision, precL, iLb);
		}

		/// <summary>
		/// Gibbs message to 'sample'
		/// </summary>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SampleConditional(Vector Mean, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			result.SetMeanAndPrecision(Mean, Precision);
			return result;
		}

		/// <summary>
		/// Gibbs message to 'mean'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian MeanConditional(Vector Sample, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleConditional(Sample, Precision, result);
		}

		/// <summary>
		/// Gibbs message to 'precision'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static Wishart PrecisionConditional(Vector Sample, Vector Mean, Wishart result, Vector diff)
		{
			if (result == default(Wishart)) result = new Wishart(Sample.Count);
			diff.SetToDifference(Sample, Mean);
			const double SQRT_HALF = 0.70710678118654752440084436210485;
			diff.Scale(SQRT_HALF);
			result.Rate.SetToOuter(diff, diff);
			result.Shape = 0.5 * (result.Dimension + 2);
			return result;
		}
		/// <summary>
		/// Gibbs message to 'precision'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static Wishart PrecisionConditional(Vector Sample, Vector Mean, Wishart result)
		{
			Vector workspace = Vector.Zero(Sample.Count);
			return PrecisionConditional(Sample, Mean, result, workspace);
		}

		//-- EP -----------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			PositiveDefiniteMatrix Precision)
		{
			return VectorGaussian.GetLogProb(SampleMean, MeanMean, Precision.Inverse() + SampleVariance + MeanVariance);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector Sample,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			PositiveDefiniteMatrix Precision)
		{
			return VectorGaussian.GetLogProb(Sample, MeanMean, Precision.Inverse() + MeanVariance);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector sample, Vector mean, PositiveDefiniteMatrix precision)
		{
			return LogAverageFactor(sample, mean, precision);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector sample,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			PositiveDefiniteMatrix precision)
		{
			return LogAverageFactor(sample, MeanMean, MeanVariance, precision);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,precision) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(VectorGaussian sample, Vector mean, PositiveDefiniteMatrix precision)
		{
			return 0.0;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,mean) p(sample,mean) factor(sample,mean,precision) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		[Skip]
		public static double LogEvidenceRatio([SkipIfUniform] VectorGaussian sample, [SkipIfUniform] VectorGaussian mean, PositiveDefiniteMatrix precision)
		{
			return 0.0;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SampleAverageConditional(Vector Mean, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleConditional(Mean, Precision, result);
		}
		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(mean) p(mean) factor(sample,mean,precision)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		public static VectorGaussian SampleAverageConditional([SkipIfUniform] VectorGaussian Mean, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			if (Mean.IsPointMass) return SampleConditional(Mean.Point, Precision, result);
			if (result == default(VectorGaussian)) result = new VectorGaussian(Mean.Dimension);
			// R = Prec/(Prec + Mean.Prec)
			PositiveDefiniteMatrix R = Precision + Mean.Precision;
			R.SetToProduct(Precision, R.Inverse());
			for (int i = 0; i < Mean.Dimension; i++) {
				if (double.IsPositiveInfinity(Mean.Precision[i, i])) R[i, i] = 1;
			}
			result.Precision.SetToProduct(R, Mean.Precision);
			for (int i = 0; i < Mean.Dimension; i++) {
				if (double.IsPositiveInfinity(Mean.Precision[i, i])) {
					for (int j = 0; j < Mean.Dimension; j++) {
						result.Precision[i, j] = 0;
						result.Precision[j, i] = 0;
					}
					result.Precision[i, i] = 1;
				}
			}
			result.MeanTimesPrecision.SetToProduct(R, Mean.MeanTimesPrecision);
			return result;
		}
		[Skip]
		public static VectorGaussian SampleAverageConditionalInit(Vector Mean)
		{
			return VectorGaussian.Uniform(Mean.Count);
		}
		[Skip]
		public static VectorGaussian SampleAverageConditionalInit([IgnoreDependency] VectorGaussian Mean)
		{
			return new VectorGaussian(Mean.Dimension);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian MeanAverageConditional(Vector Sample, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleConditional(Sample, Precision, result);
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(sample) p(sample) factor(sample,mean,precision)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		public static VectorGaussian MeanAverageConditional([SkipIfUniform] VectorGaussian Sample, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleAverageConditional(Sample, Precision, result);
		}

		/// <summary>
		/// EP message to 'precision'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static Wishart PrecisionAverageConditional(Vector Sample, Vector Mean, Wishart result)
		{
			return PrecisionConditional(Sample, Mean, result);
		}

		//-- VMP ----------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="precisionMean">Buffer 'precisionMean'.</param>
		/// <param name="precisionMeanLogDet">Buffer 'precisionMeanLogDet'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,mean,precision) p(sample,mean,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(
			[Proper] VectorGaussian sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			[Proper] VectorGaussian mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			[Proper] Wishart precision,
			[Fresh] PositiveDefiniteMatrix precisionMean,
			[Fresh] double precisionMeanLogDet)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, MeanMean, MeanVariance, precision, precisionMean, precisionMeanLogDet);
			if (mean.IsPointMass)
				return AverageLogFactor(sample, SampleMean, SampleVariance, mean.Point, precision, precisionMean, precisionMeanLogDet);
			if (precision.IsPointMass)
				return AverageLogFactor(sample, SampleMean, SampleVariance, mean, MeanMean, MeanVariance, precision.Point);

			return ComputeAverageLogFactor(SampleMean, SampleVariance, MeanMean, MeanVariance, precisionMeanLogDet, precisionMean);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="precisionMean">Buffer 'precisionMean'.</param>
		/// <param name="precisionMeanLogDet">Buffer 'precisionMeanLogDet'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(precision) p(precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector sample, Vector mean,
			[Proper] Wishart precision,
			[Fresh] PositiveDefiniteMatrix precisionMean,
			[Fresh] double precisionMeanLogDet)
		{
			if (precision.IsPointMass)
				return AverageLogFactor(sample, mean, precision.Point);
			else
				return ComputeAverageLogFactor(sample, mean, precisionMeanLogDet, precisionMean);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Vector sample, Vector mean, PositiveDefiniteMatrix precision)
		{
			return ComputeAverageLogFactor(sample, mean, precision.LogDeterminant(ignoreInfinity: true), precision);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double AverageLogFactor(
			[Proper] VectorGaussian sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			Vector mean, PositiveDefiniteMatrix precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision);
			else
				return ComputeAverageLogFactor(SampleMean, SampleVariance, mean, precision.LogDeterminant(ignoreInfinity: true), precision);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean) p(mean) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector sample,
			[Proper] VectorGaussian mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			PositiveDefiniteMatrix precision)
		{
			return AverageLogFactor(mean, MeanMean, MeanVariance, sample, precision);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="precisionMean">Buffer 'precisionMean'.</param>
		/// <param name="precisionMeanLogDet">Buffer 'precisionMeanLogDet'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean,precision) p(mean,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(Vector sample,
			[Proper] VectorGaussian mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			[Proper] Wishart precision,
			[Fresh] PositiveDefiniteMatrix precisionMean,
			[Fresh] double precisionMeanLogDet)
		{
			return AverageLogFactor(mean, MeanMean, MeanVariance, sample, precision, precisionMean, precisionMeanLogDet);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="precisionMean">Buffer 'precisionMean'.</param>
		/// <param name="precisionMeanLogDet">Buffer 'precisionMeanLogDet'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,precision) p(sample,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(
			[Proper] VectorGaussian sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			Vector mean,
			[Proper] Wishart precision,
			[Fresh] PositiveDefiniteMatrix precisionMean,
			[Fresh] double precisionMeanLogDet)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision, precisionMean, precisionMeanLogDet);
			if (precision.IsPointMass)
				return AverageLogFactor(sample, SampleMean, SampleVariance, mean, precision.Point);

			return ComputeAverageLogFactor(SampleMean, SampleVariance, mean, precisionMeanLogDet, precisionMean);
		}
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,mean) p(sample,mean) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(
			[Proper] VectorGaussian sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			[Proper] VectorGaussian mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			PositiveDefiniteMatrix precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, MeanMean, MeanVariance, precision);
			if (mean.IsPointMass)
				return AverageLogFactor(sample, SampleMean, SampleVariance, mean.Point, precision);

			return ComputeAverageLogFactor(SampleMean, SampleVariance, MeanMean, MeanVariance, precision.LogDeterminant(ignoreInfinity: true), precision);
		}
		/// <summary>
		/// Helper method for computing average log factor
		/// </summary>
		/// <param name="SampleMean">Mean of incoming message from 'sample'</param>
		/// <param name="SampleVariance">Variance of incoming message from 'sample'</param>
		/// <param name="MeanMean">Mean of incoming message from 'mean'</param>
		/// <param name="MeanVariance">Variance of incoming message from 'mean'</param>
		/// <param name="precision_Elogx">Expected log value of the incoming message from 'precision'</param>
		/// <param name="precision_Ex">Expected value of incoming message from 'precision'</param>
		/// <returns>Computed average log factor</returns>
		private static double ComputeAverageLogFactor(
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			double precision_Elogx, PositiveDefiniteMatrix precision_Ex)
		{
			int dim = SampleMean.Count;
			int nonzeroDims = 0;
			double precTimesVariance = 0.0;
			double precTimesDiff = 0.0;
			for (int i = 0; i < dim; i++) {
				if (double.IsPositiveInfinity(precision_Ex[i, i])) {
					if (SampleMean[i] != MeanMean[i] || SampleVariance[i, i]+MeanVariance[i,i] > 0) return double.NegativeInfinity;
				} else {
					nonzeroDims++;
					double sum = 0.0;
					for (int j = 0; j < dim; j++) {
						sum += precision_Ex[i, j]*(SampleMean[j] - MeanMean[j]);
						precTimesVariance += precision_Ex[i, j]*(SampleVariance[i, j] + MeanVariance[i,j]);
					}
					precTimesDiff += sum*(SampleMean[i] - MeanMean[i]);
				}
			}
			return -nonzeroDims * MMath.LnSqrt2PI + 0.5 * (precision_Elogx - precTimesVariance - precTimesDiff);
		}
		/// <summary>
		/// Helper method for computing average log factor
		/// </summary>
		/// <param name="SampleMean">Mean of incoming sample message</param>
		/// <param name="SampleVariance">Variance of incoming sample message</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision_Elogx">Expected log value of the incoming message from 'precision'</param>
		/// <param name="precision_Ex">Expected value of incoming message from 'precision'</param>
		/// <returns>Computed average log factor</returns>
		private static double ComputeAverageLogFactor(
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			Vector mean,
			double precision_Elogx, PositiveDefiniteMatrix precision_Ex)
		{
			int dim = mean.Count;
			int nonzeroDims = 0;
			double precTimesVariance = 0.0;
			double precTimesDiff = 0.0;
			for (int i = 0; i < dim; i++) {
				if (double.IsPositiveInfinity(precision_Ex[i, i])) {
					if (SampleMean[i] != mean[i] || SampleVariance[i,i] > 0) return double.NegativeInfinity;
				} else {
					nonzeroDims++;
					double sum = 0.0;
					for (int j = 0; j < dim; j++) {
						sum += precision_Ex[i, j]*(SampleMean[j] - mean[j]);
						precTimesVariance += precision_Ex[i, j]*SampleVariance[j, i];
					}
					precTimesDiff += sum*(SampleMean[i]-mean[i]);
				}
			}
			return -nonzeroDims * MMath.LnSqrt2PI + 0.5 * (precision_Elogx - precTimesVariance - precTimesDiff);

		}
		/// <summary>
		/// Helper method for computing average log factor
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision_Elogx">Expected log value of the incoming message from 'precision'</param>
		/// <param name="precision_Ex">Expected value of incoming message from 'precision'</param>
		/// <returns>Computed average log factor</returns>
		private static double ComputeAverageLogFactor(Vector sample, Vector mean, double precision_Elogx, PositiveDefiniteMatrix precision_Ex)
		{
			int dim = mean.Count;
			int nonzeroDims = 0;
			double precTimesDiff = 0.0;
			for (int i = 0; i < dim; i++) {
				if (double.IsPositiveInfinity(precision_Ex[i, i])) {
					if (sample[i] != mean[i]) return double.NegativeInfinity;
				} else {
					nonzeroDims++;
					double sum = 0.0;
					for (int j = 0; j < dim; j++) {
						sum += precision_Ex[i, j]*(sample[j] - mean[j]);
					}
					precTimesDiff += sum*(sample[i]-mean[i]);
				}
			}
			return -nonzeroDims * MMath.LnSqrt2PI + 0.5 * (precision_Elogx - precTimesDiff);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SampleAverageLogarithm(Vector Mean, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleConditional(Mean, Precision, result);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="Precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="PrecisionMean">Buffer 'PrecisionMean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(mean,precision) p(mean,precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="Precision"/> is not a proper distribution</exception>
		public static VectorGaussian SampleAverageLogarithm([Proper] VectorGaussian Mean, [Fresh] Vector MeanMean, [Proper] Wishart Precision, [Fresh] PositiveDefiniteMatrix PrecisionMean, VectorGaussian result)
		{
			return SampleAverageLogarithm(MeanMean, PrecisionMean, result);
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static VectorGaussian SampleAverageLogarithm([Proper] VectorGaussian mean, [Fresh] Vector MeanMean, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			if (result == default(VectorGaussian)) result = new VectorGaussian(MeanMean.Count);
			result.Precision.SetTo(Precision);
			result.MeanTimesPrecision.SetToProduct(result.Precision, MeanMean);
			return result;
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="Precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="PrecisionMean">Buffer 'PrecisionMean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(precision) p(precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Precision"/> is not a proper distribution</exception>
		public static VectorGaussian SampleAverageLogarithm(Vector Mean, [Proper] Wishart Precision, [Fresh] PositiveDefiniteMatrix PrecisionMean, VectorGaussian result)
		{
			return SampleAverageLogarithm(Mean, PrecisionMean, result);
		}

		[Skip]
		public static VectorGaussian SampleAverageLogarithmInit(Vector Mean)
		{
			return VectorGaussian.Uniform(Mean.Count);
		}
		[Skip]
		public static VectorGaussian SampleAverageLogarithmInit([IgnoreDependency] VectorGaussian Mean)
		{
			return new VectorGaussian(Mean.Dimension);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian MeanAverageLogarithm(Vector Sample, PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleConditional(Sample, Precision, result);
		}

		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="Precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="PrecisionMean">Buffer 'PrecisionMean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(sample,precision) p(sample,precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="Precision"/> is not a proper distribution</exception>
		public static VectorGaussian MeanAverageLogarithm([Proper] VectorGaussian Sample, [Fresh] Vector SampleMean, [Proper] Wishart Precision, [Fresh] PositiveDefiniteMatrix PrecisionMean, VectorGaussian result)
		{
			return SampleAverageLogarithm(Sample, SampleMean, Precision, PrecisionMean, result);
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="Precision">Constant value for 'precision'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		public static VectorGaussian MeanAverageLogarithm([Proper] VectorGaussian Sample,
			[Fresh] Vector SampleMean,
			PositiveDefiniteMatrix Precision, VectorGaussian result)
		{
			return SampleAverageLogarithm(Sample, SampleMean, Precision, result);
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Precision">Incoming message from 'precision'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="PrecisionMean">Buffer 'PrecisionMean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(precision) p(precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Precision"/> is not a proper distribution</exception>
		public static VectorGaussian MeanAverageLogarithm(Vector Sample, [Proper] Wishart Precision, [Fresh] PositiveDefiniteMatrix PrecisionMean, VectorGaussian result)
		{
			return SampleAverageLogarithm(Sample, Precision, PrecisionMean, result);
		}

		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static Wishart PrecisionAverageLogarithm(Vector Sample, Vector Mean, Wishart result)
		{
			return PrecisionConditional(Sample, Mean, result);
		}
		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(sample,mean) p(sample,mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		public static Wishart PrecisionAverageLogarithm(
			[Proper] VectorGaussian Sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			[Proper] VectorGaussian Mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			Wishart result)
		{
			if (Sample.IsPointMass) return PrecisionAverageLogarithm(Sample.Point, Mean, MeanMean, MeanVariance, result);
			if (Mean.IsPointMass) return PrecisionAverageLogarithm(Sample, SampleMean, SampleVariance, Mean.Point, result);
			// The formula is exp(int_x int_mean p(x) p(mean) log N(x;mean,1/prec)) =
			// exp(-0.5 prec E[(x-mean)^2] + 0.5 log(prec)) =
			// Gamma(prec; 0.5, 0.5*E[(x-mean)^2])
			// E[(x-mean)^2] = E[x^2] - 2 E[x] E[mean] + E[mean^2] = var(x) + (E[x]-E[mean])^2 + var(mean)
			if (result == default(Wishart)) result = new Wishart(Sample.Dimension);
			result.Shape = 0.5 * (result.Dimension + 2);
			Vector diff = SampleMean - MeanMean;
			result.Rate.SetToOuter(diff, diff);
			result.Rate.SetToSum(result.Rate, SampleVariance);
			result.Rate.SetToSum(result.Rate, MeanVariance);
			result.Rate.Scale(0.5);
			return result;
		}
		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="Sample">Constant value for 'sample'.</param>
		/// <param name="Mean">Incoming message from 'mean'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanMean">Buffer 'MeanMean'.</param>
		/// <param name="MeanVariance">Buffer 'MeanVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Mean"/> is not a proper distribution</exception>
		public static Wishart PrecisionAverageLogarithm(Vector Sample,
			[Proper] VectorGaussian Mean,
			[Fresh] Vector MeanMean,
			[Fresh] PositiveDefiniteMatrix MeanVariance,
			Wishart result)
		{
			if (Mean.IsPointMass) return PrecisionAverageLogarithm(Sample, Mean.Point, result);
			// The formula is exp(int_x int_mean p(x) p(mean) log N(x;mean,1/prec)) =
			// exp(-0.5 prec E[(x-mean)^2] + 0.5 log(prec)) =
			// Gamma(prec; 0.5, 0.5*E[(x-mean)^2])
			// E[(x-mean)^2] = E[x^2] - 2 E[x] E[mean] + E[mean^2] = var(x) + (E[x]-E[mean])^2 + var(mean)
			if (result == default(Wishart)) result = new Wishart(Sample.Count);
			result.Shape = 0.5 * (result.Dimension + 2);
			Vector diff = Sample - MeanMean;
			result.Rate.SetToOuter(diff, diff);
			result.Rate.SetToSum(result.Rate, MeanVariance);
			result.Rate.Scale(0.5);
			return result;
		}
		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="Sample">Incoming message from 'sample'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="SampleMean">Buffer 'SampleMean'.</param>
		/// <param name="SampleVariance">Buffer 'SampleVariance'.</param>
		/// <param name="Mean">Constant value for 'mean'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Sample"/> is not a proper distribution</exception>
		public static Wishart PrecisionAverageLogarithm(
			[Proper] VectorGaussian Sample,
			[Fresh] Vector SampleMean,
			[Fresh] PositiveDefiniteMatrix SampleVariance,
			Vector Mean, Wishart result)
		{
			return PrecisionAverageLogarithm(Mean, Sample, SampleMean, SampleVariance, result);
		}
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="VectorGaussian.SampleFromMeanAndVariance"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(VectorGaussian), "SampleFromMeanAndVariance")]
	[Quality(QualityBand.Stable)]
	public static class VectorGaussianFromMeanAndVarianceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Vector sample, Vector mean, PositiveDefiniteMatrix variance)
		{
			VectorGaussian to_sample = SampleAverageConditional(mean, variance);
			return to_sample.GetLogProb(sample);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(Vector sample, Vector mean, PositiveDefiniteMatrix variance) { return LogAverageFactor(sample, mean, variance); }
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(sample,mean,variance))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		public static double AverageLogFactor(Vector sample, Vector mean, PositiveDefiniteMatrix variance) { return LogAverageFactor(sample, mean, variance); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,variance))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(VectorGaussian sample, [Fresh] VectorGaussian to_sample)
		{
			return to_sample.GetLogAverageOf(sample);
		}

		public static double LogEvidenceRatio(Vector sample, VectorGaussian mean, PositiveDefiniteMatrix variance)
		{
			return SampleAverageConditional(sample, variance).GetLogAverageOf(mean);
		}
		[Skip]
		public static double LogEvidenceRatio(VectorGaussian sample, Vector mean, PositiveDefiniteMatrix variance) { return 0.0; }

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SampleAverageLogarithm(Vector mean, PositiveDefiniteMatrix variance)
		{
			return VectorGaussian.FromMeanAndVariance(mean, variance);
		}
		public static VectorGaussian MeanAverageLogarithm(Vector sample, PositiveDefiniteMatrix variance)
		{
			return SampleAverageLogarithm(sample, variance);
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="variance">Constant value for 'variance'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static VectorGaussian SampleAverageConditional(Vector mean, PositiveDefiniteMatrix variance)
		{
			return VectorGaussian.FromMeanAndVariance(mean, variance);
		}
		public static VectorGaussian MeanAverageConditional(Vector sample, PositiveDefiniteMatrix variance)
		{
			return SampleAverageConditional(sample, variance);
		}
	}
}
