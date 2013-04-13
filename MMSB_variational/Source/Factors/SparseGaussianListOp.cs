// (C) Copyright 2011 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;
using MicrosoftResearch.Infer.Collections;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="SparseGaussianList.Sample(IList<double>, IList<double>)"/> given random arguments to the function.
	/// </summary>
	[FactorMethod(new string[] { "sample", "mean", "precision" }, typeof(SparseGaussianList), "Sample", typeof(IList<double>), typeof(IList<double>))]
	[Quality(QualityBand.Experimental)]
	public static class SparseGaussianListOp
	{
		[Skip]
		public static SparseGaussianList SampleAverageConditionalInit([IgnoreDependency] IList<double> mean)
		{
			return SparseGaussianList.FromSize(mean.Count);
		}
		[Skip]
		public static SparseGaussianList SampleAverageConditionalInit([IgnoreDependency] SparseGaussianList mean)
		{
			return (SparseGaussianList)mean.Clone();
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(mean) p(mean) factor(sample,mean,precision)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageConditional(IList<double> mean, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageConditional(m, p));
			return result;
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static SparseGaussianList MeanAverageConditional(IList<double> sample, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageConditional(s, p));
			return result;
		}

		/// <summary>
		/// EP message to 'precision'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing EP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static SparseGammaList PrecisionAverageConditional(IList<double> sample, IList<double> mean, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, (s, m) => GaussianOp.PrecisionAverageConditional(s, m));
			return result;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(mean) p(mean) factor(sample,mean,precision)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageConditional([SkipIfUniform] SparseGaussianList mean, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageConditional(m, p));
			return result;
		}

		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(sample) p(sample) factor(sample,mean,precision)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageConditional([SkipIfUniform] SparseGaussianList sample, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageConditional(s, p));
			return result;
		}

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(mean,precision) p(mean,precision) factor(sample,mean,precision)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageConditional(SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.SampleAverageConditional(s, m, p));
			return result;
		}
		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'sample' as the random arguments are varied.
		/// The formula is <c>proj[p(sample) sum_(precision) p(precision) factor(sample,mean,precision)]/p(sample)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageConditional(SparseGaussianList sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction<Gaussian, double, Gamma>(sample, mean, precision, (s, m, p) => GaussianOp.SampleAverageConditional(s, m, p));
			return result;
		}
		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(sample,precision) p(sample,precision) factor(sample,mean,precision)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageConditional([SkipIfUniform] SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.MeanAverageConditional(s, m, p));
			return result;
		}
		/// <summary>
		/// EP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'mean' as the random arguments are varied.
		/// The formula is <c>proj[p(mean) sum_(precision) p(precision) factor(sample,mean,precision)]/p(mean)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageConditional(IList<double> sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision,  SparseGaussianList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.MeanAverageConditional(s, m, p));
			return result;
		}

		/// <summary>
		/// EP message to 'precision'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'precision' as the random arguments are varied.
		/// The formula is <c>proj[p(precision) sum_(mean) p(mean) factor(sample,mean,precision)]/p(precision)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageConditional(IList<double> sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.PrecisionAverageConditional(s, m, p));
			return result;
		}
		/// <summary>
		/// EP message to 'precision'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'precision' as the random arguments are varied.
		/// The formula is <c>proj[p(precision) sum_(sample) p(sample) factor(sample,mean,precision)]/p(precision)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageConditional([SkipIfUniform] SparseGaussianList sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.PrecisionAverageConditional(s, m, p));
			return result;
		}

		/// <summary>
		/// EP message to 'precision'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'precision' as the random arguments are varied.
		/// The formula is <c>proj[p(precision) sum_(sample,mean) p(sample,mean) factor(sample,mean,precision)]/p(precision)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageConditional([SkipIfUniform] SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, precision, (s, m, p) => GaussianOp.PrecisionAverageConditional(s, m, p));
			return result;
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
		public static double LogAverageFactor(IList<double> sample, IList<double> mean, IList<double> precision)
		{
			Func<double, double, double, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,mean) p(sample,mean) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, IList<double> precision)
		{
			Func<Gaussian, Gaussian, double, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample) p(sample) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] SparseGaussianList sample, IList<double> mean, IList<double> precision)
		{
			Func<Gaussian, double, double, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean) p(mean) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double LogAverageFactor(IList<double> sample, [SkipIfUniform] SparseGaussianList mean, IList<double> precision)
		{
			Func<double, Gaussian, double, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(precision) p(precision) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogAverageFactor(IList<double> sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<double, double, Gamma, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean,precision) p(mean,precision) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogAverageFactor(IList<double> sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<double, Gaussian, Gamma, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,precision) p(sample,precision) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] SparseGaussianList sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<Gaussian, double, Gamma, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,mean,precision) p(sample,mean,precision) factor(sample,mean,precision))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<Gaussian, Gaussian, Gamma, double> f = GaussianOp.LogAverageFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
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
		public static double LogEvidenceRatio(IList<double> sample, IList<double> mean, IList<double> precision)
		{
			Func<double, double, double, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,mean) p(sample,mean) factor(sample,mean,precision) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(SparseGaussianList sample, SparseGaussianList mean, IList<double> precision)
		{
			Func<Gaussian, Gaussian, double, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
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
		public static double LogEvidenceRatio(SparseGaussianList sample, IList<double> mean, IList<double> precision)
		{
			Func<Gaussian, double, double, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean) p(mean) factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(IList<double> sample, SparseGaussianList mean, IList<double> precision)
		{
			Func<double, Gaussian, double, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(precision) p(precision) factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(IList<double> sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<double, double, Gamma, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(mean,precision) p(mean,precision) factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(IList<double> sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision)
		{
			Func<double, Gaussian, Gamma, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,mean,precision) p(sample,mean,precision) factor(sample,mean,precision) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] SparseGaussianList sample, [SkipIfUniform] SparseGaussianList mean, [SkipIfUniform] SparseGammaList precision, [Fresh] SparseGaussianList to_sample)
		{
			Func<Gaussian, Gaussian, Gamma, Gaussian, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision, to_sample).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_sample">Outgoing message to 'sample'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(sample,precision) p(sample,precision) factor(sample,mean,precision) / sum_sample p(sample) messageTo(sample))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] SparseGaussianList sample, IList<double> mean, [SkipIfUniform] SparseGammaList precision, [Fresh] SparseGaussianList to_sample)
		{
			Func<Gaussian, double, Gamma, Gaussian, double> f = GaussianOp.LogEvidenceRatio;
			return f.Map(sample, mean, precision, to_sample).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		//-- VMP ------------------------------------------------------------------------------------------------

		[Skip]
		public static SparseGaussianList SampleAverageLogarithmInit([IgnoreDependency] IList<double> mean)
		{
			return SparseGaussianList.FromSize(mean.Count);
		}
		[Skip]
		public static SparseGaussianList SampleAverageLogarithmInit([IgnoreDependency] SparseGaussianList mean)
		{
			return (SparseGaussianList)mean.Clone();
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static SparseGaussianList SampleAverageLogarithm(IList<double> mean, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageLogarithm(m, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageLogarithm([Proper] SparseGaussianList mean, IList<double> precision,  SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageLogarithm(m, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static SparseGaussianList MeanAverageLogarithm(IList<double> sample, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageLogarithm(s, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageLogarithm([Proper] SparseGaussianList sample, IList<double> precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageLogarithm(s, p));
			return result;
		}

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(mean,precision) p(mean,precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageLogarithm([Proper] SparseGaussianList mean, [Proper] SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageLogarithm(m, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(sample,precision) p(sample,precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageLogarithm([Proper] SparseGaussianList sample, [Proper]SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageLogarithm(s, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'sample'.
		/// The formula is <c>exp(sum_(precision) p(precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList SampleAverageLogarithm(IList<double> mean, [Proper] SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(mean, precision, (m, p) => GaussianOp.SampleAverageLogarithm(m, p));
			return result;
		}
		/// <summary>
		/// VMP message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'mean'.
		/// The formula is <c>exp(sum_(precision) p(precision) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static SparseGaussianList MeanAverageLogarithm(IList<double> sample, [Proper]SparseGammaList precision, SparseGaussianList result)
		{
			result.SetToFunction(sample, precision, (s, p) => GaussianOp.MeanAverageLogarithm(s, p));
			return result;
		}

		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static SparseGammaList PrecisionAverageLogarithm(IList<double> sample, IList<double> mean, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, (s, m) => GaussianOp.PrecisionAverageLogarithm(s, m));
			return result;
		}
		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(sample,mean) p(sample,mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageLogarithm([Proper]SparseGaussianList sample, [Proper]SparseGaussianList mean, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, (s, m) => GaussianOp.PrecisionAverageLogarithm(s, m));
			return result;
		}

		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing VMP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(sample) p(sample) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageLogarithm([Proper]SparseGaussianList sample, IList<double> mean, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, (s, m) => GaussianOp.PrecisionAverageLogarithm(s, m));
			return result;
		}
		/// <summary>
		/// VMP message to 'precision'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'precision'.
		/// The formula is <c>exp(sum_(mean) p(mean) log(factor(sample,mean,precision)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static SparseGammaList PrecisionAverageLogarithm(IList<double> sample, [Proper]SparseGaussianList mean, SparseGammaList result)
		{
			result.SetToFunction(sample, mean, (s, m) => GaussianOp.PrecisionAverageLogarithm(s, m));
			return result;
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,mean,precision) p(sample,mean,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseGaussianList sample, [Proper] SparseGaussianList mean, [Proper] SparseGammaList precision)
		{
			Func<Gaussian, Gaussian, Gamma, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(precision) p(precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(IList<double> sample, IList<double> mean, [Proper] SparseGammaList precision)
		{
			Func<double, double, Gamma, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
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
		public static double AverageLogFactor(IList<double> sample, IList<double> mean, IList<double> precision)
		{
			Func<double, double, double, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample) p(sample) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseGaussianList sample, IList<double> mean, IList<double> precision)
		{
			Func<Gaussian, double, double, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean) p(mean) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor(IList<double> sample, [Proper] SparseGaussianList mean, IList<double> precision)
		{
			Func<double, Gaussian, double, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(mean,precision) p(mean,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor(IList<double> sample, [Proper] SparseGaussianList mean, [Proper] SparseGammaList precision)
		{
			Func<double, Gaussian, Gamma, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Incoming message from 'precision'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,precision) p(sample,precision) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="precision"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseGaussianList sample, IList<double> mean, [Proper] SparseGammaList precision)
		{
			Func<Gaussian, double, Gamma, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}

		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="mean">Incoming message from 'mean'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Average of the factor's log-value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>sum_(sample,mean) p(sample,mean) log(factor(sample,mean,precision))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="sample"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="mean"/> is not a proper distribution</exception>
		public static double AverageLogFactor([Proper] SparseGaussianList sample, [Proper] SparseGaussianList mean, IList<double> precision)
		{
			Func<Gaussian, Gaussian, double, double> f = GaussianOp.AverageLogFactor;
			return f.Map(sample, mean, precision).EnumerableReduce(0.0, (res, elem) => res + elem, (res, elem, cnt) => res + elem * cnt);
		}
	}
}
