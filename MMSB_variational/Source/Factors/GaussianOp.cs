// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using MicrosoftResearch.Infer.Utils;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Gaussian"/> and <see cref="Gaussian.Sample(double,double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Gaussian), "Sample", typeof(double), typeof(double))]
	[FactorMethod(new string[] { "sample", "mean", "precision" }, typeof(Factor), "Gaussian")]
	[Quality(QualityBand.Mature)]
	public static class GaussianOp
	{
		public static TruncatedGaussian SampleAverageConditional(double mean, double precision, TruncatedGaussian result)
		{
			return TruncatedGaussian.FromGaussian(Gaussian.FromMeanAndPrecision(mean, precision));
		}
		public static TruncatedGaussian MeanAverageConditional(double sample, double precision, TruncatedGaussian result)
		{
			return SampleAverageConditional(sample, precision, result);
		}

		//-- Gibbs ---------------------------------------------------------------------------------------------

		/// <summary>
		/// Evidence message for Gibbs
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>Logarithm of the factor's value at the given arguments</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		public static double LogFactorValue(double sample, double mean, double precision)
		{
			return Gaussian.GetLogProb(sample, mean, 1.0 / precision);
		}

		/// <summary>
		/// Gibbs message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing Gibbs message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian SampleConditional(double mean, double precision)
		{
			return Gaussian.FromMeanAndPrecision(mean, precision);
		}
		/// <summary>
		/// Gibbs message to 'mean'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing Gibbs message to the 'mean' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'mean' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian MeanConditional(double sample, double precision)
		{
			return SampleConditional(sample, precision);
		}
		/// <summary>
		/// Gibbs message to 'precision'
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <returns>The outgoing Gibbs message to the 'precision' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'precision' conditioned on the given values.
		/// </para></remarks>
		public static Gamma PrecisionConditional(double sample, double mean)
		{
			Gamma result = new Gamma();
			double diff = sample - mean;
			result.Rate = 0.5 * diff * diff;
			result.Shape = 1.5;
			return result;
		}

		//-- EP ------------------------------------------------------------------------------------------------
		/// <summary>
		/// Static flag to force a proper distribution
		/// </summary>
		public static bool ForceProper;

		/// <summary>
		/// EP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing EP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian SampleAverageConditional(double mean, double precision)
		{
			return SampleConditional(mean, precision);
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
		public static Gaussian MeanAverageConditional(double sample, double precision)
		{
			return MeanConditional(sample, precision);
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
		public static Gamma PrecisionAverageConditional(double sample, double mean)
		{
			return PrecisionConditional(sample, mean);
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
		public static Gaussian SampleAverageConditional([SkipIfUniform] Gaussian mean, double precision)
		{
			if (mean.IsPointMass) return SampleConditional(mean.Point, precision);
			// if (precision < 0) throw new ArgumentException("The constant precision given to the Gaussian factor is negative", "precision");
			if (precision == 0) {
				return Gaussian.Uniform();
			} else if (double.IsPositiveInfinity(precision)) {
				return mean;
			} else {
				if (mean.Precision <= -precision) throw new ImproperMessageException(mean);
				// The formula is int_mean N(x;mean,1/prec) p(mean) = N(x; mm, mv + 1/prec)
				// sample.Precision = inv(mv + inv(prec)) = mprec*prec/(prec + mprec)
				// sample.MeanTimesPrecision = sample.Precision*mm = R*(mprec*mm)
				// R = Prec/(Prec + mean.Prec)
				// This code works for mean.IsUniform() since then mean.Precision = 0, mean.MeanTimesPrecision = 0
				Gaussian result = new Gaussian();
				double R = precision / (precision + mean.Precision);
				result.Precision = R * mean.Precision;
				result.MeanTimesPrecision = R * mean.MeanTimesPrecision;
				return result;
			}

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
		public static Gaussian MeanAverageConditional([SkipIfUniform] Gaussian sample, double precision)
		{
			return SampleAverageConditional(sample, precision);
		}

		/// <summary>
		/// Number of quadrature nodes to use for computing the messages.
		/// Reduce this number to save time in exchange for less accuracy.
		/// </summary>
		public static int QuadratureNodeCount = 50;

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
		public static Gaussian SampleAverageConditional(Gaussian sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			Gaussian result = new Gaussian();
			if (precision.IsPointMass) {
				return SampleAverageConditional(mean, precision.Point);
			} else if (sample.IsUniform()) {
				// for large vx, Z =approx N(mx; mm, vx+vm+E[1/prec])
				double mm,mv;
				mean.GetMeanAndVariance(out mm, out mv);
				// NOTE: this error may happen because sample didn't receive any message yet under the schedule.
				// Need to make the scheduler smarter to avoid this.
				if(precision.Shape <= 1.0) throw new ArgumentException("The posterior has infinite variance due to precision distributed as "+precision+" (shape <= 1).  Try using a different prior for the precision, with shape > 1.");
				return Gaussian.FromMeanAndVariance(mm, mv + precision.GetMeanInverse());
			} else if (mean.IsUniform() || precision.IsUniform()) {
				result.SetToUniform();
			} else if (sample.IsPointMass) {
				// The correct answer here is not uniform, but rather a limit.  
				// However it doesn't really matter what we return since multiplication by a point mass 
				// always yields a point mass.
				result.SetToUniform();
			} else if (!precision.IsProper()) {
				throw new ImproperMessageException(precision);
			} else {
				// The formula is int_prec int_mean N(x;mean,1/prec) p(x) p(mean) p(prec) =
				// int_prec N(x; mm, mv + 1/prec) p(x) p(prec) =
				// int_prec N(x; new xm, new xv) N(xm; mm, mv + xv + 1/prec) p(prec)
				// Let R = Prec/(Prec + mean.Prec)
				// new xv = inv(R*mean.Prec + sample.Prec)
				// new xm = xv*(R*mean.PM + sample.PM)

				// In the case where sample and mean are improper distributions, 
				// we must only consider values of prec for which (new xv > 0).
				// This happens when R*mean.Prec > -sample.Prec
				// As a function of Prec, R*mean.Prec has a singularity at Prec=-mean.Prec
				// This function is greater than a threshold when Prec is sufficiently small or sufficiently large.
				// Therefore we construct an interval of Precs to exclude from the integration.
				double xm, xv, mm, mv;
				sample.GetMeanAndVarianceImproper(out xm, out xv);
				mean.GetMeanAndVarianceImproper(out mm, out mv);
				double lowerBound = 0;
				double upperBound = Double.PositiveInfinity;
				bool precisionIsBetween = true;
				if (mean.Precision >= 0) {
					if (sample.Precision < -mean.Precision) throw new ImproperMessageException(sample);
					//lowerBound = -mean.Precision * sample.Precision / (mean.Precision + sample.Precision);
					lowerBound = -1.0 / (xv + mv);
				} else {  // mean.Precision < 0
					if (sample.Precision < 0) {
						precisionIsBetween = true;
						lowerBound = -1.0 / (xv + mv);
						upperBound = -mean.Precision;
					} else if (sample.Precision < -mean.Precision) {
						precisionIsBetween = true;
						lowerBound = 0;
						upperBound = -mean.Precision;
					} else {
						// in this case, the precision should NOT be in this interval.
						precisionIsBetween = false;
						lowerBound = -mean.Precision;
						lowerBound = -1.0 / (xv + mv);
					}
				}
				double[] nodes = new double[QuadratureNodeCount];
				double[] weights = new double[nodes.Length];
				QuadratureNodesAndWeights(precision, nodes, weights);
				double Z = 0, rmean = 0, rvariance = 0;
				for (int i = 0; i < nodes.Length; i++) {
					double newVar, newMean;
					Assert.IsTrue(nodes[i] > 0);
					if ((nodes[i] > lowerBound && nodes[i] < upperBound) != precisionIsBetween) continue;
					// the following works even if sample is uniform. (sample.Precision == 0)
					if (mean.IsPointMass) {
						// take limit mean.Precision -> Inf
						newVar = 1.0 / (nodes[i] + sample.Precision);
						newMean = newVar * (nodes[i] * mean.Point + sample.MeanTimesPrecision);
					} else {
						// mean.Precision < Inf
						double R = nodes[i] / (nodes[i] + mean.Precision);
						newVar = 1.0 / (R * mean.Precision + sample.Precision);
						newMean = newVar * (R * mean.MeanTimesPrecision + sample.MeanTimesPrecision);
					}

					double f;
					// If p(x) is uniform, xv=Inf and the term N(xm; mm, mv + xv + 1/prec) goes away
					if (sample.IsUniform())
						f = weights[i];
					else
						f = weights[i] * Math.Exp(Gaussian.GetLogProb(xm, mm, xv + mv + 1.0 / nodes[i]));
					double fm = f * newMean;
					double fmm = f * (newVar + newMean * newMean);
					Z += f;
					rmean += fm;
					rvariance += fmm;
				}
				if (Z == 0.0) {
					//throw new Exception("Quadrature failed");
					//Console.WriteLine("Warning: Quadrature found zero mass.  Results may be inaccurate.");
					result.SetToUniform();
					return result;
				}
				double s = 1.0 / Z;
				rmean *= s;
				rvariance = rvariance * s - rmean * rmean;
				if (Double.IsInfinity(rmean)) {
					result.SetToUniform();
				} else {
					result.SetMeanAndVariance(rmean, rvariance);
					if (ForceProper) result.SetToRatioProper(result, sample);
					else result.SetToRatio(result, sample);
				}
			}
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
		public static Gaussian SampleAverageConditional(Gaussian sample, double mean, [SkipIfUniform] Gamma precision)
		{
			return SampleAverageConditional(sample, Gaussian.PointMass(mean), precision);
		}
		[Skip]
		public static Gaussian SampleAverageConditionalInit()
		{
			return Gaussian.Uniform();
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
		public static Gaussian MeanAverageConditional([SkipIfUniform] Gaussian sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			// TM: added SkipIfUniform to mean to encourage good schedules.
			// The result depends on mean, but can be non-uniform even if mean is uniform.
			return SampleAverageConditional(mean, sample, precision);
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
		public static Gaussian MeanAverageConditional(double sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			// TM: added SkipIfUniform to mean to encourage good schedules.
			// The result depends on mean, but can be non-uniform even if mean is uniform.
			return SampleAverageConditional(mean, sample, precision);
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
		public static Gamma PrecisionAverageConditional(double sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			return PrecisionAverageConditional(Gaussian.PointMass(sample), mean, precision);
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
		public static Gamma PrecisionAverageConditional([SkipIfUniform] Gaussian sample, double mean, [SkipIfUniform] Gamma precision)
		{
			return PrecisionAverageConditional(sample, Gaussian.PointMass(mean), precision);
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
		public static Gamma PrecisionAverageConditional([SkipIfUniform] Gaussian sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			if (sample.IsPointMass && mean.IsPointMass)
				return PrecisionAverageConditional(sample.Point, mean.Point);

			Gamma result = new Gamma();
			if (precision.IsPointMass) {
				// The correct answer here is not uniform, but rather a limit.  
				// However it doesn't really matter what we return since multiplication by a point mass 
				// always yields a point mass.
				result.SetToUniform();
			} else if (sample.IsUniform() || mean.IsUniform()) {
				result.SetToUniform();
			} else if (!precision.IsProper()) {
				// improper prior
				throw new ImproperMessageException(precision);
			} else {
				// use quadrature to integrate over the precision
				// see LogAverageFactor
				double xm, xv, mm, mv;
				sample.GetMeanAndVarianceImproper(out xm, out xv);
				mean.GetMeanAndVarianceImproper(out mm, out mv);
				double upperBound = Double.PositiveInfinity;
				if (xv + mv < 0) upperBound = -1.0 / (xv + mv);
				double[] nodes = new double[QuadratureNodeCount];
				double[] weights = new double[nodes.Length];
				QuadratureNodesAndWeights(precision, nodes, weights);
				double Z = 0, rmean = 0, rvariance = 0;
				double shift = 0;
				for (int i = 0; i < nodes.Length; i++) {
					double v = 1.0 / nodes[i] + xv + mv;
					if (v < 0) continue;
					double lp = Gaussian.GetLogProb(xm, mm, v);
					if (shift == 0) shift = lp;
					double f = weights[i] * Math.Exp(lp - shift);
					double fm = f * nodes[i];
					double fmm = fm * nodes[i];
					Z += f;
					rmean += fm;
					rvariance += fmm;
				}

                // Adaptive Clenshaw-Curtis quadrature: gives same results on easy integrals but 
                // still fails ExpFactorTest2
                //Converter<double,double> func = delegate(double y) {
                //    double x = Math.Exp(y); 
                //    double v = 1.0 / x + xv + mv;
                //    if (v < 0) return 0.0;
                //    return Math.Exp(Gaussian.GetLogProb(xm, mm, v) + Gamma.GetLogProb(x, precision.Shape+1, precision.Rate));
                //};
                //double Z2 = BernoulliFromLogOddsOp.AdaptiveClenshawCurtisQuadrature(func, 1, 16, 1e-10);
                //Converter<double, double> func2 = delegate(double y)
                //{
                //    return Math.Exp(y) * func(y); 
                //};
                //double rmean2 = BernoulliFromLogOddsOp.AdaptiveClenshawCurtisQuadrature(func2, 1, 16, 1e-10);
                //Converter<double, double> func3 = delegate(double y)
                //{
                //    return Math.Exp(2*y) * func(y);
                //};
                //double rvariance2 = BernoulliFromLogOddsOp.AdaptiveClenshawCurtisQuadrature(func3, 1, 16, 1e-10);
                //rmean2 = rmean2/ Z2;
                //rvariance2 = rvariance2 / Z2 - rmean2 * rmean2;

				if (Z == 0.0) {
					//throw new Exception("Quadrature failed");
					//Console.WriteLine("Warning: Quadrature found zero mass.  Results may be inaccurate.");
					result.SetToUniform();
					return result;
				}
				double s = 1.0 / Z;
				rmean *= s;
				rvariance = rvariance * s - rmean * rmean;
				if (Double.IsInfinity(rmean)) {
					result.SetToUniform();
				} else {
					result.SetMeanAndVariance(rmean, rvariance);
					if (ForceProper) result.SetToRatioProper(result, precision);
					else result.SetToRatio(result, precision);
				}
			}
#if KeepLastMessage
			if (LastPrecisionMessage != null) {
				if (Stepsize != 1 && Stepsize != 0) {
					LastPrecisionMessage.SetToPower(LastPrecisionMessage, 1 - Stepsize);
					result.SetToPower(result, Stepsize);
					result.SetToProduct(result, LastPrecisionMessage);
				}
			}
			// FIXME: this is not entirely safe since the caller could overwrite result.
			LastPrecisionMessage = result;
#endif
			return result;
		}

		/// <summary>
		/// Quadrature nodes for Gamma expectations
		/// </summary>
		/// <param name="precision">'precision' message</param>
		/// <param name="nodes">Place to put the nodes</param>
		/// <param name="weights">Place to put the weights</param>
		public static void QuadratureNodesAndWeights(Gamma precision, double[] nodes, double[] weights)
		{
#if KeepLastMessage
			if (LastPrecisionMessage != null) {
				Gamma PrecisionPosterior = precision * LastPrecisionMessage;
				Quadrature.GammaNodesAndWeights(PrecisionPosterior.Precision, PrecisionPosterior.PrecisionOverMean, nodes, weights);
				// modify the weights to include q(prec)/Ga(prec;a,b)
				for (int i = 0; i < weights.Length; i++) {
					weights[i] *= Math.Exp(precision.EvaluateLn(nodes[i]) - Gamma.EvaluateLn(nodes[i], PrecisionPosterior.Precision, PrecisionPosterior.PrecisionOverMean));
				}
				return;
			}
#endif
			Quadrature.GammaNodesAndWeights(precision.Shape - 1, precision.Rate, nodes, weights);
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
		public static double LogAverageFactor(double sample, double mean, double precision)
		{
			return Gaussian.GetLogProb(sample, mean, 1.0 / precision);
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
		public static double LogAverageFactor([SkipIfUniform] Gaussian sample, [SkipIfUniform] Gaussian mean, double precision)
		{
			return GaussianFromMeanAndVarianceOp.LogAverageFactor(sample, mean, 1.0 / precision);			
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
		public static double LogAverageFactor([SkipIfUniform] Gaussian sample, double mean, double precision)
		{
			return LogAverageFactor(sample, Gaussian.PointMass(mean), precision);
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
		public static double LogAverageFactor(double sample, [SkipIfUniform] Gaussian mean, double precision)
		{
			//if(mean.IsPointMass) return LogAverageFactor(sample,mean.Point,precision);
			return LogAverageFactor(Gaussian.PointMass(sample), mean, precision);
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
		public static double LogAverageFactor(double sample, double mean, [SkipIfUniform] Gamma precision)
		{
			if (precision.IsPointMass) return LogAverageFactor(sample, mean, precision.Point);
			if (precision.IsUniform()) return Double.PositiveInfinity;
			return TPdfLn(sample - mean, 2 * precision.Rate, 2 * precision.Shape + 1);
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
		public static double LogAverageFactor(double sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			if (mean.IsPointMass) return LogAverageFactor(sample, mean.Point, precision);
			if (precision.IsPointMass) return LogAverageFactor(sample, mean, precision.Point);
			if (precision.IsUniform()) return Double.PositiveInfinity;
			return LogAverageFactor(Gaussian.PointMass(sample), mean, precision);
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
		public static double LogAverageFactor([SkipIfUniform] Gaussian sample, double mean, [SkipIfUniform] Gamma precision)
		{
			if (precision.IsPointMass) return LogAverageFactor(sample, mean, precision.Point);
			if (precision.IsUniform()) return Double.PositiveInfinity;
			if (sample.IsPointMass) return LogAverageFactor(sample.Point, mean, precision);
			return LogAverageFactor(sample, Gaussian.PointMass(mean), precision);
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
		public static double LogAverageFactor([SkipIfUniform] Gaussian sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			if (precision.IsPointMass) return LogAverageFactor(sample, mean, precision.Point);
			if (precision.IsUniform()) return Double.PositiveInfinity;
			if (sample.IsPointMass && mean.IsPointMass) return LogAverageFactor(sample.Point, mean.Point, precision);
			if (sample.IsUniform() || mean.IsUniform()) return 0.0;
			// this code works even if sample and mean are point masses, but not if any variable is uniform.
			double xm, xv, mm, mv;
			sample.GetMeanAndVariance(out xm, out xv);
			mean.GetMeanAndVariance(out mm, out mv);
			// use quadrature to integrate over the precision
			double[] nodes = new double[QuadratureNodeCount];
			double[] weights = new double[nodes.Length];
			QuadratureNodesAndWeights(precision, nodes, weights);
			for (int i = 0; i < nodes.Length; i++) {
				weights[i] = Math.Log(weights[i]) + Gaussian.GetLogProb(xm, mm, xv + mv + 1.0 / nodes[i]);
			}
			return MMath.LogSumExp(weights);
		}

		/// <summary>
		/// Logarithm of Student T density.
		/// </summary>
		/// <param name="x">sample</param>
		/// <param name="v">variance parameter</param>
		/// <param name="n">degrees of freedom plus 1</param>
		/// <returns></returns>
		public static double TPdfLn(double x, double v, double n)
		{
			return MMath.GammaLn(n * 0.5) - MMath.GammaLn((n - 1) * 0.5) - 0.5 * Math.Log(v * Math.PI) - 0.5 * n * Math.Log(1 + x * x / v);
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
		public static double LogEvidenceRatio(double sample, double mean, double precision)
		{
			return LogAverageFactor(sample, mean, precision);
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
		public static double LogEvidenceRatio(Gaussian sample, Gaussian mean, double precision)
		{
			return 0.0;
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
		public static double LogEvidenceRatio(Gaussian sample, double mean, double precision)
		{
			return 0.0;
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
		public static double LogEvidenceRatio(double sample, [SkipIfUniform] Gaussian mean, double precision)
		{
			return LogAverageFactor(sample, mean, precision);
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
		public static double LogEvidenceRatio(double sample, double mean, [SkipIfUniform] Gamma precision)
		{
			return LogAverageFactor(sample, mean, precision);
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
		public static double LogEvidenceRatio(double sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision)
		{
			return LogAverageFactor(sample, mean, precision);
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
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian sample, [SkipIfUniform] Gaussian mean, [SkipIfUniform] Gamma precision, [Fresh] Gaussian to_sample)
		{
			if (precision.IsPointMass) return LogEvidenceRatio(sample, mean, precision.Point);
			//Gaussian to_Sample = SampleAverageConditional(sample, mean, precision);
			return LogAverageFactor(sample, mean, precision)
				- sample.GetLogAverageOf(to_sample);
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
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian sample, double mean, [SkipIfUniform] Gamma precision, [Fresh] Gaussian to_sample)
		{
			if (precision.IsPointMass) return LogEvidenceRatio(sample, mean, precision.Point);
			//Gaussian to_Sample = SampleAverageConditional(sample, mean, precision);
			return LogAverageFactor(sample, mean, precision)
				- sample.GetLogAverageOf(to_sample);
		}

		//-- VMP ------------------------------------------------------------------------------------------------

		/// <summary>
		/// VMP message to 'sample'
		/// </summary>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision">Constant value for 'precision'.</param>
		/// <returns>The outgoing VMP message to the 'sample' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'sample' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian SampleAverageLogarithm(double mean, double precision)
		{
			return SampleConditional(mean, precision);
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
		public static Gaussian SampleAverageLogarithm([Proper] Gaussian mean, double precision)
		{
			if (precision < 0.0) throw new ArgumentException("precision < 0 (" + precision + ")");
			Gaussian result = new Gaussian();
			result.Precision = precision;
			result.MeanTimesPrecision = precision * mean.GetMean();
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
		public static Gaussian MeanAverageLogarithm(double sample, double precision)
		{
			return MeanConditional(sample, precision);
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
		public static Gaussian MeanAverageLogarithm([Proper] Gaussian sample, double precision)
		{
			return SampleAverageLogarithm(sample, precision);
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
		public static Gaussian SampleAverageLogarithm([Proper] Gaussian mean, [Proper] Gamma precision)
		{
			// The formula is exp(int_prec int_mean p(mean) p(prec) log N(x;mean,1/prec)) =
			// exp(-0.5 E[prec*(x-mean)^2] + const.) =
			// exp(-0.5 E[prec] (x^2 - 2 x E[mean]) + const.) =
			// N(x; E[mean], 1/E[prec])
			Gaussian result = new Gaussian();
			result.Precision = precision.GetMean();
			result.MeanTimesPrecision = result.Precision * mean.GetMean();
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
		public static Gaussian MeanAverageLogarithm([Proper] Gaussian sample, [Proper]Gamma precision)
		{
			return SampleAverageLogarithm(sample, precision);
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
		public static Gaussian SampleAverageLogarithm(double mean, [Proper] Gamma precision)
		{
			Gaussian result = new Gaussian();
			result.Precision = precision.GetMean();
			result.MeanTimesPrecision = result.Precision * mean;
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
		public static Gaussian MeanAverageLogarithm(double sample, [Proper]Gamma precision)
		{
			return SampleAverageLogarithm(sample, precision);
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
		public static Gamma PrecisionAverageLogarithm(double sample, double mean)
		{
			return PrecisionConditional(sample, mean);
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
		public static Gamma PrecisionAverageLogarithm([Proper]Gaussian sample, [Proper]Gaussian mean)
		{
			if (sample.IsUniform()) throw new ImproperMessageException(sample);
			if (mean.IsUniform()) throw new ImproperMessageException(mean);
			// The formula is exp(int_x int_mean p(x) p(mean) log N(x;mean,1/prec)) =
			// exp(-0.5 prec E[(x-mean)^2] + 0.5 log(prec)) =
			// Gamma(prec; 0.5, 0.5*E[(x-mean)^2])
			// E[(x-mean)^2] = E[x^2] - 2 E[x] E[mean] + E[mean^2] = var(x) + (E[x]-E[mean])^2 + var(mean)
			Gamma result = new Gamma();
			result.Shape = 1.5;
			double mx, vx, mm, vm;
			sample.GetMeanAndVariance(out mx, out vx);
			mean.GetMeanAndVariance(out mm, out vm);
			double diff = mx - mm;
			result.Rate = 0.5 * (vx + diff * diff + vm);
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
		public static Gamma PrecisionAverageLogarithm([Proper]Gaussian sample, double mean)
		{
			if (sample.IsUniform()) throw new ImproperMessageException(sample);
			Gamma result = new Gamma();
			result.Shape = 1.5;
			double mx, vx;
			sample.GetMeanAndVariance(out mx, out vx);
			double diff = mx - mean;
			result.Rate = 0.5 * (vx + diff * diff);
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
		public static Gamma PrecisionAverageLogarithm(double sample, [Proper]Gaussian mean)
		{
			return PrecisionAverageLogarithm(mean, sample);
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
		public static double AverageLogFactor([Proper] Gaussian sample, [Proper] Gaussian mean, [Proper] Gamma precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision);
			if (mean.IsPointMass)
				return AverageLogFactor(sample, mean.Point, precision);
			if (precision.IsPointMass)
				return AverageLogFactor(sample, mean, precision.Point);

			return ComputeAverageLogFactor(sample, mean, precision.GetMeanLog(), precision.GetMean());
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
		public static double AverageLogFactor(double sample, double mean, [Proper] Gamma precision)
		{
			if (precision.IsPointMass)
				return AverageLogFactor(sample, mean, precision.Point);
			else {
				double diff = sample - mean;
				return -MMath.LnSqrt2PI + 0.5 * (precision.GetMeanLog() - precision.GetMean() * diff * diff);
			}
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
		public static double AverageLogFactor(double sample, double mean, double precision)
		{
			double diff = sample - mean;
			if (double.IsPositiveInfinity(precision)) return (diff == 0.0) ? 0.0 : double.NegativeInfinity;
			if (precision == 0.0) return 0.0;
			return -MMath.LnSqrt2PI + 0.5 * (Math.Log(precision) - precision * diff * diff);
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
		public static double AverageLogFactor([Proper] Gaussian sample, double mean, double precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision);
			else if (double.IsPositiveInfinity(precision))
				return sample.GetLogProb(mean);
			else if (precision == 0.0)
				return 0.0;
			else
				return ComputeAverageLogFactor(sample, mean, Math.Log(precision), precision);
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
		public static double AverageLogFactor(double sample, [Proper] Gaussian mean, double precision)
		{
			return AverageLogFactor(mean, sample, precision);
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
		public static double AverageLogFactor(double sample, [Proper] Gaussian mean, [Proper] Gamma precision)
		{
			return AverageLogFactor(mean, sample, precision);
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
		public static double AverageLogFactor([Proper] Gaussian sample, double mean, [Proper] Gamma precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision);
			if (precision.IsPointMass)
				return AverageLogFactor(sample, mean, precision.Point);

			return ComputeAverageLogFactor(sample, mean, precision.GetMeanLog(), precision.GetMean());
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
		public static double AverageLogFactor([Proper] Gaussian sample, [Proper] Gaussian mean, double precision)
		{
			if (sample.IsPointMass)
				return AverageLogFactor(sample.Point, mean, precision);
			else if (mean.IsPointMass)
				return AverageLogFactor(sample, mean.Point, precision);
			else if (double.IsPositiveInfinity(precision))
				return sample.GetLogAverageOf(mean);
			else if (precision == 0.0) 
				return 0.0;
			else
				return ComputeAverageLogFactor(sample, mean, Math.Log(precision), precision);
		}

		/// <summary>
		/// Helper method for computing average log factor
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Incoming message from 'mean'.</param>
		/// <param name="precision_Elogx">Expected log value of the incoming message from 'precision'</param>
		/// <param name="precision_Ex">Expected value of incoming message from 'precision'</param>
		/// <returns>Computed average log factor</returns>
		private static double ComputeAverageLogFactor(Gaussian sample, Gaussian mean, double precision_Elogx, double precision_Ex)
		{
			if (precision_Ex == 0.0) throw new ArgumentException("precision == 0");
			if (double.IsPositiveInfinity(precision_Ex)) throw new ArgumentException("precision is infinite");
			double sampleMean, sampleVariance = 0;
			double meanMean, meanVariance = 0;
			sample.GetMeanAndVariance(out sampleMean, out sampleVariance);
			mean.GetMeanAndVariance(out meanMean, out meanVariance);
			double diff = sampleMean - meanMean;
			return -MMath.LnSqrt2PI + 0.5 * (precision_Elogx
											- precision_Ex * (diff * diff + sampleVariance + meanVariance));
		}

		/// <summary>
		/// Helper method for computing average log factor
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="mean">Constant value for 'mean'.</param>
		/// <param name="precision_Elogx">Expected log value of the incoming message from 'precision'</param>
		/// <param name="precision_Ex">Expected value of incoming message from 'precision'</param>
		/// <returns>Computed average log factor</returns>
		private static double ComputeAverageLogFactor(Gaussian sample, double mean, double precision_Elogx, double precision_Ex)
		{
			if (double.IsPositiveInfinity(precision_Ex)) throw new ArgumentException("precision is infinite");
			double sampleMean, sampleVariance;
			sample.GetMeanAndVariance(out sampleMean, out sampleVariance);
			double diff = sampleMean - mean;
			return -MMath.LnSqrt2PI + 0.5 * (precision_Elogx
											- precision_Ex * (diff * diff + sampleVariance));
		}
	}
}
