// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Max"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Max", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class MaxGaussianOp
	{
		/// <summary>
		/// Static flag to force a proper distribution
		/// </summary>
		public static bool ForceProper;

		public static double LogAverageFactor(double max, double a, double b)
		{
			return (max == Factor.Max(a,b)) ? 0.0 : Double.NegativeInfinity;
		}
		public static double LogEvidenceRatio(double max, double a, double b) { return LogAverageFactor(max, a, b); }
		public static double AverageLogFactor(double max, double a, double b) { return LogAverageFactor(max, a, b); }

		// logw1 = N(mx;m1,vx+v1) phi((mx1 - m2)/sqrt(vx1+v2))
		// a1 = N(mx1;m2,vx1+v2)/phi
		internal static void ComputeStats(Gaussian max, Gaussian a, Gaussian b, out double logz,
			out double logw1, out double a1, out double vx1, out double mx1,
			out double logw2, out double a2, out double vx2, out double mx2)
		{
			double m1, v1, m2, v2;
			a.GetMeanAndVariance(out m1, out v1);
			b.GetMeanAndVariance(out m2, out v2);
			if (max.IsPointMass) {
				vx1 = 0.0;
				mx1 = max.Point;
				vx2 = 0.0;
				mx2 = max.Point;
				if (b.IsPointMass) {
					if (b.Point > max.Point) throw new AllZeroException();
					else if (b.Point == max.Point) {
						// the factor reduces to the constraint (max.Point > a)
						logw1 = Double.NegativeInfinity;
						logw2 = MMath.NormalCdfLn((max.Point - m1)/Math.Sqrt(v1));
						logz = logw2;
						a1 = 0;
						a2 = Math.Exp(Gaussian.GetLogProb(max.Point, m1, v1) - logw2);
						return;
					} else {
						// b.Point < max.Point
						// the factor reduces to the constraint (a == max.Point)
						throw new NotImplementedException();
					}
				} else if (a.IsPointMass) throw new NotImplementedException();
			} else {
				if (a.IsPointMass) {
					vx1 = 0.0;
					mx1 = a.Point;
				} else {
					vx1 = 1.0 / (max.Precision + a.Precision);
					mx1 = vx1 * (max.MeanTimesPrecision + a.MeanTimesPrecision);
				}
				if (b.IsPointMass) {
					vx2 = 0.0;
					mx2 = b.Point;
				} else {
					vx2 = 1.0 / (max.Precision + b.Precision);
					mx2 = vx2 * (max.MeanTimesPrecision + b.MeanTimesPrecision);
				}
			}
			logw1 = max.GetLogAverageOf(a);
			double logPhi1 = MMath.NormalCdfLn((mx1 - m2) / Math.Sqrt(vx1 + v2));
			logw1 += logPhi1;

			logw2 = max.GetLogAverageOf(b);
			double logPhi2 = MMath.NormalCdfLn((mx2 - m1) / Math.Sqrt(vx2 + v1));
			logw2 += logPhi2;

			logz = MMath.LogSumExp(logw1, logw2);

			double logN1 = Gaussian.GetLogProb(mx1, m2, vx1 + v2);
			a1 = Math.Exp(logN1 - logPhi1);
			double logN2 = Gaussian.GetLogProb(mx2, m1, vx2 + v1);
			a2 = Math.Exp(logN2 - logPhi2);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,a,b) p(max,a,b) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Gaussian max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			double logw1, a1, vx1, mx1;
			double logw2, a2, vx2, mx2;
			double logz;
			ComputeStats(max, a, b, out logz, out logw1, out a1, out vx1, out mx1,
				out logw2, out a2, out vx2, out mx2);
			return logz;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,b) p(max,b) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Gaussian max, double a, [Proper] Gaussian b)
		{
			return LogAverageFactor(max, Gaussian.PointMass(a), b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,a) p(max,a) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogAverageFactor([SkipIfUniform] Gaussian max, [Proper] Gaussian a, double b)
		{
			return LogAverageFactor(max, a, Gaussian.PointMass(b));
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			return LogAverageFactor(Gaussian.PointMass(max), a, b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double max, double a, [Proper] Gaussian b)
		{
			return LogAverageFactor(Gaussian.PointMass(max), Gaussian.PointMass(a), b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(max,a,b))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogAverageFactor(double max, [Proper] Gaussian a, double b)
		{
			return LogAverageFactor(Gaussian.PointMass(max), a, Gaussian.PointMass(b));
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,a,b) p(max,a,b) factor(max,a,b) / sum_max p(max) messageTo(max))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			Gaussian to_max = MaxAverageConditional(max, a, b);
			return LogAverageFactor(max, a, b) 
				- to_max.GetLogAverageOf(max);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,b) p(max,b) factor(max,a,b) / sum_max p(max) messageTo(max))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian max, double a, [Proper] Gaussian b)
		{
			Gaussian to_max = MaxAverageConditional(max, a, b);
			return LogAverageFactor(max, a, b) 
				- to_max.GetLogAverageOf(max);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(max,a) p(max,a) factor(max,a,b) / sum_max p(max) messageTo(max))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian max, [Proper] Gaussian a, double b)
		{
			Gaussian to_max = MaxAverageConditional(max, a, b);
			return LogAverageFactor(max, a, b) 
				- to_max.GetLogAverageOf(max);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(max,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(double max, [Proper] Gaussian a, [Proper] Gaussian b) { return LogAverageFactor(max, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(max,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(double max, double a, [Proper] Gaussian b) { return LogAverageFactor(max, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(max,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static double LogEvidenceRatio(double max, [Proper] Gaussian a, double b) { return LogAverageFactor(max, a, b); }

		[Skip]
		public static Gaussian MaxAverageConditionalInit()
		{
			return Gaussian.Uniform();
		}

		/// <summary>
		/// EP message to 'max'
		/// </summary>
		/// <param name="max">Incoming message from 'max'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'max' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'max' as the random arguments are varied.
		/// The formula is <c>proj[p(max) sum_(b) p(b) factor(max,a,b)]/p(max)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian MaxAverageConditional(Gaussian max, double a, [Proper] Gaussian b)
		{
			return MaxAverageConditional(max, Gaussian.PointMass(a), b);
		}
		/// <summary>
		/// EP message to 'max'
		/// </summary>
		/// <param name="max">Incoming message from 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'max' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'max' as the random arguments are varied.
		/// The formula is <c>proj[p(max) sum_(a) p(a) factor(max,a,b)]/p(max)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian MaxAverageConditional(Gaussian max, [Proper] Gaussian a, double b)
		{
			return MaxAverageConditional(max, a, Gaussian.PointMass(b));
		}

		/// <summary>
		/// EP message to 'max'
		/// </summary>
		/// <param name="max">Incoming message from 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'max' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'max' as the random arguments are varied.
		/// The formula is <c>proj[p(max) sum_(a,b) p(a,b) factor(max,a,b)]/p(max)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian MaxAverageConditional(Gaussian max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			if (max.IsPointMass) return Gaussian.Uniform();
			// the following code works correctly even if max is uniform or improper.
			if (!a.IsProper()) throw new ImproperMessageException(a);
			if (!b.IsProper()) throw new ImproperMessageException(b);
			double m1, v1, m2, v2;
			a.GetMeanAndVariance(out m1, out v1);
			b.GetMeanAndVariance(out m2, out v2);

			double logw1, a1, vx1, mx1;
			double logw2, a2, vx2, mx2;
			double logz;
			ComputeStats(max, a, b, out logz, out logw1, out a1, out vx1, out mx1,
				out logw2, out a2, out vx2, out mx2);
			double w1 = Math.Exp(logw1 - logz);
			double w2 = Math.Exp(logw2 - logz);
			// the posterior is a mixture model with weights exp(logw1-logz), exp(logw2-logz) and distributions
			// N(x; mx1, vx1) phi((x - m2)/sqrt(v2)) / phi((mx1 - m2)/sqrt(vx1 + v2))
			// the moments of the posterior are computed via the moments of these two components.
			double mc1 = mx1 + a1 * vx1;
			double mc2 = mx2 + a2 * vx2;
			double m = w1 * mc1 + w2 * mc2;
			double r1 = (mx1 - m2) / (vx1 + v2);
			double r2 = (mx2 - m1) / (vx2 + v1);
			double beta1 = a1 * (a1 + r1);
			double beta2 = a2 * (a2 + r2);
			double vc1 = vx1 * (1 - vx1 * beta1);
			double vc2 = vx2 * (1 - vx2 * beta2);
			double diff = mc1 - mc2;
			double v = w1 * vc1 + w2 * vc2 + w1 * w2 * diff * diff;
			Gaussian result = new Gaussian(m, v);
			if (ForceProper) result.SetToRatioProper(result, max);
			else result.SetToRatio(result, max);
			if (Double.IsNaN(result.Precision) || Double.IsNaN(result.MeanTimesPrecision)) throw new Exception("result is NaN");
			return result;
		}

#if false
		public static void ComputeStats(Gaussian max, Gaussian a, Gaussian b, out double logz,
			out double logw1, out double logPhi1, out double logu1, out double vx1, out double mx1,
			out double logw2, out double logPhi2, out double logu2, out double vx2, out double mx2)
		{
			double mx, vx, ma, va, mb, vb;
			max.GetMeanAndVariance(out mx, out vx);
			a.GetMeanAndVariance(out ma, out va);
			b.GetMeanAndVariance(out mb, out vb);
			if (false) {
				vx1 = 1.0 / (1.0 / vx + 1.0 / va);
				mx1 = vx1 * (mx / vx + ma / va);
				vx2 = 1.0 / (1.0 / vx + 1.0 / vb);
				mx2 = vx2 * (mx / vx + mb / vb);
			} else {
				if(max.IsPointMass || a.IsPointMass || b.IsPointMass) throw new NotImplementedException();
				vx1 = 1.0/(max.Precision + a.Precision);
				mx1 = vx1*(max.MeanTimesPrecision + a.MeanTimesPrecision);
				vx2 = 1.0/(max.Precision + b.Precision);
				mx2 = vx2*(max.MeanTimesPrecision + b.MeanTimesPrecision);
			}
			logw1 = max.GetLogAverageOf(a);
			logPhi1 = MMath.NormalCdfLn((mx1 - mb) / Math.Sqrt(vx1 + vb));
			logu1 = Gaussian.GetLogProb(mx1, mb, vx1 + vb);

			logw2 = max.GetLogAverageOf(b);
			logPhi2 = MMath.NormalCdfLn((mx2 - ma) / Math.Sqrt(vx2 + va));
			logu2 = Gaussian.GetLogProb(mx2, ma, vx2 + va);

			logz = MMath.LogSumExp(logw1 + logPhi1, logw2 + logPhi2);
		}
#endif

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(max,b) p(max,b) factor(max,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			if (a.IsPointMass) return Gaussian.Uniform();
			if (max.IsUniform()) return Gaussian.Uniform();
			if (!a.IsProper()) throw new ImproperMessageException(a);
			if (!b.IsProper()) throw new ImproperMessageException(b);
			double m1, v1, m2, v2;
			a.GetMeanAndVariance(out m1, out v1);
			b.GetMeanAndVariance(out m2, out v2);

#if true
			double logw1, a1, vx1, mx1;
			double logw2, a2, vx2, mx2;
			double logz;
			ComputeStats(max, a, b, out logz, out logw1, out a1, out vx1, out mx1,
				out logw2, out a2, out vx2, out mx2);
			double w1 = Math.Exp(logw1 - logz);
			double w2 = Math.Exp(logw2 - logz);
			// the posterior is a mixture model with weights exp(logw1-logz), exp(logw2-logz) and distributions
			// N(a; mx1, vx1) phi((a - m2)/sqrt(v2)) / phi((mx1 - m2)/sqrt(vx1 + v2))
			// N(a; m1, v1) phi((mx2 - a)/sqrt(vx2)) / phi((mx2 - m1)/sqrt(vx2 + v1))
			// the moments of the posterior are computed via the moments of these two components.
			double mc1 = mx1 + a1 * vx1;
			a2 = -a2;
			double mc2 = m1 + a2 * v1;
			double m = w1 * mc1 + w2 * mc2;
			double r1 = (mx1 - m2) / (vx1 + v2);
			double r2 = (mx2 - m1) / (vx2 + v1);
			double beta1 = a1 * (a1 + r1);
			double beta2 = a2 * (a2 - r2);
			double vc1 = vx1 * (1 - vx1 * beta1);
			double vc2 = v1 * (1 - v1 * beta2);
			double diff = mc1 - mc2;
			double v = w1 * vc1 + w2 * vc2 + w1 * w2 * diff * diff;
			Gaussian result = new Gaussian(m, v);
			if (ForceProper) result.SetToRatioProper(result, a);
			else result.SetToRatio(result, a);
#else
			double logw1, logPhi1, logu1, vx1, mx1;
			double logw2, logPhi2, logu2, vx2, mx2;
			double logz;
			ComputeStats(max, a, b, out logz, out logw1, out logPhi1, out logu1, out vx1, out mx1,
				out logw2, out logPhi2, out logu2, out vx2, out mx2);
			if (a.IsPointMass || max.IsPointMass || b.IsPointMass) throw new NotImplementedException();
			//double ra = (mx-ma)/(vx+va);
			double ra = (max.MeanTimesPrecision * a.Precision - a.MeanTimesPrecision * max.Precision) / (a.Precision + max.Precision);
			//double rb = (mx-mb)/(vx+vb);
			double rb = (max.MeanTimesPrecision * b.Precision - b.MeanTimesPrecision * max.Precision) / (b.Precision + max.Precision);
			double vxinva = a.Precision / (a.Precision + max.Precision);
			double vxinvb = b.Precision / (b.Precision + max.Precision);
			double mx, vx;
			max.GetMeanAndVariance(out mx, out vx);
			double invxa = 1.0 / (vx + v1);
			double invxb = 1.0 / (vx + v2);
			double alpha = Math.Exp(logw1 + logPhi1 - logz) * ra + Math.Exp(logw1 + logu1 - logz) * vxinva - Math.Exp(logw2 + logu2 - logz);
			double dvx = Math.Exp(logw1 + logPhi1 - logz) * (ra * ra - invxa) + Math.Exp(logw1 + logu1 - logz) * (
				2 * ra * vxinva - (mx1 - m2) / (vx1 + v2) * vxinva * vxinva) - Math.Exp(logw2 + logu2 - logz) * (mx2 - m1) / (vx2 + v1);
			double beta = alpha * alpha - dvx;

			double tau = a.MeanTimesPrecision;
			double prec = a.Precision;
			double weight = beta / (prec - beta);
			if (ForceProper && weight < 0) weight = 0;
			Gaussian result = new Gaussian();
			// eq (31) in EP quickref; same as inv(inv(beta)-inv(prec))
			result.Precision = prec * weight;
			// eq (30) in EP quickref times above and simplified
			result.MeanTimesPrecision = weight * (tau + alpha) + alpha;
#endif
			if (Double.IsNaN(result.Precision) || Double.IsNaN(result.MeanTimesPrecision)) throw new Exception("result is NaN");
			return result;
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(b) p(b) factor(max,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional(double max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			return AAverageConditional(Gaussian.PointMass(max), a, b);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(max) p(max) factor(max,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian max, [Proper] Gaussian a, double b)
		{
			return AAverageConditional(max, a, Gaussian.PointMass(b));
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional(double max, [Proper] Gaussian a, double b)
		{
			return AAverageConditional(Gaussian.PointMass(max), a, Gaussian.PointMass(b));
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(max,a) p(max,a) factor(max,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			return AAverageConditional(max, b, a);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(a) p(a) factor(max,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="a"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional(double max, [Proper] Gaussian a, [Proper] Gaussian b)
		{
			return AAverageConditional(max, b, a);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="max">Incoming message from 'max'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(max) p(max) factor(max,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="max"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian max, double a, [Proper] Gaussian b)
		{
			return AAverageConditional(max, b, a);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="max">Constant value for 'max'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="b"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional(double max, double a, [Proper] Gaussian b)
		{
			return AAverageConditional(max, b, a);
		}
	}
}
