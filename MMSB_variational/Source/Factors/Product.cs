// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing WrappedGaussian messages for <see cref="Factor.Product(double,double)"/> and <see cref="Factor.Ratio"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[FactorMethod(new string[] { "A", "Product", "B" }, typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Experimental)]
	public static class WrappedGaussianProductOp
	{
		public static WrappedGaussian AAverageConditional([SkipIfUniform] WrappedGaussian Product, double B, WrappedGaussian result)
		{
			result.Period = Product.Period/B;
			result.Gaussian = GaussianProductOp.AAverageConditional(Product.Gaussian, B);
			result.Normalize();
			return result;
		}
		public static WrappedGaussian BAverageConditional([SkipIfUniform] WrappedGaussian Product, double A, WrappedGaussian result)
		{
			return AAverageConditional(Product, A, result);
		}
		public static WrappedGaussian AAverageConditional(double Product, double B, WrappedGaussian result)
		{
			if (B == 0) {
				if (Product != 0) throw new AllZeroException();
				result.SetToUniform();
			} else result.Point = Product / B;
			result.Normalize();
			return result;
		}

		// ----------------------------------------------------------------------------------------------------------------------
		// VMP
		// ----------------------------------------------------------------------------------------------------------------------

		[Skip]
		public static double AverageLogFactor(WrappedGaussian product) { return 0.0; }
		public static WrappedGaussian ProductAverageLogarithm([SkipIfUniform] WrappedGaussian A, double B, WrappedGaussian result)
		{
			double m, v;
			A.Gaussian.GetMeanAndVariance(out m, out v);
			result.Gaussian.SetMeanAndVariance(B*m, B*B*v);
			double period = B*A.Period;
			if (period != result.Period) {
				double ratio = period / result.Period;
				double intRatio = Math.Round(ratio);
				if (Math.Abs(ratio - intRatio) > result.Period*1e-4) throw new ArgumentException("B*A.Period ("+period+") is not a multiple of result.Period ("+result.Period+")");
				// if period is a multiple of result.Period, then wrapping to result.Period is equivalent to first wrapping to period, then to result.Period.
			}
			result.Normalize();
			return result;
		}
		public static WrappedGaussian ProductAverageLogarithm(double A, [SkipIfUniform] WrappedGaussian B, WrappedGaussian result)
		{
			return ProductAverageLogarithm(B, A, result);
		}
		public static WrappedGaussian AAverageLogarithm([SkipIfUniform] WrappedGaussian Product, double B, WrappedGaussian result)
		{
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B, result);
			return AAverageConditional(Product, B, result);
		}
		public static WrappedGaussian BAverageLogarithm([SkipIfUniform] WrappedGaussian Product, double A, WrappedGaussian result)
		{
			return AAverageLogarithm(Product, A, result);
		}
		public static WrappedGaussian AAverageLogarithm(double Product, double B, WrappedGaussian result)
		{
			return AAverageConditional(Product, B, result);
		}
	}
	/// <summary>
	/// Provides outgoing Gaussian messages for <see cref="Factor.Product(double,double)"/> and <see cref="Factor.Ratio"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double), Default=true)]
	[FactorMethod(new string[] { "A", "Product", "B" }, typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class GaussianProductOp
	{
		/// <summary>
		/// The number of quadrature nodes used to compute the messages.
		/// Reduce this number to save time in exchange for less accuracy.
		/// Must be an odd number.
		/// </summary>
		public static int QuadratureNodeCount = 1001; // must be odd to avoid A=0

		/// <summary>
		/// Force proper messages
		/// </summary>
		public static bool ForceProper;

		[Skip]
		public static Gaussian ProductAverageConditionalInit(Gaussian A, Gaussian B)
		{
			return Gaussian.Uniform();
		}

		/// <summary>
		/// EP message to 'product'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[p(product) sum_(a,b) p(a,b) factor(product,a,b)]/p(product)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gaussian ProductAverageConditional(Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			if (A.IsPointMass) return ProductAverageConditional(A.Point, B);
			if (B.IsPointMass) return ProductAverageConditional(A, B.Point);
			if (Product.IsPointMass) return Gaussian.Uniform();
			if (Product.IsUniform()) return GaussianProductVmpOp.ProductAverageLogarithm(A, B);
			double mA, vA;
			A.GetMeanAndVariance(out mA, out vA);
			double mB, vB;
			B.GetMeanAndVariance(out mB, out vB);
			double mProduct, vProduct;
			Product.GetMeanAndVariance(out mProduct, out vProduct);
			// algorithm: quadrature on A from -1 to 1, plus quadrature on 1/A from -1 to 1.
			double z = 0, sumX = 0, sumX2 = 0;
			for (int i = 0; i <= QuadratureNodeCount; i++) {
				double a = (2.0 * i) / QuadratureNodeCount - 1;
				double logfA = Gaussian.GetLogProb(mProduct, a * mB, vProduct + a * a * vB) + Gaussian.GetLogProb(a, mA, vA);
				double fA = Math.Exp(logfA);

				z += fA;
				double b = (mB * vProduct + a * mProduct * vB) / (vProduct + a * a * vB);
				double b2 = b * b + (vProduct * vB) / (vProduct + a * a * vB);
				double x = a * b;
				double x2 = a * a * b2;
				sumX += x * fA;
				sumX2 += x2 * fA;

				double invA = a;
				a = 1.0 / invA;
				double logfInvA = Gaussian.GetLogProb(mProduct * invA, mB, vProduct * invA * invA + vB) + Gaussian.GetLogProb(a, mA, vA) - Math.Log(Math.Abs(invA + Double.Epsilon));
				double fInvA = Math.Exp(logfInvA);
				z += fInvA;
				b = (mB * vProduct + a * mProduct * vB) / (vProduct + a * a * vB);
				b2 = b * b + (vProduct * vB) / (vProduct + a * a * vB);
				x = a * b;
				x2 = a * a * b2;
				sumX += x * fInvA;
				sumX2 += x2 * fInvA;
			}
			double mean = sumX / z;
			double var = sumX2 / z - mean * mean;
			Gaussian result = Gaussian.FromMeanAndVariance(mean, var);
			if (ForceProper) result.SetToRatioProper(result, Product);
			else result.SetToRatio(result, Product);
			return result;
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(product,b) p(product,b) factor(product,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, Gaussian A, [SkipIfUniform] Gaussian B)
		{
			if (B.IsPointMass) return AAverageConditional(Product, B.Point);
			if (A.IsPointMass || Product.IsUniform()) return Gaussian.Uniform();
			Gaussian result = new Gaussian();
			// algorithm: quadrature on A from -1 to 1, plus quadrature on 1/A from -1 to 1.
			double mProduct, vProduct;
			Product.GetMeanAndVariance(out mProduct, out vProduct);
			double mA, vA;
			A.GetMeanAndVariance(out mA, out vA);
			double mB, vB;
			B.GetMeanAndVariance(out mB, out vB);
			double z = 0, sumA = 0, sumA2 = 0;
			for (int i = 0; i <= QuadratureNodeCount; i++) {
				double a = (2.0 * i) / QuadratureNodeCount - 1;
				double logfA = Gaussian.GetLogProb(mProduct, a * mB, vProduct + a * a * vB) + Gaussian.GetLogProb(a, mA, vA);
				double fA = Math.Exp(logfA);
				z += fA;
				sumA += a * fA;
				sumA2 += a * a * fA;

				double invA = a;
				a = 1.0 / invA;
				double logfInvA = Gaussian.GetLogProb(mProduct * invA, mB, vProduct * invA * invA + vB) + Gaussian.GetLogProb(a, mA, vA) - Math.Log(Math.Abs(invA + Double.Epsilon));
				double fInvA = Math.Exp(logfInvA);
				z += fInvA;
				sumA += a * fInvA;
				sumA2 += a * a * fInvA;
			}
			double mean = sumA / z;
			double var = sumA2 / z - mean * mean;
			result.SetMeanAndVariance(mean, var);
			if (ForceProper) result.SetToRatioProper(result, A);
			else result.SetToRatio(result, A);
			return result;
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(b) p(b) factor(product,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gaussian AAverageConditional(double Product, Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return AAverageConditional(Gaussian.PointMass(Product), A, B);
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(product,a) p(product,a) factor(product,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, Gaussian B)
		{
			return AAverageConditional(Product, B, A);
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(a) p(a) factor(product,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[Quality(QualityBand.Experimental)]
		public static Gaussian BAverageConditional(double Product, [SkipIfUniform] Gaussian A, Gaussian B)
		{
			return AAverageConditional(Product, B, A);
		}

		/// <summary>
		/// EP message to 'product'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[p(product) sum_(b) p(b) factor(product,a,b)]/p(product)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gaussian ProductAverageConditional(double A, [SkipIfUniform] Gaussian B)
		{
			return GaussianProductVmpOp.ProductAverageLogarithm(A, B);
		}
		/// <summary>
		/// EP message to 'product'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[p(product) sum_(a) p(a) factor(product,a,b)]/p(product)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, double B)
		{
			return ProductAverageConditional(B, A);
		}

		/// <summary>
		/// EP message to 'product'
		/// </summary>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'product' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian ProductAverageConditional(double a, double b)
		{
			return Gaussian.PointMass(a * b);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(product) p(product) factor(product,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, double B)
		{
			Gaussian result = new Gaussian();
			if (Product.IsPointMass) return AAverageConditional(Product.Point, B);
			// (m - ab)^2/v = (a^2 b^2 - 2abm + m^2)/v
			// This code works correctly even if B=0 or Product is uniform. 
			result.Precision = B * B * Product.Precision;
			result.MeanTimesPrecision = B * Product.MeanTimesPrecision;
			return result;
		}
		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageConditional(double Product, double B)
		{
			if (B == 0) {
				if (Product != 0) throw new AllZeroException();
				return Gaussian.Uniform();
			} else return Gaussian.PointMass(Product / B);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(product) p(product) factor(product,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, double A)
		{
			return AAverageConditional(Product, A);
		}
		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing EP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian BAverageConditional(double Product, double A)
		{
			return AAverageConditional(Product, A);
		}

#if false
    public static double AverageValueLn([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			double pm, pv;
			Product.GetMeanAndVariance(out pm, out pv);
			// find the variable with most precision.
			Gaussian Wide, Thin;
			if (A.Precision > B.Precision) {
				Wide = B;
				Thin = A;
			} else {
				Wide = A;
				Thin = B;
			}
			if (Thin.IsPointMass) {
				double am, av;
				Wide.GetMeanAndVariance(out am, out av);
				double b = Thin.Point;
				return Gaussian.EvaluateLn(pm, b * am, pv + b * b * av);
			} else {
				// use quadrature to integrate over B
				throw new NotImplementedException();
			}
		}
#endif
	}

	/// <summary>
	/// This class allows EP to process the product factor as if running VMP, as required by Stern's algorithm.
	/// </summary>
	/// <remarks>
	/// This algorithm comes from "Matchbox: Large Scale Online Bayesian Recommendations" by David Stern, Ralf Herbrich, and Thore Graepel, WWW 2009.
	/// </remarks>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class GaussianProductOp_SHG09
	{
		// IgnoreDependency is needed to get good schedules
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B, [IgnoreDependency] Gaussian to_A, [IgnoreDependency] Gaussian to_B)
		{
			//return ProductAverageConditional(A.GetMean(), B, to_B);
			return GaussianProductVmpOp.ProductAverageLogarithm(A*to_A, B*to_B);
		}
		public static Gaussian ProductAverageConditional(double A, [SkipIfUniform] Gaussian B, [IgnoreDependency] Gaussian to_B)
		{
			return GaussianProductVmpOp.ProductAverageLogarithm(A, B*to_B);
		}
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, double B, [IgnoreDependency] Gaussian to_A)
		{
			return GaussianProductVmpOp.ProductAverageLogarithm(A*to_A, B);
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, [Proper/*, Fresh*/] Gaussian B, [IgnoreDependency] Gaussian to_B)
		{
			return GaussianProductVmpOp.AAverageLogarithm(Product, B*to_B);
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, double B)
		{
			return GaussianProductVmpOp.AAverageLogarithm(Product, B);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, [Proper/*, Fresh*/] Gaussian A, [IgnoreDependency] Gaussian to_A)
		{
			//return BAverageConditional(Product, A.GetMean());
			return AAverageConditional(Product, A, to_A);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, double A)
		{
			return AAverageConditional(Product, A);
		}
		[Skip]
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return 0.0;
		}
	}

	/// <summary>
	/// This class allows EP to process the product factor as a linear factor.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Experimental)]
	[Buffers("weights")]
	public static class GaussianProductOp3
	{
		public static Vector Weights(double A, Gaussian B)
		{
			Vector weights = Vector.Zero(4);
			weights[1] = A;
			return weights;
		}
		public static Vector Weights(Gaussian A, double B)
		{
			Vector weights = Vector.Zero(4);
			weights[0] = B;
			return weights;
		}
		public static Vector Weights(Gaussian A, Gaussian B, Gaussian to_A, Gaussian to_B)
		{
			if (A.IsPointMass) return Weights(A.Point, B);
			if (B.IsPointMass) return Weights(A, B.Point);
			A *= to_A;
			B *= to_B;
			double ma, va, mb, vb;
			A.GetMeanAndVariance(out ma, out va);
			B.GetMeanAndVariance(out mb, out vb);
			double ma2 = va + ma*ma;
			double mb2 = vb + mb*mb;
			Vector w = Vector.Zero(3);
			w[0] = ma2*mb;
			w[1] = mb2*ma;
			w[2] = ma*mb;
			PositiveDefiniteMatrix M = new PositiveDefiniteMatrix(3, 3);
			M[0, 0] = ma2;
			M[0, 1] = ma*mb;
			M[0, 2] = ma;
			M[1, 0] = ma*mb;
			M[1, 1] = mb2;
			M[1, 2] = mb;
			M[2, 0] = ma;
			M[2, 1] = mb;
			M[2, 2] = 1;
			w = w.PredivideBy(M);
			Vector weights = Vector.Zero(4);
			weights[0] = w[0];
			weights[1] = w[1];
			weights[2] = w[2];
			weights[3] = ma2*mb2 - w[0]*ma2*mb - w[1]*mb2*ma - w[2]*ma*mb;
			if (weights[3] < 0) weights[3] = 0;
			if (false) {
				// debugging
				GaussianEstimator est = new GaussianEstimator();
				for (int i = 0; i < 10000; i++) {
					double sa = A.Sample();
					double sb = B.Sample();
					double f = sa*sb;
					double g = sa*weights[0] + sb*weights[1] + weights[2];
					est.Add(f-g);
				}
				Console.WriteLine(weights);
				Console.WriteLine(est.GetDistribution(new Gaussian()));
			}
			return weights;
		}
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B, [Fresh] Vector weights)
		{
			// factor is product = N(w[0]*a + w[1]*b + w[2], w[3])
			double v = weights[3];
			Gaussian m = DoublePlusOp.SumAverageConditional(GaussianProductOp.ProductAverageConditional(weights[0], A), GaussianProductOp.ProductAverageConditional(weights[1], B));
			m = DoublePlusOp.SumAverageConditional(m, weights[2]);
			return GaussianFromMeanAndVarianceOp.SampleAverageConditional(m, v);
		}
		public static Gaussian ProductAverageConditional(double A, [SkipIfUniform] Gaussian B, [Fresh] Vector weights)
		{
			return ProductAverageConditional(Gaussian.PointMass(A), B, weights);
		}
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, double B, [Fresh] Vector weights)
		{
			return ProductAverageConditional(A, Gaussian.PointMass(B), weights);
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian B, [Fresh] Vector weights)
		{
			// factor is product = N(w[0]*a + w[1]*b + w[2], w[3])
			double v = weights[3];
			Gaussian sum_B = GaussianFromMeanAndVarianceOp.MeanAverageConditional(Product, v);
			sum_B = DoublePlusOp.AAverageConditional(sum_B, weights[2]);
			Gaussian scale_B = DoublePlusOp.AAverageConditional(sum_B, GaussianProductOp.ProductAverageConditional(weights[1], B));
			return GaussianProductOp.AAverageConditional(scale_B, weights[0]);
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, double B, [Fresh] Vector weights)
		{
			return AAverageConditional(Product, Gaussian.PointMass(B), weights);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [Fresh] Vector weights)
		{
			// factor is product = N(w[0]*a + w[1]*b + w[2], w[3])
			double v = weights[3];
			Gaussian sum_B = GaussianFromMeanAndVarianceOp.MeanAverageConditional(Product, v);
			sum_B = DoublePlusOp.AAverageConditional(sum_B, weights[2]);
			Gaussian scale_B = DoublePlusOp.AAverageConditional(sum_B, GaussianProductOp.ProductAverageConditional(weights[0], A));
			return GaussianProductOp.AAverageConditional(scale_B, weights[1]);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, double A, [Fresh] Vector weights)
		{
			return BAverageConditional(Product, Gaussian.PointMass(A), weights);
		}
		[Skip]
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return 0.0;
		}
	}

	/// <summary>
	/// This class allows EP to process the product factor using an approximation to the integral Z.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Experimental)]
	public static class GaussianProductOp4
	{
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return GaussianProductVmpOp.ProductAverageLogarithm(A, B);
		}
		public static Gaussian ProductAverageConditional(double A, [SkipIfUniform] Gaussian B)
		{
			return ProductAverageConditional(Gaussian.PointMass(A), B);
		}
		public static Gaussian ProductAverageConditional([SkipIfUniform] Gaussian A, double B)
		{
			return ProductAverageConditional(A, Gaussian.PointMass(B));
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, Gaussian A, [SkipIfUniform] Gaussian B)
		{
			if (A.IsPointMass) return Gaussian.Uniform();
			double ma, va;
			A.GetMeanAndVariance(out ma, out va);
			double mb, vb;
			B.GetMeanAndVariance(out mb, out vb);
			double mx, vx;
			Product.GetMeanAndVariance(out mx, out vx);
			double diff = mx - ma*mb;
			double prec = 1/(vx + va*vb + va*mb*mb + vb*ma*ma);
			//if (prec < 1e-14) return Gaussian.Uniform();
			double alpha = prec*(-vb*ma + prec*diff*diff*ma*vb + diff*mb);
			double beta = alpha*alpha + prec*(1 - diff*diff*prec)*(vb+mb*mb);
			//if (beta == 0) return Gaussian.Uniform();
			if (double.IsNaN(alpha) || double.IsNaN(beta)) throw new Exception("alpha is nan");
			double r = beta/(A.Precision - beta);
			Gaussian result = new Gaussian();
			result.Precision = r*A.Precision;
			result.MeanTimesPrecision = r*(alpha + A.MeanTimesPrecision) + alpha;
			//Gaussian result = new Gaussian(ma + alpha/beta, 1/beta - va);
			if (double.IsNaN(result.Precision) || double.IsNaN(result.MeanTimesPrecision)) throw new Exception("result is nan");
			return result;
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, Gaussian A, double B)
		{
			return AAverageConditional(Product, A, Gaussian.PointMass(B));
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, Gaussian B)
		{
			return AAverageConditional(Product, B, A);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, double A, Gaussian B)
		{
			return BAverageConditional(Product, B, Gaussian.PointMass(A));
		}
		[Skip]
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return 0.0;
		}
	}

	/// <summary>
	/// This class allows EP to process the product factor using a log-normal approximation to the input distributions
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Experimental)]
	public static class GaussianProductOp5
	{
		public static Gaussian GetExpMoments(Gaussian x)
		{
			double m,v;
			x.GetMeanAndVariance(out m, out v);
			return new Gaussian(Math.Exp(m+v/2), Math.Exp(2*m+v)*(Math.Exp(v)-1));
		}
		public static Gaussian GetLogMoments(Gaussian x)
		{
			double m,v;
			x.GetMeanAndVariance(out m, out v);
			double lv = Math.Log(v/(m*m) + 1);
			double lm = Math.Log(Math.Abs(m)) - lv/2;
			return new Gaussian(lm, lv);
		}
		public static Gaussian ProductAverageConditional(Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			Gaussian logA = GetLogMoments(A);
			Gaussian logB = GetLogMoments(B);
			Gaussian logMsg = DoublePlusOp.SumAverageConditional(logA, logB);
			if (false) {
				Gaussian logProduct = GetLogMoments(Product);
				Gaussian logPost = logProduct * logMsg;
				return GetExpMoments(logPost) / Product;
			} else {
				return GetExpMoments(logMsg);
			}
			//return GaussianProductVmpOp.ProductAverageLogarithm(A, B);
		}
		public static Gaussian ProductAverageConditional(Gaussian Product, double A, [SkipIfUniform] Gaussian B)
		{
			return ProductAverageConditional(Product, Gaussian.PointMass(A), B);
		}
		public static Gaussian ProductAverageConditional(Gaussian Product, [SkipIfUniform] Gaussian A, double B)
		{
			return ProductAverageConditional(Product, A, Gaussian.PointMass(B));
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, Gaussian A, [SkipIfUniform] Gaussian B)
		{
			if (A.IsPointMass) return Gaussian.Uniform();
			Gaussian logA = GetLogMoments(A);
			Gaussian logB = GetLogMoments(B);
			Gaussian logProduct = GetLogMoments(Product);
			Gaussian logMsg = DoublePlusOp.AAverageConditional(logProduct, logB);
			return GetExpMoments(logMsg);
		}
		public static Gaussian AAverageConditional([SkipIfUniform] Gaussian Product, Gaussian A, double B)
		{
			return AAverageConditional(Product, A, Gaussian.PointMass(B));
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, Gaussian B)
		{
			return AAverageConditional(Product, B, A);
		}
		public static Gaussian BAverageConditional([SkipIfUniform] Gaussian Product, double A, Gaussian B)
		{
			return BAverageConditional(Product, B, Gaussian.PointMass(A));
		}
		[Skip]
		public static double LogEvidenceRatio([SkipIfUniform] Gaussian Product, [SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			return 0.0;
		}
	}

	/// <summary>
	/// Provides Gaussian evidence messages for <see cref="Factor.Product(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class GaussianProductEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="to_product">Outgoing message to 'product'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,a) p(product,a) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian product, Gaussian a, double b, [Fresh] Gaussian to_product)
		{
			//Gaussian to_product = GaussianProductOp.ProductAverageConditional(a, b);
			return to_product.GetLogAverageOf(product);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="to_product">Outgoing message to 'product'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,b) p(product,b) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian product, double a, Gaussian b, [Fresh] Gaussian to_product)
		{
			return LogAverageFactor(product, b, a, to_product);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double product, Gaussian a, double b)
		{
			Gaussian to_product = GaussianProductOp.ProductAverageConditional(a, b);
			return to_product.GetLogProb(product);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double product, double a, Gaussian b)
		{
			return LogAverageFactor(product, b, a);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double product, double a, double b)
		{
			return (product == Factor.Product(a, b)) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(product,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double product, double a, double b)
		{
			return LogAverageFactor(product, a, b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product) p(product) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian product, double a, double b)
		{
			return product.GetLogProb(Factor.Product(a, b));
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,a) p(product,a) factor(product,a,b) / sum_product p(product) messageTo(product))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian product, Gaussian a, double b) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,b) p(product,b) factor(product,a,b) / sum_product p(product) messageTo(product))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian product, double a, Gaussian b) { return LogEvidenceRatio(product, b, a); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(product,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double product, Gaussian a, double b) { return LogAverageFactor(product, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Constant value for 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(product,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double product, double a, Gaussian b) { return LogEvidenceRatio(product, b, a); }

	}

	/// <summary>
	/// Provides Gaussian evidence messages for <see cref="Factor.Ratio(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class GaussianRatioEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <param name="to_ratio">Outgoing message to 'ratio'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,a) p(ratio,a) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian ratio, Gaussian a, double b, [Fresh] Gaussian to_ratio)
		{
			//Gaussian to_ratio = GaussianProductOp.AAverageConditional(a, b);
			return to_ratio.GetLogAverageOf(ratio);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <param name="to_ratio">Outgoing message to 'ratio'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,b) p(ratio,b) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian ratio, double a, Gaussian b, [Fresh] Gaussian to_ratio)
		{
			return LogAverageFactor(ratio, b, a, to_ratio);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double ratio, Gaussian a, double b)
		{
			Gaussian to_ratio = GaussianProductOp.AAverageConditional(a, b);
			return to_ratio.GetLogProb(ratio);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double ratio, double a, Gaussian b)
		{
			return LogAverageFactor(ratio, b, a);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double ratio, double a, double b)
		{
			return (ratio == Factor.Ratio(a, b)) ? 0.0 : Double.NegativeInfinity;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(ratio,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double ratio, double a, double b)
		{
			return LogAverageFactor(ratio, a, b);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio) p(ratio) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian ratio, double a, double b)
		{
			return ratio.GetLogProb(Factor.Ratio(a, b));
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,a) p(ratio,a) factor(ratio,a,b) / sum_ratio p(ratio) messageTo(ratio))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian ratio, Gaussian a, double b) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,b) p(ratio,b) factor(ratio,a,b) / sum_ratio p(ratio) messageTo(ratio))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian ratio, double a, Gaussian b) { return LogEvidenceRatio(ratio, b, a); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a) p(a) factor(ratio,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double ratio, Gaussian a, double b) { return LogAverageFactor(ratio, a, b); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(b) p(b) factor(ratio,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double ratio, double a, Gaussian b) { return LogEvidenceRatio(ratio, b, a); }

	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Product(double, double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class GaussianProductVmpOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(product,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor(Gaussian product) { return 0.0; }

		internal const string NotSupportedMessage = "Variational Message Passing does not support a Product factor with fixed output and two random inputs.";

		/// <summary>
		/// VMP message to 'product'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a,b) p(a,b) factor(product,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gaussian ProductAverageLogarithm([SkipIfUniform] Gaussian A, [SkipIfUniform] Gaussian B)
		{
			Gaussian result = new Gaussian();
			// p(x|a,b) = N(E[a]*E[b], E[b]^2*var(a) + E[a]^2*var(b) + var(a)*var(b))
			double ma, va, mb, vb;
			A.GetMeanAndVariance(out ma, out va);
			B.GetMeanAndVariance(out mb, out vb);
			// Uses John Winn's rule for deterministic factors.
			// Strict variational inference would set the variance to 0.
			result.SetMeanAndVariance(ma * mb, mb * mb * va + ma * ma * vb + va * vb);
			return result;
		}

		public static Gaussian ProductDeriv(Gaussian Product, [SkipIfUniform, Stochastic] Gaussian A, [SkipIfUniform, Stochastic] Gaussian B, Gaussian to_A, Gaussian to_B)
		{
			if (A.IsPointMass) return ProductDeriv(Product, A.Point, B, to_B);
			if (B.IsPointMass) return ProductDeriv(Product, A, B.Point, to_A);
			double ma, va, mb, vb;
			A.GetMeanAndVariance(out ma, out va);
			B.GetMeanAndVariance(out mb, out vb);
			//Console.WriteLine("ma = {0}, va = {1}, mb = {2}, vb = {3}", ma, va, mb, vb);
			double ma0, va0, mb0, vb0;
			(A/to_A).GetMeanAndVariance(out ma0, out va0);
			(B/to_B).GetMeanAndVariance(out mb0, out vb0);
			Gaussian to_A2 = AAverageLogarithm(Product, B);
			double va2 = 1/(1/va0 + to_A2.Precision);
			double ma2 = va2*(ma0/va0 + to_A2.MeanTimesPrecision);
			Gaussian to_B2 = BAverageLogarithm(Product, A);
			double vb2 = 1/(1/vb0 + to_B2.Precision);
			double mb2 = vb2*(mb0/vb0 + to_B2.MeanTimesPrecision);
			double dva2 = 0;
			double dma2 = va2*mb + dva2*ma2/va2;
			double dvb2 = 0;
			// this doesn't seem to help
			//dvb2 = -vb2*vb2*2*ma2*dma2*Product.Precision;
			double dmb2 = vb2*ma + dvb2*mb2/vb2;
			double pPrec2 = 1/(va2*vb2 + va2*mb2*mb2 + vb2*ma2*ma2);
			double dpPrec2 = -(dva2*vb2 + va2*dvb2 + dva2*mb2*mb2 + va2*2*mb2*dmb2 + dvb2*ma2*ma2 + vb2*2*ma2*dma2)*pPrec2*pPrec2;
			double pMeanTimesPrec2 = ma2*mb2*pPrec2;
			double pMeanTimesPrec2Deriv = dma2*mb2*pPrec2 + ma2*dmb2*pPrec2 + ma2*mb2*dpPrec2;
			return Gaussian.FromNatural(pMeanTimesPrec2Deriv-1, 0);
		}

		[Skip]
		public static Gaussian ProductDeriv(Gaussian Product, [SkipIfUniform, Stochastic] Gaussian A, double B, Gaussian to_A)
		{
			return Gaussian.Uniform();
		}
		[Skip]
		public static Gaussian ProductDeriv(Gaussian Product, double A, [SkipIfUniform, Stochastic] Gaussian B, Gaussian to_B)
		{
			return ProductDeriv(Product, B, A, to_B);
		}

		/// <summary>
		/// VMP message to 'product'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(product,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gaussian ProductAverageLogarithm(double A, [SkipIfUniform] Gaussian B)
		{
			if (B.IsPointMass) return Gaussian.PointMass(A*B.Point);
			if (A == 0) return Gaussian.PointMass(A);
			// m = A*mb
			// v = A*A*vb
			// 1/v = (1/vb)/(A*A)
			// m/v = (mb/vb)/A
			return Gaussian.FromNatural(B.MeanTimesPrecision/A, B.Precision/(A*A));
		}
		/// <summary>
		/// VMP message to 'product'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'product' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'product' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(product,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Gaussian ProductAverageLogarithm([SkipIfUniform] Gaussian A, double B)
		{
			return ProductAverageLogarithm(B, A);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'product' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_product p(product) factor(product,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform] Gaussian Product, [Proper] Gaussian B)
		{
			if (B.IsPointMass) return AAverageLogarithm(Product, B.Point);
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B);
			double mb, vb;
			B.GetMeanAndVariance(out mb, out vb);
			// note this is exact if B is a point mass (vb=0).
			Gaussian result = new Gaussian();
			result.Precision = Product.Precision * (vb + mb*mb);
			result.MeanTimesPrecision = Product.MeanTimesPrecision * mb;
			return result;
		}


		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(product,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(GaussianProductVmpOp.NotSupportedMessage)]
		public static Gaussian AAverageLogarithm(double Product, [Proper] Gaussian B)
		{
			// Throw an exception rather than return a meaningless point mass.
			throw new NotSupportedException(GaussianProductVmpOp.NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'product' integrated out.
		/// The formula is <c>sum_product p(product) factor(product,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform] Gaussian Product, double B)
		{
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B);
			return GaussianProductOp.AAverageConditional(Product, B);
		}
		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageLogarithm(double Product, double B)
		{
			return GaussianProductOp.AAverageConditional(Product, B);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'product' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_product p(product) factor(product,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform] Gaussian Product, [Proper] Gaussian A)
		{
			return AAverageLogarithm(Product, A);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(product,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[NotSupported(GaussianProductVmpOp.NotSupportedMessage)]
		public static Gaussian BAverageLogarithm(double Product, [Proper] Gaussian A)
		{
			return AAverageLogarithm(Product, A);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'product' integrated out.
		/// The formula is <c>sum_product p(product) factor(product,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		public static Gaussian BAverageLogarithm([SkipIfUniform] Gaussian Product, double A)
		{
			return AAverageLogarithm(Product, A);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian BAverageLogarithm(double Product, double A)
		{
			return AAverageLogarithm(Product, A);
		}
	}

	/// <summary>
	/// Provides VMP Gaussian evidence messages for <see cref="Factor.Ratio(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Mature)]
	public static class GaussianRatioVmpOp
	{
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(ratio,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para><para>
		/// Variational Message Passing does not support a Ratio factor with fixed output or random denominator
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		internal const string NotSupportedMessage = "Variational Message Passing does not support a Ratio factor with fixed output or random denominator.";

		/// <summary>
		/// VMP message to 'ratio'
		/// </summary>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'ratio' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'ratio' as the random arguments are varied.
		/// The formula is <c>proj[sum_(b) p(b) factor(ratio,a,b)]</c>.
		/// </para><para>
		/// Variational Message Passing does not support a Ratio factor with fixed output or random denominator
		/// </para></remarks>
		[NotSupported(GaussianRatioVmpOp.NotSupportedMessage)]
		public static Gaussian RatioAverageLogarithm(Gaussian B)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'ratio'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'ratio' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'ratio' as the random arguments are varied.
		/// The formula is <c>proj[sum_(a) p(a) factor(ratio,a,b)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Gaussian RatioAverageLogarithm([SkipIfUniform] Gaussian A, double B)
		{
			return GaussianProductOp.AAverageConditional(A, B);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(ratio,a,b)))</c>.
		/// </para><para>
		/// Variational Message Passing does not support a Ratio factor with fixed output or random denominator
		/// </para></remarks>
		[NotSupported(GaussianRatioVmpOp.NotSupportedMessage)]
		public static Gaussian AAverageLogarithm(Gaussian B)
		{
			throw new NotSupportedException(NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'ratio' integrated out.
		/// The formula is <c>sum_ratio p(ratio) factor(ratio,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="ratio"/> is not a proper distribution</exception>
		public static Gaussian AAverageLogarithm([SkipIfUniform] Gaussian ratio, double B)
		{
			return GaussianProductOp.ProductAverageConditional(ratio, B);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="ratio">Constant value for 'ratio'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'a' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian AAverageLogarithm(double ratio, double B)
		{
			return GaussianProductOp.ProductAverageConditional(ratio, B);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para><para>
		/// Variational Message Passing does not support a Ratio factor with fixed output or random denominator
		/// </para></remarks>
		[NotSupported(GaussianRatioVmpOp.NotSupportedMessage)]
		public static Gaussian BAverageLogarithm()
		{
			throw new NotSupportedException(NotSupportedMessage);
		}
	}
}
