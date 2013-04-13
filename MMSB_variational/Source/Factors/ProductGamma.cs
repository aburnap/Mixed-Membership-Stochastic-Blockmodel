// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing Gamma messages for <see cref="Factor.Product(double, double)"/> and <see cref="Factor.Ratio"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[FactorMethod(new string[] { "A", "Product", "B" }, typeof(Factor), "Ratio", typeof(double), typeof(double))]
  [Quality(QualityBand.Preview)]
	public static class GammaProductOp
	{
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
		public static Gamma ProductAverageConditional(double A, [SkipIfUniform] Gamma B)
		{
			return GammaProductVmpOp.ProductAverageLogarithm(A, B);
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
		public static Gamma ProductAverageConditional([SkipIfUniform] Gamma A, double B)
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
		public static Gamma ProductAverageConditional(double a, double b)
		{
			return Gamma.PointMass(a * b);
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
		public static Gamma AAverageConditional([SkipIfUniform] Gamma Product, double B)
		{
			if (Product.IsPointMass) return AAverageConditional(Product.Point, B);
			// (ab)^(shape-1) exp(-rate*(ab))
			return Gamma.FromShapeAndRate(Product.Shape, Product.Rate * B);
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
		public static Gamma AAverageConditional(double Product, double B)
		{
			Gamma result = new Gamma();
			if (B == 0) {
				if (Product != 0) throw new AllZeroException();
				result.SetToUniform();
			} else if((Product > 0) != (B > 0)) throw new AllZeroException("Product and argument do not have the same sign");
			else result.Point = Product / B;
			return result;
		}
#if false
		// This implementation allows GibbsRatio() to pass
		public static Gamma AAverageConditional(double Product, [SkipIfUniform] Gamma B)
		{
			throw new NotImplementedException();
		}
		public static Gamma BAverageConditional(double Product, [SkipIfUniform] Gamma A)
		{
			return AAverageConditional(Product, A);
		}
#else
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
		public static Gamma AAverageConditional(double Product, Gamma A, [SkipIfUniform] Gamma B)
		{
			// Z = int_x int_y delta(c - xy) Ga(x;ax,bx) Ga(y;ay,by) dx dy
			//   = int_x Ga(x;ax,bx) Ga(c/x;ay,by)/x dx
			//   = bx^ax by^ay /Gamma(ax)/Gamma(ay) 2 (bx/by)^(-(ax-ay)/2) BesselK(ax-ay, 2 sqrt(bx by))
			throw new NotImplementedException(); // BesselK is not implemented yet
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
		public static Gamma BAverageConditional(double Product, [SkipIfUniform] Gamma A, Gamma B)
		{
			return AAverageConditional(Product, B, A);
		}
#endif

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
		public static Gamma BAverageConditional([SkipIfUniform] Gamma Product, double A)
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
		public static Gamma BAverageConditional(double Product, double A)
		{
			return AAverageConditional(Product, A);
		}
	}

	/// <summary>
	/// Provides Gamma evidence messages for <see cref="Factor.Product(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class GammaProductEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,a) p(product,a) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma product, Gamma a, double b)
		{
			Gamma to_product = GammaProductOp.ProductAverageConditional(a, b);
			return to_product.GetLogAverageOf(product);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="product">Incoming message from 'product'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(product,b) p(product,b) factor(product,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma product, double a, Gamma b)
		{
			return LogAverageFactor(product, b, a);
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
		public static double LogAverageFactor(double product, Gamma a, double b)
		{
			Gamma to_product = GammaProductOp.ProductAverageConditional(a, b);
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
		public static double LogAverageFactor(double product, double a, Gamma b)
		{
			return LogAverageFactor(product, b, a);
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
		public static double LogEvidenceRatio(Gamma product, Gamma a, double b) { return 0.0; }
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
		public static double LogEvidenceRatio(Gamma product, double a, Gamma b) { return LogEvidenceRatio(product, b, a); }
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
		public static double LogEvidenceRatio(double product, Gamma a, double b) { return LogAverageFactor(product, a, b); }
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
		public static double LogEvidenceRatio(double product, double a, Gamma b) { return LogEvidenceRatio(product, b, a); }
	}

	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Product(double, double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class GammaProductVmpOp
	{
		[Skip]
		public static double AverageLogFactor(Gamma product) { return 0.0; }

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
		public static Gamma ProductAverageLogarithm([SkipIfUniform] Gamma A, [SkipIfUniform] Gamma B)
		{
			// E[x] = E[a]*E[b]
			// E[log(x)] = E[log(a)]+E[log(b)]
			return Gamma.FromMeanAndMeanLog(A.GetMean() * B.GetMean(), A.GetMeanLog() + B.GetMeanLog());
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
		public static Gamma ProductAverageLogarithm(double A, [SkipIfUniform] Gamma B)
		{
			if (B.IsPointMass) return Gamma.PointMass(A * B.Point);
			return Gamma.FromShapeAndRate(B.Shape, B.Rate/A);
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
		public static Gamma ProductAverageLogarithm([SkipIfUniform] Gamma A, double B)
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
		public static Gamma AAverageLogarithm([SkipIfUniform] Gamma Product, [Proper] Gamma B)
		{
			if (B.IsPointMass) return AAverageLogarithm(Product, B.Point);
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B);
			// (ab)^(shape-1) exp(-rate*(ab))
			return Gamma.FromShapeAndRate(Product.Shape, Product.Rate * B.GetMean());
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
		public static Gamma AAverageLogarithm(double Product, [Proper] Gamma B)
		{
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
		public static Gamma AAverageLogarithm([SkipIfUniform] Gamma Product, double B)
		{
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B);
			return GammaProductOp.AAverageConditional(Product, B);
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
		public static Gamma AAverageLogarithm(double Product, double B)
		{
			return Gamma.PointMass(Product/B);
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
		public static Gamma BAverageLogarithm([SkipIfUniform] Gamma Product, [Proper] Gamma A)
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
		public static Gamma BAverageLogarithm(double Product, [Proper] Gamma A)
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
		public static Gamma BAverageLogarithm([SkipIfUniform] Gamma Product, double A)
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
		public static Gamma BAverageLogarithm(double Product, double A)
		{
			return AAverageLogarithm(Product, A);
		}
	}

	/// <summary>
	/// Provides Gamma evidence messages for <see cref="Factor.Ratio(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class GammaRatioEvidenceOp
	{
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="b">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,a) p(ratio,a) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma ratio, Gamma a, double b)
		{
			Gamma to_ratio = GammaProductOp.AAverageConditional(a, b);
			return to_ratio.GetLogAverageOf(ratio);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="ratio">Incoming message from 'ratio'.</param>
		/// <param name="a">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(ratio,b) p(ratio,b) factor(ratio,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gamma ratio, double a, Gamma b)
		{
			return LogAverageFactor(ratio, b, a);
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
		public static double LogAverageFactor(double ratio, Gamma a, double b)
		{
			Gamma to_ratio = GammaProductOp.AAverageConditional(a, b);
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
		public static double LogAverageFactor(double ratio, double a, Gamma b)
		{
			return LogAverageFactor(ratio, b, a);
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
		public static double LogEvidenceRatio(Gamma ratio, Gamma a, double b) { return 0.0; }
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
		public static double LogEvidenceRatio(Gamma ratio, double a, Gamma b) { return LogEvidenceRatio(ratio, b, a); }
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
		public static double LogEvidenceRatio(double ratio, Gamma a, double b) { return LogAverageFactor(ratio, a, b); }
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
		public static double LogEvidenceRatio(double ratio, double a, Gamma b) { return LogEvidenceRatio(ratio, b, a); }

	}

	/// <summary>
	/// Provides VMP Gamma evidence messages for <see cref="Factor.Ratio(double,double)"/>,
	/// given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Ratio", typeof(double), typeof(double))]
	[Quality(QualityBand.Preview)]
	public static class GammaRatioVmpOp
	{
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
		public static Gamma RatioAverageLogarithm(Gamma B)
		{
			throw new NotSupportedException(GaussianRatioVmpOp.NotSupportedMessage);
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
		public static Gamma RatioAverageLogarithm([SkipIfUniform] Gamma A, double B)
		{
			return GammaProductOp.AAverageConditional(A, B);
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
		public static Gamma AAverageLogarithm(Gamma B)
		{
			throw new NotSupportedException(GaussianRatioVmpOp.NotSupportedMessage);
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
		public static Gamma AAverageLogarithm([SkipIfUniform] Gamma ratio, double B)
		{
			return GammaProductOp.ProductAverageConditional(ratio, B);
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
		public static Gamma AAverageLogarithm(double ratio, double B)
		{
			return GammaProductOp.ProductAverageConditional(ratio, B);
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
		public static Gamma BAverageLogarithm()
		{
			throw new NotSupportedException(GaussianRatioVmpOp.NotSupportedMessage);
		}
	}
}
