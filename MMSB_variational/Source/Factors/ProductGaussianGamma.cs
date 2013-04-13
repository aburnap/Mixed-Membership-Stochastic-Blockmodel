// (C) Copyright 2009-2010 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.Product(double, double)"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "Product", typeof(double), typeof(double))]
    [Quality(QualityBand.Preview)]
	public static class ProductGaussianGammaVmpOp
	{
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
		public static Gaussian ProductAverageLogarithm([SkipIfUniform] Gaussian A, [SkipIfUniform] Gamma B)
		{
			double ma = A.GetMean(), mb = B.GetMean();
			double va = A.GetVariance(), vb = B.GetVariance();
			if (Double.IsPositiveInfinity(va) || Double.IsPositiveInfinity(vb)) return Gaussian.Uniform();
			return Gaussian.FromMeanAndVariance(ma * mb, va*vb + va * mb*mb + vb * ma*ma);
		}
#if false
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
		public static Gaussian ProductAverageLogarithm(double A, [SkipIfUniform] Gamma B)
		{
			double mb, vb;
			B.GetMeanAndVariance(out mb, out vb);
			if (Double.IsPositiveInfinity(vb)) return Gaussian.Uniform();
			return Gaussian.FromMeanAndVariance(A * mb, A * A * vb);
		}
#endif

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
		public static Gaussian AAverageLogarithm([SkipIfUniform] Gaussian Product, [Proper] Gamma B)
		{
			if (B.IsPointMass) return GaussianProductVmpOp.AAverageLogarithm(Product, B.Point);
			if (Product.IsPointMass) return AAverageLogarithm(Product.Point, B);
			if (!B.IsProper()) throw new ImproperMessageException(B);
			double mb, vb;
			B.GetMeanAndVariance(out mb, out vb);
			// catch uniform case to avoid 0*Inf
			if (Product.IsUniform()) return Product;
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
		public static Gaussian AAverageLogarithm(double Product, [Proper] Gamma B)
		{
			// Throw an exception rather than return a meaningless point mass.
			throw new NotSupportedException(GaussianProductVmpOp.NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_B">Previous outgoing message to 'B'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'product' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_product p(product) factor(product,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gamma BAverageLogarithm([SkipIfUniform] Gaussian Product, [Proper] Gaussian A, [Proper] Gamma B, Gamma to_B)
		{
			if(B.IsPointMass) return Gamma.Uniform();
			if(Product.IsPointMass) return BAverageLogarithm(Product.Point, A);
			if (!B.IsProper()) throw new ImproperMessageException(B);
			// catch uniform case to avoid 0*Inf
			if (Product.IsUniform()) return Gamma.Uniform();
			double rateb = B.Rate;
			double shapeb = B.Shape;
			double ma, va;
			A.GetMeanAndVariance(out ma, out va);

			double c1 = ma*Product.MeanTimesPrecision;
			double c2 = -0.5*(ma*ma+va)*Product.Precision;
			double u = c2/(1 - shapeb*MMath.Trigamma(shapeb));
			double vb = shapeb/(rateb*rateb);
			double shape = 1 + vb*u;
			double rate = (u - 2*c2*(shapeb+1))/rateb - c1;
			//double rate2 = (1.0/(1 - shapeb*MMath.Trigamma(shapeb)) - 2*(shapeb+1))*c2/rateb - c1;
			//double rate3 = (ma*ma+va)*Product.Precision/rateb - ma*Product.MeanTimesPrecision;
			//if (rate < 0) throw new Exception("rate ("+rate+") < 0");
			Gamma msg = Gamma.FromShapeAndRate(shape, rate);
			double step = Rand.Double()*0.25;
			step = 1;
			if (step==1.0) return msg;
			return (msg^step) * (to_B^(1-step));
		}
#if false
		public static Gamma BAverageLogarithm([SkipIfUniform] Gaussian Product, [Proper] Gaussian A, [Proper] Gamma B, Gamma result)
		{
			if (B.IsPointMass) return Gamma.Uniform();
			if (Product.IsPointMass) return BAverageLogarithm(Product.Point, A);
			if (!B.IsProper()) throw new ImproperMessageException(B);
			// catch uniform case to avoid 0*Inf
			if (Product.IsUniform()) return Gamma.Uniform();
			double rateb = B.Rate;
			double shapeb = B.Shape;
			double ma, va;
			A.GetMeanAndVariance(out ma, out va);

			// encourage resulting (shape,rate) to be close to B
			// unfortunately, this doesn't force rate>0 - only scales result
			double penalty = 1;
			Vector Bvec = new Vector(2);
			Bvec[0] = result.Shape-1;
			Bvec[1] = result.Rate;
			Vector g = GradientVector(Product, A, B);
			Matrix W = WeightMatrix(B);
			W.SetToTranspose(W);
			PositiveDefiniteMatrix W2 = W.Outer();
			//W2.SetToSum(W2, PositiveDefiniteMatrix.IdentityScaledBy(W2.Rows, penalty));
			Vector Wg = W*g;
			Bvec = W2*Bvec;
			Bvec.SetToSum(1.0, Wg, penalty, Bvec);
			return Gamma.FromShapeAndRate(Bvec[0]/W2[0, 0] + 1, 0.0);

			Bvec.PredivideBy(W2);
			Bvec.Scale(1.0/(1+penalty));
			return Gamma.FromShapeAndRate(Bvec[0]+1, Bvec[1]);
		}
		public static Vector GradientVector(Gaussian Product, Gaussian A, Gamma B)
		{
			double rateb = B.Rate;
			double scaleb = 1.0/rateb;
			double shapeb = B.Shape;
			double ma, va;
			A.GetMeanAndVariance(out ma, out va);

			Matrix C = new Matrix(2, 2);
			C[0, 0] = scaleb;
			C[0, 1] = (2*shapeb+1)*scaleb*scaleb;
			C[1, 0] = -shapeb*scaleb*scaleb;
			C[1, 1] = -2*shapeb*(shapeb+1)*scaleb*scaleb*scaleb;
			Vector v = new Vector(2);
			v[0] = ma*Product.MeanTimesPrecision;
			v[1] = -0.5*(ma*ma+va)*Product.Precision;
			return C*v;
		}
		public static Matrix WeightMatrix(Gamma B)
		{
			Matrix A = new Matrix(2, 2);
			A[0, 0] = MMath.Trigamma(B.Shape);
			A[0, 1] = -1.0/B.Rate;
			A[1, 0] = A[0, 1];
			A[1, 1] = B.Shape/(B.Rate*B.Rate);
			return A;
		}
#endif
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Constant value for 'product'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(product,a,b)))</c>.
		/// </para></remarks>
		[NotSupported(GaussianProductVmpOp.NotSupportedMessage)]
		public static Gamma BAverageLogarithm(double Product, Gaussian A)
		{
			// Throw an exception rather than return a meaningless point mass.
			throw new NotSupportedException(GaussianProductVmpOp.NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="Product">Incoming message from 'product'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="to_B">Previous outgoing message to 'B'.</param>
		/// <returns>The outgoing VMP message to the 'b' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'product' integrated out.
		/// The formula is <c>sum_product p(product) factor(product,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="Product"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gamma BAverageLogarithm([SkipIfUniform] Gaussian Product, double A, [Proper] Gamma B, Gamma to_B)
		{
			return BAverageLogarithm(Product, Gaussian.PointMass(A), B, to_B);
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
			return GammaProductOp.BAverageConditional(Product, A);
		}
	}

}