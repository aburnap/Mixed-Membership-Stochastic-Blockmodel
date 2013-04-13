// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.InnerProduct"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Vector), "InnerProduct")]
	[FactorMethod(typeof(Factor), "InnerProduct", typeof(Vector), typeof(Vector))]
	[Buffers("AVariance", "BVariance", "AMean", "BMean")]
	[Quality(QualityBand.Mature)]
	public static class InnerProductOp
	{
		//-- VMP ---------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(innerProduct,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		const string NotSupportedMessage = "Variational Message Passing does not support an InnerProduct factor with fixed output.";

		/// <summary>
		/// VMP message to 'innerProduct'
		/// </summary>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <returns>The outgoing VMP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'innerProduct' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian InnerProductAverageLogarithm([Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance)
		{
			Gaussian result = new Gaussian();
			// p(x|a,b) = N(E[a]'*E[b], E[b]'*var(a)*E[b] + E[a]'*var(b)*E[a] + trace(var(a)*var(b)))
			// Uses John Winn's rule for deterministic factors.
			// Strict variational inference would set the variance to 0.
			result.SetMeanAndVariance(AMean.Inner(BMean), AVariance.QuadraticForm(BMean) + BVariance.QuadraticForm(AMean) + AVariance.Inner(BVariance));
			return result;
		}

		/// <summary>
		/// VMP message to 'innerProduct'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <returns>The outgoing VMP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'innerProduct' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian InnerProductAverageLogarithm(Vector A, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance)
		{
			Gaussian result = new Gaussian();
			// Uses John Winn's rule for deterministic factors.
			// Strict variational inference would set the variance to 0.
			// p(x) = N(a' E[b], a' var(b) a)
			result.SetMeanAndVariance(A.Inner(BMean), BVariance.QuadraticForm(A));
			return result;
		}
		[Skip]
		public static Gaussian InnerProductAverageLogarithmInit()
		{
			return new Gaussian();
		}
		/// <summary>
		/// VMP message to 'innerProduct'
		/// </summary>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing VMP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'innerProduct' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian InnerProductAverageLogarithm([Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, Vector B)
		{
			return InnerProductAverageLogarithm(B, AMean, AVariance);
		}

		/// <summary>
		/// Initialise the buffer 'BVariance'
		/// </summary>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Initial value of buffer 'BVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix BVarianceInit([IgnoreDependency] VectorGaussian B)
		{
			return new PositiveDefiniteMatrix(B.Dimension, B.Dimension);
		}
		/// <summary>
		/// Update the buffer 'BVariance'
		/// </summary>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix BVariance([Proper] VectorGaussian B, PositiveDefiniteMatrix result)
		{
			return B.GetVariance(result);
		}

		/// <summary>
		/// Initialise the buffer 'BMean'
		/// </summary>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Initial value of buffer 'BMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector BMeanInit([IgnoreDependency] VectorGaussian B)
		{
			return Vector.Zero(B.Dimension);
		}
		/// <summary>
		/// Update the buffer 'BMean'
		/// </summary>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Vector BMean([Proper] VectorGaussian B, [Fresh] PositiveDefiniteMatrix BVariance, Vector result)
		{
			return B.GetMean(result, BVariance);
		}

		/// <summary>
		/// Initialise the buffer 'AVariance'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>Initial value of buffer 'AVariance'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static PositiveDefiniteMatrix AVarianceInit([IgnoreDependency] VectorGaussian A)
		{
			return new PositiveDefiniteMatrix(A.Dimension, A.Dimension);
		}
		/// <summary>
		/// Update the buffer 'AVariance'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix AVariance([Proper] VectorGaussian A, PositiveDefiniteMatrix result)
		{
			return A.GetVariance(result);
		}

		/// <summary>
		/// Initialise the buffer 'AMean'
		/// </summary>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <returns>Initial value of buffer 'AMean'</returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		[Skip]
		public static Vector AMeanInit([IgnoreDependency] VectorGaussian A)
		{
			return Vector.Zero(A.Dimension);
		}
		/// <summary>
		/// Update the buffer 'AMean'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// 
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static Vector AMean([Proper] VectorGaussian A, [Fresh] PositiveDefiniteMatrix AVariance, Vector result)
		{
			return A.GetMean(result, AVariance);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// Because the factor is deterministic, 'innerProduct' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(b) p(b) log(sum_innerProduct p(innerProduct) factor(innerProduct,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static VectorGaussian AAverageLogarithm([SkipIfUniform] Gaussian innerProduct, [SkipIfUniform] VectorGaussian B, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance, VectorGaussian result)
		{
			if (innerProduct.IsPointMass) return AAverageLogarithm(innerProduct.Point, B, result);
			if (result == default(VectorGaussian)) result = new VectorGaussian(B.Dimension);
			// E[log N(x; ab, 0)] = -0.5 E[(x-ab)^2]/0 = -0.5 (E[x^2] - 2 E[x] a' E[b] + trace(aa' E[bb']))/0
			// message to a = N(a; E[x]*inv(var(b)+E[b]E[b]')*E[b], var(x)*inv(var(b)+E[b]E[b]'))
			// result.Precision = (var(b)+E[b]*E[b]')/var(x)
			// result.MeanTimesPrecision = E[x]/var(x)*E[b] = E[b]*X.MeanTimesPrecision
			// note this is exact if B is a point mass (vb=0).
			result.Precision.SetToSumWithOuter(BVariance, 1, BMean, BMean);
			result.Precision.SetToProduct(result.Precision, innerProduct.Precision);
			result.MeanTimesPrecision.SetToProduct(BMean, innerProduct.MeanTimesPrecision);
			return result;
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'a'.
		/// The formula is <c>exp(sum_(b) p(b) log(factor(innerProduct,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(InnerProductOp.NotSupportedMessage)]
		public static VectorGaussian AAverageLogarithm(double innerProduct, [SkipIfUniform] VectorGaussian B, VectorGaussian result)
		{
			throw new NotSupportedException(InnerProductOp.NotSupportedMessage);
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' with 'innerProduct' integrated out.
		/// The formula is <c>sum_innerProduct p(innerProduct) factor(innerProduct,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		public static VectorGaussian AAverageLogarithm([SkipIfUniform] Gaussian innerProduct, Vector B, VectorGaussian result)
		{
			if (innerProduct.IsPointMass) return AAverageLogarithm(innerProduct.Point, B, result);
			if (result == default(VectorGaussian)) result = new VectorGaussian(B.Count);
			result.Precision.SetToOuter(B, B);
			result.Precision.Scale(innerProduct.Precision);
			result.MeanTimesPrecision.SetToProduct(B, innerProduct.MeanTimesPrecision);
			return result;
		}

		/// <summary>
		/// VMP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(InnerProductOp.NotSupportedMessage)]
		public static VectorGaussian AAverageLogarithm(double innerProduct, Vector B, VectorGaussian result)
		{
			// This case could be supported if we had low-rank VectorGaussian distributions.
			throw new NotSupportedException(InnerProductOp.NotSupportedMessage);
			if (result == default(VectorGaussian)) result = new VectorGaussian(B.Count);
			result.Point = result.Point;
			result.Point.SetToProduct(B, innerProduct / B.Inner(B));
			return result;
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// Because the factor is deterministic, 'innerProduct' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(a) p(a) log(sum_innerProduct p(innerProduct) factor(innerProduct,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static VectorGaussian BAverageLogarithm([SkipIfUniform] Gaussian innerProduct, [SkipIfUniform] VectorGaussian A, [Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, VectorGaussian result)
		{
			return AAverageLogarithm(innerProduct, A, AMean, AVariance, result);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'b'.
		/// The formula is <c>exp(sum_(a) p(a) log(factor(innerProduct,a,b)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[NotSupported(InnerProductOp.NotSupportedMessage)]
		public static VectorGaussian BAverageLogarithm(double innerProduct, [SkipIfUniform] VectorGaussian A, VectorGaussian result)
		{
			throw new NotSupportedException(InnerProductOp.NotSupportedMessage);
			return AAverageLogarithm(innerProduct, A, result);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' with 'innerProduct' integrated out.
		/// The formula is <c>sum_innerProduct p(innerProduct) factor(innerProduct,a,b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		public static VectorGaussian BAverageLogarithm([SkipIfUniform] Gaussian innerProduct, Vector A, VectorGaussian result)
		{
			return AAverageLogarithm(innerProduct, A, result);
		}

		/// <summary>
		/// VMP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(InnerProductOp.NotSupportedMessage)]
		public static VectorGaussian BAverageLogarithm(double innerProduct, Vector A, VectorGaussian result)
		{
			throw new NotSupportedException(InnerProductOp.NotSupportedMessage);
			return AAverageLogarithm(innerProduct, A, result);
		}

		// ----------------------- AverageConditional ------------------------------

		const string LowRankNotSupportedMessage = "A InnerProduct factor with fixed output is not yet implemented for Expectation Propagation.";
		const string BothRandomNotSupportedMessage = "An InnerProduct factor between two VectorGaussian variables is not yet implemented for Expectation Propagation.  Try using Variational Message Passing.";

		/// <summary>
		/// EP message to 'innerProduct'
		/// </summary>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>The outgoing EP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'innerProduct' as the random arguments are varied.
		/// The formula is <c>proj[p(innerProduct) sum_(a,b) p(a,b) factor(innerProduct,a,b)]/p(innerProduct)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static Gaussian InnerProductAverageConditional([SkipIfUniform] VectorGaussian A, [SkipIfUniform] VectorGaussian B)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(innerProduct,b) p(innerProduct,b) factor(innerProduct,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static VectorGaussian AAverageConditional(Gaussian innerProduct, VectorGaussian A, [SkipIfUniform] VectorGaussian B, VectorGaussian result)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(innerProduct,a) p(innerProduct,a) factor(innerProduct,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static VectorGaussian BAverageConditional(Gaussian innerProduct, [SkipIfUniform] VectorGaussian A, VectorGaussian B, VectorGaussian result)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'innerProduct'
		/// </summary>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <returns>The outgoing EP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'innerProduct' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian InnerProductAverageConditional([SkipIfUniform] Vector A, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance)
		{
			return InnerProductAverageLogarithm(A, BMean, BVariance);
		}
		[Skip]
		public static Gaussian InnerProductAverageConditionalInit()
		{
			return new Gaussian();
		}

		/// <summary>
		/// EP message to 'innerProduct'
		/// </summary>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>The outgoing EP message to the 'innerProduct' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'innerProduct' conditioned on the given values.
		/// </para></remarks>
		public static Gaussian InnerProductAverageConditional([Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, Vector B)
		{
			return InnerProductAverageConditional(B, AMean, AVariance);
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'a' as the random arguments are varied.
		/// The formula is <c>proj[p(a) sum_(innerProduct) p(innerProduct) factor(innerProduct,a,b)]/p(a)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		public static VectorGaussian AAverageConditional([SkipIfUniform] Gaussian innerProduct, Vector B, VectorGaussian result)
		{
			if (innerProduct.IsPointMass) return AAverageConditional(innerProduct.Point, B, result);
			if (result == default(VectorGaussian)) result = new VectorGaussian(B.Count);
			// (m - a'b)^2/v = (a'bb'a - 2a'bm + m^2)/v
			result.Precision.SetToOuter(B, B);
			result.Precision.Scale(innerProduct.Precision);
			result.MeanTimesPrecision.SetToProduct(B, innerProduct.MeanTimesPrecision);
			return result;
		}

		/// <summary>
		/// EP message to 'a'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'a' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(InnerProductOp.LowRankNotSupportedMessage)]
		public static VectorGaussian AAverageConditional(double innerProduct, Vector B, VectorGaussian result)
		{
			// a'*b == ip, therefore:
			// E[a]'*b == ip
			// b'*var(a)*b == 0
			// inv(var(a)) = Inf*bb'
			// E[a] = ip*b/(b'b)
			throw new NotImplementedException(LowRankNotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'b' as the random arguments are varied.
		/// The formula is <c>proj[p(b) sum_(innerProduct) p(innerProduct) factor(innerProduct,a,b)]/p(b)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="innerProduct"/> is not a proper distribution</exception>
		public static VectorGaussian BAverageConditional([SkipIfUniform] Gaussian innerProduct, Vector A, VectorGaussian result)
		{
			return AAverageConditional(innerProduct, A, result);
		}

		/// <summary>
		/// EP message to 'b'
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'b' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(InnerProductOp.LowRankNotSupportedMessage)]
		public static VectorGaussian BAverageConditional(double innerProduct, Vector A, VectorGaussian result)
		{
			throw new NotImplementedException(LowRankNotSupportedMessage);
			return AAverageConditional(innerProduct, A, result);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(innerProduct,a,b) p(innerProduct,a,b) factor(innerProduct,a,b) / sum_innerProduct p(innerProduct) messageTo(innerProduct))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static double LogEvidenceRatio(Gaussian innerProduct, VectorGaussian A, VectorGaussian B)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(innerProduct,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static double LogEvidenceRatio(double innerProduct, VectorGaussian A, VectorGaussian B)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="b">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(innerProduct,b) p(innerProduct,b) factor(innerProduct,a,b) / sum_innerProduct p(innerProduct) messageTo(innerProduct))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian innerProduct, Vector A, VectorGaussian b)
		{
			return 0.0;
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(innerProduct,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double innerProduct, Vector A, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance)
		{
			return LogAverageFactor(innerProduct, A, BMean, BVariance);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="a">Incoming message from 'a'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(innerProduct,a) p(innerProduct,a) factor(innerProduct,a,b) / sum_innerProduct p(innerProduct) messageTo(innerProduct))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian innerProduct, VectorGaussian a, Vector B)
		{
			return LogEvidenceRatio(innerProduct, B, a);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(innerProduct,a,b))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double innerProduct, [Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, Vector B)
		{
			return LogEvidenceRatio(innerProduct, B, AMean, AVariance);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(innerProduct,a,b) p(innerProduct,a,b) factor(innerProduct,a,b))</c>.
		/// </para></remarks>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static double LogAverageFactor(Gaussian innerProduct, VectorGaussian A, VectorGaussian B)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
			Gaussian to_innerProduct = InnerProductAverageConditional(A, B);
			return to_innerProduct.GetLogAverageOf(innerProduct);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Incoming message from 'a'.</param>
		/// <param name="B">Incoming message from 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(a,b) p(a,b) factor(innerProduct,a,b))</c>.
		/// </para></remarks>
		[NotSupported(InnerProductOp.BothRandomNotSupportedMessage)]
		public static double LogAverageFactor(double innerProduct, VectorGaussian A, VectorGaussian B)
		{
			throw new NotSupportedException(InnerProductOp.BothRandomNotSupportedMessage);
			Gaussian to_innerProduct = InnerProductAverageConditional(A, B);
			return to_innerProduct.GetLogProb(innerProduct);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Incoming message from 'innerProduct'.</param>
		/// <param name="to_innerProduct">Outgoing message to 'innerProduct'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(innerProduct) p(innerProduct) factor(innerProduct,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(Gaussian innerProduct, [Fresh] Gaussian to_innerProduct)
		{
			return to_innerProduct.GetLogAverageOf(innerProduct);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="A">Constant value for 'a'.</param>
		/// <param name="BMean">Buffer 'BMean'.</param>
		/// <param name="BVariance">Buffer 'BVariance'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(innerProduct,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double innerProduct, Vector A, [Fresh] Vector BMean, [Fresh] PositiveDefiniteMatrix BVariance)
		{
			Gaussian to_innerProduct = InnerProductAverageConditional(A, BMean, BVariance);
			return to_innerProduct.GetLogProb(innerProduct);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="innerProduct">Constant value for 'innerProduct'.</param>
		/// <param name="AMean">Buffer 'AMean'.</param>
		/// <param name="AVariance">Buffer 'AVariance'.</param>
		/// <param name="B">Constant value for 'b'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(innerProduct,a,b))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double innerProduct, [Fresh] Vector AMean, [Fresh] PositiveDefiniteMatrix AVariance, Vector B)
		{
			return LogAverageFactor(innerProduct, B, AMean, AVariance);
		}

	}
}
