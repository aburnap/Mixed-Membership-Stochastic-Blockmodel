// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using System.Linq;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using GaussianArray = MicrosoftResearch.Infer.Distributions.DistributionStructArray<MicrosoftResearch.Infer.Distributions.Gaussian, double>;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.InnerProductPartialCovariance"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "InnerProductPartialCovariance")]
	[Buffers("MeanOfB", "CovarianceOfB"/*, "Ebbt"*/)]
	[Quality(QualityBand.Preview)]
	public static class InnerProductPartialCovarianceOp
	{
		//-- VMP -------------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(X,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// Update the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Skip]
		public static PositiveDefiniteMatrix EbbtInit([IgnoreDependency] VectorGaussian B)
		{
			return new PositiveDefiniteMatrix(B.Dimension, B.Dimension);
		}

		/// <summary>
		/// Update the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix Ebbt(PositiveDefiniteMatrix CovarianceOfB, Vector MeanOfB, PositiveDefiniteMatrix result)
		{
			result.SetToSumWithOuter(CovarianceOfB, 1, MeanOfB, MeanOfB);
			return result;
		}

		///<summary>
		///Update the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Skip]
		public static PositiveDefiniteMatrix EaatInit([IgnoreDependency] GaussianArray A)
		{
			return new PositiveDefiniteMatrix(A.Count, A.Count);
		}

		/// <summary>
		/// Update the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix Eaat(GaussianArray A, PositiveDefiniteMatrix result)
		{
			int inner = A.Count;
			var VarianceOfA = Vector.Zero(inner);
			var MeanOfA = Vector.Zero(inner);
			for (int k = 0; k < inner; k++) {
				MeanOfA[k] = A[k].GetMean();
				VarianceOfA[k] = A[k].GetVariance();
			}
			result.SetToDiagonal(VarianceOfA);
			result.SetToSumWithOuter(result, 1, MeanOfA, MeanOfA);
			return result;
		}

		/// <summary>
		/// Initialise the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[Skip]
		public static PositiveDefiniteMatrix CovarianceOfBInit([IgnoreDependency] VectorGaussian B)
		{
			return new PositiveDefiniteMatrix(B.Dimension, B.Dimension);
		}

		/// <summary>
		/// Update the buffer 'CovarianceOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Storage for result.</param>
		/// <returns>New value of buffer 'CovarianceOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static PositiveDefiniteMatrix CovarianceOfB([Proper] VectorGaussian B, PositiveDefiniteMatrix result)
		{
			return B.GetVariance(result);
		}

		/// <summary>
		/// Initialise the buffer 'MeanOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'..</param>
		/// <returns>New value of buffer 'MeanOfB'</returns>
		[Skip]
		public static Vector MeanOfBInit([IgnoreDependency] VectorGaussian B)
		{
			return Vector.Zero(B.Dimension);
		}

		/// <summary>
		/// Update the buffer 'MeanOfB'
		/// </summary>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="CovarianceOfB">Buffer 'CovarianceOfB'.</param>
		/// <param name="result">Storage for result</param>
		/// <returns>New value of buffer 'MeanOfB'</returns>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Vector MeanOfB([Proper] VectorGaussian B, PositiveDefiniteMatrix CovarianceOfB, Vector result)
		{
			return B.GetMean(result, CovarianceOfB);
		}

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="X">Incoming message from 'X'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="MeanOfB">Buffer 'MeanOfB'.</param>
		/// <param name="CovarianceOfB">Buffer 'CovarianceOfB'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'A'.
		/// Because the factor is deterministic, 'X' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(B) p(B) log(sum_X p(X) factor(X,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="X"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray AAverageLogarithm([SkipIfUniform] Gaussian X, GaussianArray A, [SkipIfUniform] VectorGaussian B, Vector MeanOfB, PositiveDefiniteMatrix CovarianceOfB, /*PositiveDefiniteMatrix Ebbt,*/ GaussianArray result)
		{
			int inner = MeanOfB.Count;
			if (result == null) result = new GaussianArray(inner);
			// E[log N(x[i,j]; a[i,:]*b[:,j], 0)] = -0.5 E[(x[i,j]- sum_k a[i,k]*b[k,j])^2]/0 
			// = -0.5 (E[x[i,j]^2] - 2 E[x[i,j]] a[i,k] E[b[k,j]] + a[i,k] a[i,k2] E(b[k,j] b[k2,j]))/0
			// a[i,k] * (-2 E[x[i,j]] E[b[k,j]] + sum_{k2 not k} E[a[i,k2]] E(b[k,j] b[k2,j]))
			// a[i,k]^2 * E(b[k,j]^2)
			// message to a[i,k] = N(a; inv(prec[i,k])*(sum_j E[b[k,j]]*res[i,j,k]/var(x[i,j])), inv(prec[i,k]))
			// where res[i,j,k] = E[x[i,j]] - sum_{k2 not k} E[a[i,k2]] E[b[k2,j]]
			// prec[i,k] = sum_j E(b[k,j]^2)/var(x[i,j])
			// result.Precision = prec[i,k]
			// result.MeanTimesPrecision = sum_j E[b[k,j]]*res[i,j,k]/var(x[i,j]) 
			//                           = sum_j E[b[k,j]]*(X.MeanTimesPrecision - X.precision*(sum_{k2 not k}))

			var Ebbt = new PositiveDefiniteMatrix(inner, inner);
			// should we be caching this too? 
			Ebbt.SetToSumWithOuter(CovarianceOfB, 1, MeanOfB, MeanOfB);
			//var ma = A.Select(z => z.GetMean()).ToArray(); 
			var ma = new double[inner];
			for (int k = 0; k < inner; k++)
				ma[k] = A[k].GetMean();
			for (int k = 0; k < inner; k++) {
				double prec = 0.0;
				double pm = 0.0;
				prec += Ebbt[k, k] * X.Precision;
				double sum = 0.0;
				for (int r = 0; r < inner; r++)
					if (r != k)
						sum += Ebbt[r, k] * ma[r];
				pm += MeanOfB[k] * X.MeanTimesPrecision - X.Precision * sum;

				Gaussian rk = result[k];
				rk.Precision = prec;
				if (prec < 0)
					throw new ApplicationException("improper message");
				rk.MeanTimesPrecision = pm;
				result[k] = rk;
			}

			return result;
		}

		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="X">Incoming message from 'X'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'B'.
		/// Because the factor is deterministic, 'X' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(A) p(A) log(sum_X p(X) factor(X,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="X"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static VectorGaussian BAverageLogarithm([SkipIfUniform, Proper] Gaussian X, [SkipIfAllUniform] GaussianArray A, VectorGaussian result)
		{
			int K = A.Count;
			var va = Vector.Zero(K);
			var ma = Vector.Zero(K);
			for (int k = 0; k < K; k++) {
				double m, v;
				A[k].GetMeanAndVariance(out m, out v);
				ma[k] = m;
				va[k] = v;
			}

			result.Precision.SetToDiagonal(X.Precision * va);
			result.Precision.SetToSumWithOuter(result.Precision, X.Precision, ma, ma);
			result.MeanTimesPrecision.SetToProduct(ma, X.MeanTimesPrecision);
			//if (!result.Precision.IsPositiveDefinite())
			//    throw new ApplicationException("improper message");

			return result;
		}

		/// <summary>
		/// VMP message to 'X'
		/// </summary>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If all elements are uniform, the result will be uniform.</param>
		/// <param name="MeanOfB">Buffer 'MeanOfB'.</param>
		/// <param name="CovarianceOfB">Buffer 'CovarianceOfB'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'X' as the random arguments are varied.
		/// The formula is <c>proj[sum_(A,B) p(A,B) factor(X,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static Gaussian XAverageLogarithm([SkipIfAllUniform] GaussianArray A, [SkipIfAllUniform] VectorGaussian B, Vector MeanOfB, PositiveDefiniteMatrix CovarianceOfB)
		{
			int K = MeanOfB.Count;
			// p(x|a,b) = N(E[a]'*E[b], E[b]'*var(a)*E[b] + E[a]'*var(b)*E[a] + trace(var(a)*var(b)))
			var ma = Vector.Zero(K);
			var va = Vector.Zero(K);
			for (int k = 0; k < K; k++) {
				double m, v;
				A[k].GetMeanAndVariance(out m, out v);
				ma[k] = m;
				va[k] = v;
			}
			// Uses John Winn's rule for deterministic factors.
			// Strict variational inference would set the variance to 0.
			var mbj2 = Vector.Zero(K);
			mbj2.SetToFunction(MeanOfB, x => x * x);
			// slooow
			Gaussian result = new Gaussian();
			result.SetMeanAndVariance(ma.Inner(MeanOfB), va.Inner(mbj2) + CovarianceOfB.QuadraticForm(ma) + va.Inner(CovarianceOfB.Diagonal()));
			if (result.Precision < 0)
				throw new ApplicationException("improper message");

			return result;
		}

		[Skip]
		public static Gaussian XAverageLogarithmInit()
		{
			return new Gaussian();
		}
	}
}
