// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using GaussianArray2D = MicrosoftResearch.Infer.Distributions.DistributionStructArray2D<MicrosoftResearch.Infer.Distributions.Gaussian, double>;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.MatrixMultiply"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "MatrixMultiply")]
	[Quality(QualityBand.Preview)]
	public static class MatrixMultiplyOp
	{
		//-- VMP -------------------------------------------------------------------------------------------------------------
		/// <summary>
		/// Evidence message for VMP
		/// </summary>
		/// <returns>Zero</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(factor(matrixMultiply,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for VMP.
		/// </para></remarks>
		[Skip]
		public static double AverageLogFactor() { return 0.0; }

		/// <summary>
		/// VMP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[sum_(A,B) p(A,B) factor(matrixMultiply,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D MatrixMultiplyAverageLogarithm([SkipIfUniform] DistributionArray2D<Gaussian> A, [SkipIfUniform] DistributionArray2D<Gaussian> B, GaussianArray2D result)
		{
			if (result == null) result = new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), B.GetLength(1));
			// x[i,j] = sum_k a[i,k]*b[k,j]
			// E(x[i,j]) = sum_k E(a[i,k])*E(b[k,j])
			// var(x[i,j]) = sum_k E(a[i,k])^2*var(b[k,j]) + var(a[i,k])*E(b[k,j])^2 + var(a[i,k])*var(b[k,j])
			int rows = result.GetLength(0);
			int cols = result.GetLength(1);
			int inner = A.GetLength(1);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					double mean = 0; double var = 0;
					for (int k = 0; k < inner; k++) {
						double Am, Av, Bm, Bv;
						A[i, k].GetMeanAndVariance(out Am, out Av);
						B[k, j].GetMeanAndVariance(out Bm, out Bv);
						mean += Am * Bm;
						var += Av * Bm * Bm + Bv * Am * Am + Av * Bv;
					}
					Gaussian rij = result[i, j];
					rij.SetMeanAndVariance(mean, var);
					result[i, j] = rij;
				}
			}
			return result;
		}
		[Skip]
		public static GaussianArray2D MatrixMultiplyAverageLogarithmInit([IgnoreDependency] DistributionArray2D<Gaussian> A, [IgnoreDependency] DistributionArray2D<Gaussian> B)
		{
			return new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), B.GetLength(1));
		}
		[Skip]
		public static GaussianArray2D MatrixMultiplyAverageLogarithmInit([IgnoreDependency] double[,] A, [IgnoreDependency] DistributionArray2D<Gaussian> B)
		{
			return new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), B.GetLength(1));
		}
		[Skip]
		public static GaussianArray2D MatrixMultiplyAverageLogarithmInit([IgnoreDependency] DistributionArray2D<Gaussian> A, [IgnoreDependency] double[,] B)
		{
			return new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), B.GetLength(1));
		}
		/// <summary>
		/// VMP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[sum_(B) p(B) factor(matrixMultiply,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D MatrixMultiplyAverageLogarithm(double[,] A, [SkipIfUniform] GaussianArray2D B, GaussianArray2D result)
		{
			return MatrixMultiplyAverageConditional(A, B, result);
		}
		/// <summary>
		/// VMP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[sum_(A) p(A) factor(matrixMultiply,A,B)]</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static GaussianArray2D MatrixMultiplyAverageLogarithm([SkipIfUniform] GaussianArray2D A, double[,] B, GaussianArray2D result)
		{
			return MatrixMultiplyAverageConditional(A, B, result);
		}

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_A">Previous outgoing message to 'A'.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'A'.
		/// Because the factor is deterministic, 'matrixMultiply' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(B) p(B) log(sum_matrixMultiply p(matrixMultiply) factor(matrixMultiply,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D AAverageLogarithm([SkipIfUniform] GaussianArray2D matrixMultiply, [Stochastic] GaussianArray2D A, [SkipIfUniform] GaussianArray2D B, GaussianArray2D to_A)
		{
			GaussianArray2D result = to_A;
			if (result == null) result = new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), A.GetLength(1));
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
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = A.GetLength(1);
			double[] ab = new double[cols];
			//Gaussian temp = new Gaussian();
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					double sum = 0.0;
					for (int k = 0; k < inner; k++) {
						sum += A[i, k].GetMean() * B[k, j].GetMean();
					}
					ab[j] = sum;
				}
				for (int k = 0; k < inner; k++) {
					double prec = 0.0;
					double pm = 0.0;
					double Am = A[i, k].GetMean();
					for (int j = 0; j < cols; j++) {
						double Bm, Bv;
						B[k, j].GetMeanAndVariance(out Bm, out Bv);
						Gaussian x = matrixMultiply[i, j];
						prec += (Bv + Bm * Bm) * x.Precision;
						//pm += Bm * (x.MeanTimesPrecision - x.Precision * (ab[j] - Am * Bm));
						//Replace previous line with:
						ab[j] -= Am * Bm;
						pm += Bm * (x.MeanTimesPrecision - x.Precision * ab[j]);
					}

					Gaussian old_result = (Gaussian)result[i, k].Clone();
					Gaussian rik = result[i, k];
					rik.Precision = prec;
					rik.MeanTimesPrecision = pm;
					result[i, k] = rik;
					Gaussian partial = A[i, k] / old_result;
					Gaussian newPosterior = partial * rik;
					Am = newPosterior.GetMean();
					for (int j = 0; j < cols; j++) ab[j] += Am * B[k, j].GetMean();
				}
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="to_A">Previous outgoing message to 'A'.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' with 'matrixMultiply' integrated out.
		/// The formula is <c>sum_matrixMultiply p(matrixMultiply) factor(matrixMultiply,A,B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static GaussianArray2D AAverageLogarithm([SkipIfUniform] GaussianArray2D matrixMultiply, [Proper, Stochastic] GaussianArray2D A, double[,] B, GaussianArray2D to_A)
		{
			GaussianArray2D result = to_A;
			if (result == null) result = new DistributionStructArray2D<Gaussian, double>(A.GetLength(0), A.GetLength(1));
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = A.GetLength(1);
			double[] ab = new double[cols];
			//Gaussian temp = new Gaussian();
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					double sum = 0.0;
					for (int k = 0; k < inner; k++) {
						sum += A[i, k].GetMean() * B[k, j];
					}
					ab[j] = sum;
				}
				for (int k = 0; k < inner; k++) {
					double prec = 0.0;
					double pm = 0.0;
					double Am = A[i, k].GetMean();
					for (int j = 0; j < cols; j++) {
						double Bm = B[k, j];
						Gaussian x = matrixMultiply[i, j];
						prec += (Bm * Bm) * x.Precision;
						//pm += Bm * (x.MeanTimesPrecision - x.Precision * (ab[j] - Am * Bm));
						//Replace previous line with:
						ab[j] -= Am * Bm;
						pm += Bm * (x.MeanTimesPrecision - x.Precision * ab[j]);
					}

					Gaussian old_result = (Gaussian)result[i, k].Clone();
					Gaussian rik = result[i, k];
					rik.Precision = prec;
					rik.MeanTimesPrecision = pm;
					result[i, k] = rik;
					Gaussian partial = A[i, k] / old_result;
					Gaussian newPosterior = partial * rik;
					Am = newPosterior.GetMean();
					for (int j = 0; j < cols; j++) ab[j] += Am * B[k, j];
				}
			}
			return result;
		}

		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_B">Previous outgoing message to 'B'.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the exponential of the average log-factor value, where the average is over all arguments except 'B'.
		/// Because the factor is deterministic, 'matrixMultiply' is integrated out before taking the logarithm.
		/// The formula is <c>exp(sum_(A) p(A) log(sum_matrixMultiply p(matrixMultiply) factor(matrixMultiply,A,B)))</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D BAverageLogarithm([SkipIfUniform] GaussianArray2D matrixMultiply, [SkipIfUniform] GaussianArray2D A, [Proper, Stochastic] GaussianArray2D B, GaussianArray2D to_B)
		{
			GaussianArray2D result = to_B;
			if (result == null) result = new GaussianArray2D(B.GetLength(0), B.GetLength(1));
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = A.GetLength(1);
			double[] ab = new double[rows];
			//Gaussian temp = new Gaussian();
			for (int j = 0; j < cols; j++) {
				for (int i = 0; i < rows; i++) {
					double sum = 0.0;
					for (int k = 0; k < inner; k++) {
						sum += A[i, k].GetMean() * B[k, j].GetMean();
					}
					ab[i] = sum;
				}
				for (int k = 0; k < inner; k++) {
					double prec = 0.0;
					double pm = 0.0;
					double Bm = B[k, j].GetMean();
					for (int i = 0; i < rows; i++) {
						double Am, Av;
						A[i, k].GetMeanAndVariance(out Am, out Av);
						Gaussian x = matrixMultiply[i, j];
						prec += (Av + Am * Am) * x.Precision;
						//pm += Am * (x.MeanTimesPrecision - x.Precision * (ab[i] - Am * Bm));
						//Replace previous line with:
						ab[i] -= Am * Bm;
						pm += Am * (x.MeanTimesPrecision - x.Precision * ab[i]);
					}

					Gaussian old_result = (Gaussian)result[k, j].Clone();
					Gaussian rkj = result[k, j];
					rkj.Precision = prec;
					rkj.MeanTimesPrecision = pm;
					result[k, j] = rkj;
					Gaussian partial = B[k, j] / old_result;
					Gaussian newPosterior = partial * rkj;
					Bm = newPosterior.GetMean();
					for (int i = 0; i < rows; i++) ab[i] += Bm * A[i, k].GetMean();
				}
			}
			return result;
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="to_B">Previous outgoing message to 'B'.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' with 'matrixMultiply' integrated out.
		/// The formula is <c>sum_matrixMultiply p(matrixMultiply) factor(matrixMultiply,A,B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D BAverageLogarithm([SkipIfUniform] GaussianArray2D matrixMultiply, double[,] A, [Proper, Stochastic] GaussianArray2D B, GaussianArray2D to_B)
		{
			GaussianArray2D result = to_B;
			if (result == null) result = new GaussianArray2D(B.GetLength(0), B.GetLength(1));
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = A.GetLength(1);
			double[] ab = new double[rows];
			//Gaussian temp = new Gaussian();
			for (int j = 0; j < cols; j++) {
				for (int i = 0; i < rows; i++) {
					double sum = 0.0;
					for (int k = 0; k < inner; k++) {
						sum += A[i, k] * B[k, j].GetMean();
					}
					ab[i] = sum;
				}
				for (int k = 0; k < inner; k++) {
					double prec = 0.0;
					double pm = 0.0;
					double Bm = B[k, j].GetMean();
					for (int i = 0; i < rows; i++) {
						double Am = A[i, k];
						Gaussian x = matrixMultiply[i, j];
						prec += (Am * Am) * x.Precision;
						//pm += Am * (x.MeanTimesPrecision - x.Precision * (ab[i] - Am * Bm));
						//Replace previous line with:
						ab[i] -= Am * Bm;
						pm += Am * (x.MeanTimesPrecision - x.Precision * ab[i]);
					}

					Gaussian old_result = (Gaussian)result[k, j].Clone();
					Gaussian rkj = result[k, j];
					rkj.Precision = prec;
					rkj.MeanTimesPrecision = pm;
					result[k, j] = rkj;
					Gaussian partial = B[k, j] / old_result;
					Gaussian newPosterior = partial * rkj;
					Bm = newPosterior.GetMean();
					for (int i = 0; i < rows; i++) ab[i] += Bm * A[i, k];
				}
			}
			return result;
		}

		const string NotSupportedMessage = "Variational Message Passing does not support a MatrixMultiply factor with fixed output.";

		/// <summary>
		/// VMP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <returns>The outgoing VMP message to the 'A' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(MatrixMultiplyOp.NotSupportedMessage)]
		public static GaussianArray2D AAverageLogarithm(double[,] matrixMultiply)
		{
			throw new NotSupportedException(MatrixMultiplyOp.NotSupportedMessage);
		}
		/// <summary>
		/// VMP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <returns>The outgoing VMP message to the 'B' argument</returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(MatrixMultiplyOp.NotSupportedMessage)]
		public static GaussianArray2D BAverageLogarithm(double[,] matrixMultiply)
		{
			throw new NotSupportedException(MatrixMultiplyOp.NotSupportedMessage);
		}

		// AverageConditional ------------------------------------------------------------------------------------------------------------------------

		const string LowRankNotSupportedMessage = "A MatrixMultiply factor with fixed output is not yet implemented for Expectation Propagation.";
		const string BothRandomNotSupportedMessage = "A MatrixMultiply factor between two Gaussian arrays is not yet implemented for Expectation Propagation.  Try using Variational Message Passing.";

		/// <summary>
		/// EP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[p(matrixMultiply) sum_(A,B) p(A,B) factor(matrixMultiply,A,B)]/p(matrixMultiply)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(MatrixMultiplyOp.BothRandomNotSupportedMessage)]
		public static GaussianArray2D MatrixMultiplyAverageConditional([SkipIfUniform] GaussianArray2D A, [SkipIfUniform] GaussianArray2D B, GaussianArray2D result)
		{
			throw new NotSupportedException(MatrixMultiplyOp.BothRandomNotSupportedMessage);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(matrixMultiply,B) p(matrixMultiply,B) factor(matrixMultiply,A,B)]/p(A)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		[NotSupported(MatrixMultiplyOp.BothRandomNotSupportedMessage)]
		public static GaussianArray2D AAverageConditional(GaussianArray2D matrixMultiply, [SkipIfUniform] GaussianArray2D B, GaussianArray2D result)
		{
			throw new NotSupportedException(MatrixMultiplyOp.BothRandomNotSupportedMessage);
		}
		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(matrixMultiply,A) p(matrixMultiply,A) factor(matrixMultiply,A,B)]/p(B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		[NotSupported(MatrixMultiplyOp.BothRandomNotSupportedMessage)]
		public static GaussianArray2D BAverageConditional(GaussianArray2D matrixMultiply, [SkipIfUniform] GaussianArray2D A, GaussianArray2D result)
		{
			throw new NotSupportedException(MatrixMultiplyOp.BothRandomNotSupportedMessage);
		}

		public static double LogAverageFactor(double[,] matrixMultiply, double[,] A, double[,] B)
		{
			for (int k = 0; k < A.GetLength(0); k++) {
				for (int k2 = 0; k2 < B.GetLength(1); k2++) {
					double sum = 0;
					for (int i = 0; i < A.GetLength(1); i++) {
						sum += A[k, i] * B[i, k2];
					}
					if (matrixMultiply[k, k2] != sum) return Double.NegativeInfinity;
				}
			}
			return 0.0;
		}
		public static double LogEvidenceRatio(double[,] matrixMultiply, double[,] A, double[,] B) { return LogAverageFactor(matrixMultiply, A, B); }
		public static double AverageLogFactor(double[,] matrixMultiply, double[,] A, double[,] B) { return LogAverageFactor(matrixMultiply, A, B); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(matrixMultiply,B) p(matrixMultiply,B) factor(matrixMultiply,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(GaussianArray2D matrixMultiply, double[,] A, GaussianArray2D B)
		{
			GaussianArray2D to_matrixMultiply = MatrixMultiplyAverageConditional(A, B, null);
			return to_matrixMultiply.GetLogAverageOf(matrixMultiply);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(matrixMultiply,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double[,] matrixMultiply, double[,] A, GaussianArray2D B)
		{
			GaussianArray2D to_matrixMultiply = MatrixMultiplyAverageConditional(A, B, null);
			return to_matrixMultiply.GetLogProb(matrixMultiply);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(matrixMultiply,A) p(matrixMultiply,A) factor(matrixMultiply,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(GaussianArray2D matrixMultiply, GaussianArray2D A, double[,] B)
		{
			GaussianArray2D to_matrixMultiply = MatrixMultiplyAverageConditional(A, B, null);
			return to_matrixMultiply.GetLogAverageOf(matrixMultiply);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(matrixMultiply,A,B))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double[,] matrixMultiply, GaussianArray2D A, double[,] B)
		{
			GaussianArray2D to_matrixMultiply = MatrixMultiplyAverageConditional(A, B, null);
			return to_matrixMultiply.GetLogProb(matrixMultiply);
		}

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(matrixMultiply,B) p(matrixMultiply,B) factor(matrixMultiply,A,B) / sum_matrixMultiply p(matrixMultiply) messageTo(matrixMultiply))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(GaussianArray2D matrixMultiply, double[,] A, GaussianArray2D B) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(matrixMultiply,A) p(matrixMultiply,A) factor(matrixMultiply,A,B) / sum_matrixMultiply p(matrixMultiply) messageTo(matrixMultiply))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(GaussianArray2D matrixMultiply, GaussianArray2D A, double[,] B) { return 0.0; }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(B) p(B) factor(matrixMultiply,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double[,] matrixMultiply, double[,] A, GaussianArray2D B) { return LogAverageFactor(matrixMultiply, A, B); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="A">Incoming message from 'A'.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(A) p(A) factor(matrixMultiply,A,B))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double[,] matrixMultiply, GaussianArray2D A, double[,] B) { return LogAverageFactor(matrixMultiply, A, B); }

		/// <summary>
		/// EP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[p(matrixMultiply) sum_(B) p(B) factor(matrixMultiply,A,B)]/p(matrixMultiply)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D MatrixMultiplyAverageConditional(double[,] A, [SkipIfUniform] GaussianArray2D B, GaussianArray2D result)
		{
			if (result == null) result = new GaussianArray2D(A.GetLength(0), B.GetLength(1));
			// x[i,j] = sum_k a[i,k]*b[k,j]
			// E(x[i,j]) = sum_k a[i,k]*E(b[k,j])
			// var(x[i,j]) = sum_k a[i,k]^2*var(b[k,j])
			int rows = result.GetLength(0);
			int cols = result.GetLength(1);
			int inner = A.GetLength(1);
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < cols; j++) {
					double mean = 0; double var = 0;
					for (int k = 0; k < inner; k++) {
						double Am = A[i, k];
						double Bm, Bv;
						B[k, j].GetMeanAndVariance(out Bm, out Bv);
						mean += Am * Bm;
						var += Bv * Am * Am;
					}
					Gaussian rij = result[i, j];
					rij.SetMeanAndVariance(mean, var);
					result[i, j] = rij;
				}
			}
			return result;
		}
		/// <summary>
		/// EP message to 'matrixMultiply'
		/// </summary>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'matrixMultiply' as the random arguments are varied.
		/// The formula is <c>proj[p(matrixMultiply) sum_(A) p(A) factor(matrixMultiply,A,B)]/p(matrixMultiply)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static GaussianArray2D MatrixMultiplyAverageConditional([SkipIfUniform] GaussianArray2D A, double[,] B, GaussianArray2D result)
		{
			return MatrixMultiplyAverageConditional(B, A, result);
		}
		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Incoming message from 'A'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'A' as the random arguments are varied.
		/// The formula is <c>proj[p(A) sum_(matrixMultiply) p(matrixMultiply) factor(matrixMultiply,A,B)]/p(A)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
		public static GaussianArray2D AAverageConditional([SkipIfUniform] GaussianArray2D matrixMultiply, [SkipIfUniform] GaussianArray2D A, double[,] B, GaussianArray2D result)
		{
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = B.GetLength(0);
			if (result == null) result = new GaussianArray2D(rows, inner);
			// sum_{i,j} (m[i,j] - a[i,:]*b[:,j])^2/v[i,j] = 
			// sum_{i,j} (m[i,j]^2 - 2m[i,j]a[i,:]*b[:,j] + a[i,:]*(b[:,j] b[:,j]')*a[i,:]')/v[i,j]
			// meanTimesPrec(a[i,:]) = sum_j (m[i,j]/v[i,j]) b[:,j]
			// prec(a[i,:]) = sum_j b[:,j]*b[:,j]'/v[i,j]
			Vector bj = Vector.Zero(inner);
			Vector mean = Vector.Zero(inner);
			VectorGaussian ai = new VectorGaussian(inner);
			PositiveDefiniteMatrix variance = new PositiveDefiniteMatrix(inner, inner);
			for (int i = 0; i < rows; i++) {
				ai.Precision.SetAllElementsTo(0.0);
				ai.MeanTimesPrecision.SetAllElementsTo(0.0);
				// we are projecting from family of full covariance Gaussians to diagonal
				// covariance, so we should include the context
				for (int c = 0; c < inner; c++) {
					ai.Precision[c, c] = A[i, c].Precision;
					ai.MeanTimesPrecision[c] = A[i, c].MeanTimesPrecision;
				}
				for (int j = 0; j < cols; j++) {
					Gaussian xij = matrixMultiply[i, j];
					for (int k = 0; k < inner; k++) {
						bj[k] = B[k, j];
					}
					if (xij.IsPointMass) throw new NotImplementedException(LowRankNotSupportedMessage);
					ai.Precision.SetToSumWithOuter(ai.Precision, xij.Precision, bj, bj);
					ai.MeanTimesPrecision.SetToSum(1.0, ai.MeanTimesPrecision, xij.MeanTimesPrecision, bj);
				}
				ai.GetMeanAndVariance(mean, variance);
				for (int k = 0; k < inner; k++) {
					Gaussian rik = result[i, k];
					rik.SetMeanAndVariance(mean[k], variance[k, k]);
					result[i, k] = rik / A[i, k];
				}
			}
			return result;
		}

		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Incoming message from 'matrixMultiply'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="B">Incoming message from 'B'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'B' as the random arguments are varied.
		/// The formula is <c>proj[p(B) sum_(matrixMultiply) p(matrixMultiply) factor(matrixMultiply,A,B)]/p(B)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="matrixMultiply"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
		public static GaussianArray2D BAverageConditional([SkipIfUniform] GaussianArray2D matrixMultiply, double[,] A, [SkipIfUniform] GaussianArray2D B, GaussianArray2D result)
		{
			int rows = matrixMultiply.GetLength(0);
			int cols = matrixMultiply.GetLength(1);
			int inner = A.GetLength(1);
			if (result == null) result = new GaussianArray2D(inner, cols);
			var ai = DenseVector.Zero(inner);
			var mean = DenseVector.Zero(inner);
			PositiveDefiniteMatrix variance = new
			PositiveDefiniteMatrix(inner, inner);
			var bj = new VectorGaussian(inner);
			for (int j = 0; j < cols; j++) {
				bj.Precision.SetAllElementsTo(0);
				bj.MeanTimesPrecision.SetAllElementsTo(0);
				// we are projecting from family of full covariance Gaussians to diagonal
				// covariance, so we should include the context
				for (int c = 0; c < inner; c++) {
					bj.Precision[c, c] = B[c, j].Precision;
					bj.MeanTimesPrecision[c] = B[c, j].MeanTimesPrecision;
				}
				for (int i = 0; i < rows; i++) {
					Gaussian xij = matrixMultiply[i, j];
					for (int k = 0; k < inner; k++) {
						ai[k] = A[i, k];
					}
					if (xij.IsPointMass) throw new NotImplementedException(LowRankNotSupportedMessage);
					bj.Precision.SetToSumWithOuter(bj.Precision, xij.Precision, ai, ai);
					bj.MeanTimesPrecision.SetToSum(1.0, bj.MeanTimesPrecision, xij.MeanTimesPrecision, ai);
				}
				bj.GetMeanAndVariance(mean, variance);
				for (int k = 0; k < inner; k++) {
					Gaussian rkj = result[k, j];
					rkj.SetMeanAndVariance(mean[k], variance[k, k]);
					result[k, j] = rkj / B[k, j];
				}
			}
			return result;
		}

		/// <summary>
		/// EP message to 'A'
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="B">Constant value for 'B'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'A' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(MatrixMultiplyOp.LowRankNotSupportedMessage)]
		public static GaussianArray2D AAverageConditional(double[,] matrixMultiply, double[,] B, GaussianArray2D result)
		{
			throw new NotImplementedException(MatrixMultiplyOp.LowRankNotSupportedMessage);
		}

		/// <summary>
		/// EP message to 'B'
		/// </summary>
		/// <param name="matrixMultiply">Constant value for 'matrixMultiply'.</param>
		/// <param name="A">Constant value for 'A'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'B' conditioned on the given values.
		/// </para></remarks>
		[NotSupported(MatrixMultiplyOp.LowRankNotSupportedMessage)]
		public static GaussianArray2D BAverageConditional(double[,] matrixMultiply, double[,] A, GaussianArray2D result)
		{
			throw new NotImplementedException(MatrixMultiplyOp.LowRankNotSupportedMessage);
		}
	}
}
