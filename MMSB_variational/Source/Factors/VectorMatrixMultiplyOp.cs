// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;
using GaussianArray2D = MicrosoftResearch.Infer.Distributions.DistributionStructArray2D<MicrosoftResearch.Infer.Distributions.Gaussian, double>;
using VectorGaussianArray = MicrosoftResearch.Infer.Distributions.DistributionRefArray<MicrosoftResearch.Infer.Distributions.VectorGaussian, MicrosoftResearch.Infer.Maths.Vector>;

namespace MicrosoftResearch.Infer.Factors
{
    /// <summary>
    /// Provides outgoing messages for <see cref="Factor.VectorMatrixMultiply"/>, given random arguments to the function.
    /// </summary>
    [FactorMethod(typeof(Factor), "VectorMatrixMultiply")]
    [Quality(QualityBand.Preview)]
    public static class VectorMatrixMultiplyOp
    {
        //-- VMP -------------------------------------------------------------------------------------------------------------
        /// <summary>
        /// Evidence message for VMP.
        /// </summary>
        /// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
        /// <remarks><para>
        /// The formula for the result is <c>int log(f(x)) q(x) dx</c>
        /// where <c>x = (matrixMultiply,A,B)</c>.
        /// </para></remarks>
        [Skip]
        public static double AverageLogFactor() { return 0.0; }

        /// <summary>
        /// VMP message to 'A'.
        /// </summary>
        /// <param name="X">Incoming message from 'X'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
        /// <param name="A">Incoming message from 'A'.</param>
        /// <param name="result">Modified to contain the outgoing message.</param>
        /// <returns><paramref name="result"/></returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'A'.
        /// The formula is <c>int log(f(A,x)) q(x) dx</c> where <c>x = (X,B)</c>.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="X"/> is not a proper distribution</exception>
        /// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
        public static VectorGaussianArray AAverageLogarithm([SkipIfAllUniform] GaussianArray2D X, [SkipIfAllUniform] VectorGaussianArray B, [SkipIfAllUniform] VectorGaussianArray result)
        {
            int I = X.GetLength(0), J = X.GetLength(1), K = B[0].Dimension;
            var Ebbp = new PositiveDefiniteMatrix[J];
            var mb = new Vector[J];
            for (int j = 0; j < J; j++)
            {
                Ebbp[j] = new PositiveDefiniteMatrix(K, K);
                mb[j] = Vector.Zero(K);
                B[j].GetMeanAndVariance(mb[j], Ebbp[j]);
                Ebbp[j].SetToSumWithOuter(Ebbp[j], 1, mb[j], mb[j]);
            }
            for (int i = 0; i < I; i++)
            {
                result[i].Precision.SetAllElementsTo(0);
                result[i].MeanTimesPrecision.SetAllElementsTo(0);
                for (int j = 0; j < J; j++)
                {
                    // nb: would be more memory efficient to have a SetToAPlusCB routine
                    result[i].Precision.SetToSum(result[i].Precision, Ebbp[j] * X[i, j].Precision);
                    result[i].MeanTimesPrecision.SetToSum(result[i].MeanTimesPrecision, mb[j] * X[i, j].MeanTimesPrecision);
                }
            }
            return result;
        }

        /// <summary>
        /// VMP message to 'B'.
        /// </summary>
        /// <param name="X">Incoming message from 'X'. Must be a proper distribution.  If any element is uniform, the result will be uniform.</param>
        /// <param name="A">Incoming message from 'A'.</param>
        /// <param name="result">Modified to contain the outgoing message.</param>
        /// <returns><paramref name="result"/></returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'A'.
        /// The formula is <c>int log(f(A,x)) q(x) dx</c> where <c>x = (X,B)</c>.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="X"/> is not a proper distribution</exception>
        /// <exception cref="ImproperMessageException"><paramref name="B"/> is not a proper distribution</exception>
        public static VectorGaussianArray BAverageLogarithm([SkipIfAllUniform] GaussianArray2D X, [SkipIfAllUniform] VectorGaussianArray A, VectorGaussianArray result)
        {
            int I = X.GetLength(0), J = X.GetLength(1), K = A[0].Dimension;
            var Eaap = new PositiveDefiniteMatrix[I];
            var ma = new Vector[I];
            for (int i = 0; i < I; i++)
            {
                Eaap[i] = new PositiveDefiniteMatrix(K, K);
                ma[i] = Vector.Zero(K);
                A[i].GetMeanAndVariance(ma[i], Eaap[i]);
                Eaap[i].SetToSumWithOuter(Eaap[i], 1, ma[i], ma[i]);
            }
            for (int j = 0; j < J; j++)
            {
                result[j].Precision.SetAllElementsTo(0);
                result[j].MeanTimesPrecision.SetAllElementsTo(0);
                for (int i = 0; i < I; i++)
                {
                    // nb: would be more memory efficient to have a SetToAPlusCB routine
                    result[j].Precision.SetToSum(result[j].Precision, Eaap[i] * X[i, j].Precision);
                    result[j].MeanTimesPrecision.SetToSum(result[j].MeanTimesPrecision, ma[i] * X[i, j].MeanTimesPrecision);
                }
            }
            return result;
        }

        /// <summary>
        /// VMP message to 'X'.
        /// </summary>
        /// <param name="A">Incoming message from 'a'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
        /// <param name="B">Incoming message from 'b'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
        /// <returns>The outgoing VMP message to the 'X' argument.</returns>
        /// <remarks><para>
        /// The outgoing message is the exponential of the integral of the log-factor times incoming messages, over all arguments except 'innerProduct'.
        /// The formula is <c>int log(f(innerProduct,x)) q(x) dx</c> where <c>x = (a,b)</c>.
        /// </para></remarks>
        /// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
        /// <exception cref="ImproperMessageException"><paramref name="A"/> is not a proper distribution</exception>
        public static GaussianArray2D XAverageLogarithm([SkipIfAllUniform] VectorGaussianArray A, [SkipIfAllUniform] VectorGaussianArray B, GaussianArray2D result)
        {
            int I = A.GetLength(0), J = B.Count, K = A[0].Dimension;
            if (result == null) result = new DistributionStructArray2D<Gaussian, double>(I, J);
            // p(x|a,b) = N(E[a]'*E[b], E[b]'*var(a)*E[b] + E[a]'*var(b)*E[a] + trace(var(a)*var(b)))
            var ma = new Vector[I];
            var mb = new Vector[J];
            var va = new PositiveDefiniteMatrix[I];
            var vb = new PositiveDefiniteMatrix[J];
            for (int i = 0; i < I; i++)
            {
                ma[i] = Vector.Zero(K);
                va[i] = new PositiveDefiniteMatrix(K, K);
                A[i].GetMeanAndVariance(ma[i], va[i]);
            }
            for (int j = 0; j < J; j++)
            {
                mb[j] = Vector.Zero(K);
                vb[j] = new PositiveDefiniteMatrix(K, K);
                B[j].GetMeanAndVariance(mb[j], vb[j]);
            }
            // Uses John Winn's rule for deterministic factors.
            // Strict variational inference would set the variance to 0.
            for (int i = 0; i < I; i++)
                for (int j = 0; j < J; j++)
                {
                    var rij = result[i, j];
                    rij.SetMeanAndVariance(ma[i].Inner(mb[j]), va[i].QuadraticForm(mb[j]) + vb[j].QuadraticForm(ma[i]) + va[i].Inner(vb[j]));
                    result[i, j] = rij;
                }
            return result;
        }

    }
}
