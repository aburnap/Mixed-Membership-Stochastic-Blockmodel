// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.FunctionEvaluate"/>, given random arguments to the function.
	/// </summary>
	[FactorMethod(typeof(Factor), "FunctionEvaluate", typeof(IFunction), typeof(Vector))]
	[Quality(QualityBand.Preview)]
	public static class SparseGPOp
	{
		public static double LogAverageFactor(double y, IFunction func, Vector x)
		{
			return (y == func.Evaluate(x)) ? 0.0 : Double.NegativeInfinity;
		}
		public static double LogEvidenceRatio(double y, IFunction func, Vector x) { return LogAverageFactor(y, func, x); }
		public static double AverageLogFactor(double y, IFunction func, Vector x) { return LogAverageFactor(y, func, x); }

		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="func">Incoming message from 'func'.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <returns>Logarithm of the factor's average value across the given argument distributions</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(func) p(func) factor(y,func,x))</c>.
		/// </para></remarks>
		public static double LogAverageFactor(double y, SparseGP func, Vector x)
		{
			Gaussian to_y = YAverageConditional(func, x);
			return to_y.GetLogProb(y);
		}
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="func">Incoming message from 'func'.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(func) p(func) factor(y,func,x))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		public static double LogEvidenceRatio(double y, SparseGP func, Vector x) { return LogAverageFactor(y, func, x); }
		/// <summary>
		/// Evidence message for EP
		/// </summary>
		/// <param name="y">Incoming message from 'y'.</param>
		/// <returns>Logarithm of the factor's contribution the EP model evidence</returns>
		/// <remarks><para>
		/// The formula for the result is <c>log(sum_(y) p(y) factor(y,func,x) / sum_y p(y) messageTo(y))</c>.
		/// Adding up these values across all factors and variables gives the log-evidence estimate for EP.
		/// </para></remarks>
		[Skip]
		public static double LogEvidenceRatio(Gaussian y) { return 0.0; }

		/// <summary>
		/// EP message to 'y'
		/// </summary>
		/// <param name="func">Incoming message from 'func'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <returns>The outgoing EP message to the 'y' argument</returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'y' as the random arguments are varied.
		/// The formula is <c>proj[p(y) sum_(func) p(func) factor(y,func,x)]/p(y)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="func"/> is not a proper distribution</exception>
		public static Gaussian YAverageConditional([SkipIfUniform] SparseGP func, Vector x)
		{
			return func.Marginal(x);
		}

		/// <summary>
		/// EP message to 'func'
		/// </summary>
		/// <param name="y">Incoming message from 'y'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="func">Incoming message from 'func'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is a distribution matching the moments of 'func' as the random arguments are varied.
		/// The formula is <c>proj[p(func) sum_(y) p(y) factor(y,func,x)]/p(func)</c>.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="y"/> is not a proper distribution</exception>
		/// <exception cref="ImproperMessageException"><paramref name="func"/> is not a proper distribution</exception>
		public static SparseGP FuncAverageConditional([SkipIfUniform] Gaussian y, [SkipIfUniform] SparseGP func, Vector x, SparseGP result)
		{
			if (y.IsUniform() || func.IsUniform()) { result.SetToUniform(); return result; }
			result.FixedParameters = func.FixedParameters;
			result.IncludePrior = false;
			
			double vf = func.Variance(x);
			double my, vy;
			y.GetMeanAndVariance(out my, out vy);
			Vector kbx = func.FixedParameters.KernelOf_X_B(x);
			Vector proj = func.FixedParameters.InvKernelOf_B_B * kbx;
			double prec = 1.0/(vy + vf - func.Var_B_B.QuadraticForm(proj));
			result.InducingDist.Precision.SetToOuter(proj, proj);
			result.InducingDist.Precision.Scale(prec);
			result.InducingDist.MeanTimesPrecision.SetTo(proj);
			result.InducingDist.MeanTimesPrecision.Scale(prec*my);
			result.ClearCachedValues();
			return result;
		}
		/// <summary>
		/// EP message to 'func'
		/// </summary>
		/// <param name="y">Constant value for 'y'.</param>
		/// <param name="func">Incoming message from 'func'. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="x">Constant value for 'x'.</param>
		/// <param name="result">Modified to contain the outgoing message</param>
		/// <returns><paramref name="result"/></returns>
		/// <remarks><para>
		/// The outgoing message is the factor viewed as a function of 'func' conditioned on the given values.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="func"/> is not a proper distribution</exception>
		public static SparseGP FuncAverageConditional(double y, [SkipIfUniform] SparseGP func, Vector x, SparseGP result)
		{
			return FuncAverageConditional(Gaussian.PointMass(y), func, x, result);
		}

#if false
        /// <summary>
        /// EP message to 'x'.
        /// </summary>
        /// <param name="y">Incoming message from 'y'.</param>
        /// <param name="func">Incoming message from 'func'.</param>
        /// <param name="x">Constant value for 'x'.</param>
        /// <param name="result">Modified to contain the outgoing message.</param>
        /// <returns><paramref name="result"/></returns>
        /// <remarks><para>
        /// The outgoing message is the integral of the factor times incoming messages, over all arguments except 'x'.
        /// The formula is <c>int f(x,x) q(x) dx</c> where <c>x = (y,func)</c>.
        /// </para></remarks>
        public static VectorGaussian XAverageConditional(Gaussian y, SparseGP func, Vector x, VectorGaussian result)
        {
            // Doesn't matter what message we return as we are only supporting a point x
            result.SetToUniform();
            return result;
        }
#endif
	}
}
