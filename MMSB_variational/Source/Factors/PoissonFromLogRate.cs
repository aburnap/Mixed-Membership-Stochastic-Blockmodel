// (C) Copyright 2008 Microsoft Research Cambridge
using System;
using System.Collections.Generic;
using System.Text;
using MicrosoftResearch.Infer.Distributions;
using MicrosoftResearch.Infer.Maths;

namespace MicrosoftResearch.Infer.Factors
{
	/// <summary>
	/// Provides outgoing messages for <see cref="Factor.PoissonFromLogRate"/>, given random arguments to the function.
	/// Performs KL minimisation using gradient matching, a distributed gradient descent algorithm. 
	/// </summary>
    [FactorMethod(typeof(Factor), "PoissonFromLogRate")]
    public static class PoissonFromLogRateOp
	{
        public static double damping = 0.0; 
        internal const string NotSupportedMessage = "PoissonFromLogRate not yet implemented for Expectation Propagation.  Try Variational Message Passing.";

		/// <summary>
		/// Evidence message for VMP.
		/// </summary>
		/// <param name="sample">Constant value for 'sample'.</param>
		/// <param name="logOdds">Incoming message from 'logOdds'.</param>
		/// <returns><c>sum_x marginal(x)*log(factor(x))</c></returns>
		/// <remarks><para>
		/// The formula for the result is <c>int log(f(x)) q(x) dx</c>
		/// </para></remarks>
		public static double AverageLogFactor(int sample, [Proper, SkipIfUniform] Gaussian logRate)
		{
            double m,v;
            logRate.GetMeanAndVariance(out m,out v);
            return m * (double)sample - Math.Exp(m + v / 2) - MMath.GammaLn(sample + 1); 
		}

		
		/// <summary>
		/// Gradient matching VMP message from factor to logOdds variable
		/// </summary>
		/// <param name="sample">Incoming message from 'sample'.</param>
		/// <param name="logOdds">Incoming message. Must be a proper distribution.  If uniform, the result will be uniform.</param>
		/// <param name="result">Previous message sent, used for damping</param>
		/// <returns>The outgoing VMP message.</returns>
		/// <remarks><para>
		/// The outgoing message is the Gaussian approximation to the factor which results in the 
		/// same derivatives of the KL(q||p) divergence with respect to the parameters of the posterior
		/// as if the true factor had been used.
		/// </para></remarks>
		/// <exception cref="ImproperMessageException"><paramref name="logOdds"/> is not a proper distribution</exception>
        public static Gaussian LogRateAverageLogarithm(int sample, [Proper, SkipIfUniform] Gaussian logRate, Gaussian result)
		{
            double m,v;
            logRate.GetMeanAndVariance(out m,out v);
            Gaussian message = new Gaussian();
            message.Precision = Math.Exp(m + v / 2);
            message.MeanTimesPrecision = (m - 1) * Math.Exp(m + v / 2) + sample;
            if (damping == 0)
                return message;
            else
                return (message ^ (1 - damping)) * (result ^ damping); 
		}
	}
}
